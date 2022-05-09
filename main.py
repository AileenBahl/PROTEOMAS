# -*- coding: utf-8 -*-
"""
Created on Mon Dec 06 14:13:14 2021

@author: Aileen Bahl
"""

from readData import readData
from preprocessData import preprocessData
from downstreamAnalysis import downstreamAnalysis
from generatePlots import generatePlots
import os
import pandas as pd
import numpy as np

myBasepath = r'C:/Users/1/Desktop/PROTEOMAS'

np.random.seed(1)


def main():

    ### Create output folder
    if not os.path.exists(os.path.join(myBasepath, 'Output')):
        os.makedirs(os.path.join(myBasepath, 'Output'))

    ### Get all filenames
        myFilenames = readData.getAllFilenames(myBasepath)

    ### Extract project ids
    myProjectIds = readData.extractProjectIds(myFilenames)

    ### Continue analyzing one project after the other
    for myId in myProjectIds:

        ### Create one output folder for each project
        if not os.path.exists(os.path.join(myBasepath, 'Output', myId)):
            os.makedirs(os.path.join(myBasepath, 'Output', myId))

        ### Create a folder for all plots
        if not os.path.exists(os.path.join(myBasepath, 'Output', myId, 'Plots')):
            os.makedirs(os.path.join(myBasepath, 'Output', myId, 'Plots'))

        ### Grep all files belonging to the current project
        myFiles = [i for i in myFilenames if myId in i]

        ### Read input file (proteinGroups.txt)
        myProteinGroupsFilename = [i for i in myFiles if 'proteinGroups' in i]
        myProteinGroupsFiles = readData.readFiles(myProteinGroupsFilename[0],'proteinGroups','.txt')

        # read ConditionAssignment file if present
        myConditionAssignmentFilename = [i for i in myFiles if 'ConditionAssignment' in i]
        if len(myConditionAssignmentFilename) != 0:
            myConditionAssignmentFiles = readData.readFiles(myConditionAssignmentFilename[0],'ConditionAssignment','.csv')
        else:
            myConditionAssignmentFiles = None

        # read Comparisons file if present
        myComparisonsFilename = [i for i in myFiles if 'Comparisons' in i]
        if len(myComparisonsFilename) != 0:
            myComparisonsFiles = readData.readFiles(myComparisonsFilename[0], 'Comparisons', '.csv')
        else:
            myComparisonsFiles = None

        # read assignedK file if present
        myAssignedKFilename = [i for i in myFiles if 'AssignedK' in i]
        if len(myAssignedKFilename) != 0:
            myAssignedKFiles = readData.readFiles(myAssignedKFilename[0], 'AssignedK', '.txt')
        else:
            myAssignedKFiles = None


        ### filter proteinGroups.txt

        ### filter reverse
        myProteinGroupsFiles_reverseFiltered = preprocessData.filterUnnecessaryProteins(myProteinGroupsFiles,'Reverse')

        ### filter contaminants
        myProteinGroupsFiles_contaminantFiltered = preprocessData.filterUnnecessaryProteins(myProteinGroupsFiles_reverseFiltered, 'Contaminant')

        ### filter proteins only identified by site
        myProteinGroupsFiles_onlyIdentifiedBySiteFiltered = preprocessData.filterUnnecessaryProteins(myProteinGroupsFiles_contaminantFiltered, 'Only identified by site')
        myProteinGroupsFiles_onlyIdentifiedBySiteFiltered.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_FilteredUnnecessaryCols_'+myId+'.xlsx'))

        ### filter proteins identified by <2 peptides and <1 unique peptide
        myProteinGroupsFiles_uncertainProteinsFiltered = preprocessData.filterUncertainProteins(myProteinGroupsFiles_onlyIdentifiedBySiteFiltered)
        myProteinGroupsFiles_uncertainProteinsFiltered.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_FilteredUncertainCols_' + myId + '.xlsx'))

        ### log2 transform data
        myProteinGroupsFiles_log2 = preprocessData.log2Transform(myProteinGroupsFiles_uncertainProteinsFiltered)
        myProteinGroupsFiles_log2.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_log2_' + myId + '.xlsx'))


        ### Assign groups
        ### if ConditionAssignmentFile is given, use defined groups
        if(myConditionAssignmentFiles is not None):
            myProteinGroupsFiles_groupsAssigned = preprocessData.assignGroups(myProteinGroupsFiles_log2, myConditionAssignmentFiles)
            myProteinGroupsFiles_groupsAssigned.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_groupsAssigned_' + myId + '.xlsx'))
        ### if ConditionAssignmentFile is not given, assign group via k-means clustering
        else:
            myProteinGroupsFiles_log2_copy = myProteinGroupsFiles_log2.copy()
            ### Impute values from normal distribution as NbClust requires a complete dataset and RF imputation within each group cannot be done as no groups are defined yet
            myProteinGroupsFiles_imputedComplete = preprocessData.imputeFromNormal(myProteinGroupsFiles_log2_copy)
            ### if AssignedKFile is provided, use given k
            if (myAssignedKFiles is not None):
                myProteinGroupsFiles_groupsAssigned = preprocessData.predictGroupAssignmentFromGivenK(myProteinGroupsFiles_log2_copy, myProteinGroupsFiles_log2, myAssignedKFiles)
                myProteinGroupsFiles_groupsAssigned.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_groupsAssignedFromGivenK_' + myId + '.xlsx'))
            ### if AssignedKFile is not provieded, predict k using NbClust
            else:
                myProteinGroupsFiles_groupsAssigned = preprocessData.predictGroupAssignment(myProteinGroupsFiles_log2_copy, myProteinGroupsFiles_log2)
                myProteinGroupsFiles_groupsAssigned.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_groupsAssignedPredicted_' + myId + '.xlsx'))


        ### Keep only treatments with at least triplicate measurements (as this will be a pre-requisit for the t-test and meaningful outlier detection)
        myProteinGroupsFiles_triplicates = preprocessData.keepOnlyTriplicatesOrMore(myProteinGroupsFiles_groupsAssigned)
        myProteinGroupsFiles_triplicates.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_triplicates_' + myId + '.xlsx'))

        ### Keep only proteins with at least 70% valid values in one group (to avoid too many false-positives)
        myProteinGroupsFiles_filteredValidValues = preprocessData.filterValid(myProteinGroupsFiles_triplicates)
        myProteinGroupsFiles_filteredValidValues.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_filteredValidValues_' + myId + '.xlsx'))

        #generatePlots.plotHistogram(myProteinGroupsFiles_filteredValidValues, myId, 'filterValidValues', myBasepath)
        #generatePlots.plotBoxplot(myProteinGroupsFiles_filteredValidValues, myId, 'filterValidValues', myBasepath)

        myProteinGroupsFiles_filteredValidValues_copy = myProteinGroupsFiles_filteredValidValues.copy()


        ### impute values from Random Forest imputation (in case there are more than 70% valid values for that protein with this condition) and from downshifted, shrunk normal distribution (otherwise)
        myProteinGroupsFiles_imputed = preprocessData.imputeRFandNormal(myProteinGroupsFiles_filteredValidValues, myId)
        myProteinGroupsFiles_imputed.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_imputedRFandNormal_' + myId + '.xlsx'))

        #generatePlots.plotHistogram(myProteinGroupsFiles_imputed, myId, 'imputed', myBasepath)
        #generatePlots.plotBoxplot(myProteinGroupsFiles_imputed, myId, 'imputed', myBasepath)


        ### Detect outliers using PCAGrid
        ### Remove outliers as long as there are at least still triplicate measurements left
        myProteinGroupsFiles_outliersRemoved = preprocessData.detectOutliers(myProteinGroupsFiles_imputed)
        myProteinGroupsFiles_outliersRemoved.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_outliersRemoved_' + myId + '.xlsx'))


        ### Perform t-test
        if myProteinGroupsFiles_outliersRemoved.shape[0] > 0:
           myProteinGroupsFiles_tTest = downstreamAnalysis.performTTest(myProteinGroupsFiles_outliersRemoved, myProteinGroupsFiles_filteredValidValues_copy, myComparisonsFiles)
            for key, value in myProteinGroupsFiles_tTest.items():
                value.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_tTestResults_'+myId+'_'+key+'.xlsx'))

                #generatePlots.plotVolcano(value, myId, key, myBasepath)

            ### Extract significant proteins (FDR <= 0.05 and log2ratio >= log2(1.5))
            myProteinGroupsFiles_significant = downstreamAnalysis.extractSignificant(myProteinGroupsFiles_tTest)
            for key, value in myProteinGroupsFiles_significant.items():
                value.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_significantProteins_'+myId+'_'+key+'.xlsx'))


            ### Perform enrichment analysis by querying STRING
            myProteinGroupsFiles_STRING = downstreamAnalysis.querySTRING(myProteinGroupsFiles_significant)
            for key, value in myProteinGroupsFiles_STRING.items():
                value.to_excel(os.path.join(myBasepath, 'Output', myId, 'myProjects_STRINGResults_'+myId+'_'+key+'.xlsx'))


# run main

if __name__ == "__main__":
    main()