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
import csv


myBasepath = r'C:/Users/1/Desktop/PROTEOMAS'
myInputFile = 'Input_test'
myOutputFile = 'Output_test'


dbs = ['H']#,'kegg','reactome','wikiPathways','go']


np.random.seed(1)


def main():
    ### Create output folder
    if not os.path.exists(os.path.join(myBasepath, myOutputFile)):
        os.makedirs(os.path.join(myBasepath, myOutputFile))

    ### Get all filenames
    myFilenames = readData.getAllFilenames(myBasepath, myInputFile)

    ### Extract project ids
    myProjectIds = readData.extractProjectIds(myFilenames)

    ### Continue analyzing one project after the other
    for myId in myProjectIds:

        ### Create one output folder for each project
        if not os.path.exists(os.path.join(myBasepath, myOutputFile, myId)):
            os.makedirs(os.path.join(myBasepath, myOutputFile, myId))

        ### Create a folder for all plots
        if not os.path.exists(os.path.join(myBasepath, myOutputFile, myId, 'Plots')):
            os.makedirs(os.path.join(myBasepath, myOutputFile, myId, 'Plots'))

        ### Grep all files belonging to the current project
        myFiles = [i for i in myFilenames if myId in i]

        ### Read input file (proteinGroups.txt)
        myProteinGroupsFilename = [i for i in myFiles if 'proteinGroups' in i]
        myProteinGroupsFiles = readData.readFiles(myProteinGroupsFilename[0], 'proteinGroups', '.txt')

        # read ConditionAssignment file if present
        myConditionAssignmentFilename = [i for i in myFiles if 'ConditionAssignment' in i]
        if len(myConditionAssignmentFilename) != 0:
            myConditionAssignmentFiles = readData.readFiles(myConditionAssignmentFilename[0], 'ConditionAssignment',
                                                            '.xlsx')
        else:
            myConditionAssignmentFiles = None

        # read Comparisons file if present
        myComparisonsFilename = [i for i in myFiles if 'Comparisons' in i]
        if len(myComparisonsFilename) != 0:
            myComparisonsFiles = readData.readFiles(myComparisonsFilename[0], 'Comparisons', '.xlsx')
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
        myProteinGroupsFiles_reverseFiltered = preprocessData.filterUnnecessaryProteins(myProteinGroupsFiles, 'Reverse')

        ### filter contaminants
        myProteinGroupsFiles_contaminantFiltered = preprocessData.filterUnnecessaryProteins(
            myProteinGroupsFiles_reverseFiltered, 'Contaminant')

        ### filter proteins only identified by site
        myProteinGroupsFiles_onlyIdentifiedBySiteFiltered = preprocessData.filterUnnecessaryProteins(
            myProteinGroupsFiles_contaminantFiltered, 'Only identified by site')
        myProteinGroupsFiles_onlyIdentifiedBySiteFiltered.to_csv(
            os.path.join(myBasepath, myOutputFile, myId, 'myProjects_FilteredUnnecessaryCols_' + myId + '.csv'),
            index=False)

        ### filter proteins identified by <2 peptides and <1 unique peptide
        myProteinGroupsFiles_uncertainProteinsFiltered = preprocessData.filterUncertainProteins(
            myProteinGroupsFiles_onlyIdentifiedBySiteFiltered)
        myProteinGroupsFiles_uncertainProteinsFiltered.to_csv(
            os.path.join(myBasepath, myOutputFile, myId, 'myProjects_FilteredUncertainCols_' + myId + '.csv'),
            index=False)

        if myProteinGroupsFiles_uncertainProteinsFiltered.shape[0] > 0:
            generatePlots.plotPCA(myProteinGroupsFiles_uncertainProteinsFiltered, myId, 'beforeLog2', myBasepath, myOutputFile)

        #### log2 transform data
        myProteinGroupsFiles_log2 = preprocessData.log2Transform(myProteinGroupsFiles_uncertainProteinsFiltered)
        myProteinGroupsFiles_log2.to_csv(
            os.path.join(myBasepath, myOutputFile, myId, 'myProjects_log2_' + myId + '.csv'), index=False)

        if myProteinGroupsFiles_log2.shape[0] > 0:
            generatePlots.plotHistogram(myProteinGroupsFiles_log2, myId, 'beforeNormalization', myBasepath,
                                        myOutputFile)
            generatePlots.plotBoxplot(myProteinGroupsFiles_log2, myId, 'beforeNormalization', myBasepath, myOutputFile)

        ### Assign groups
        ### if ConditionAssignmentFile is given, use defined groups
        if (myConditionAssignmentFiles is not None):
            myProteinGroupsFiles_groupsAssigned = preprocessData.assignGroups(myProteinGroupsFiles_log2, myConditionAssignmentFiles)
            #myProteinGroupsFiles_groupsAssigned = preprocessData.assignGroups(myProteinGroupsFiles_normalized,
            #                                                                  myConditionAssignmentFiles)
            myProteinGroupsFiles_groupsAssigned.to_csv(
                os.path.join(myBasepath, myOutputFile, myId, 'myProjects_groupsAssigned_' + myId + '.csv'), index=False)
        ### if ConditionAssignmentFile is not given, assign group via k-means clustering
        else:
            myProteinGroupsFiles_log2_copy = myProteinGroupsFiles_log2.copy()
            #myProteinGroupsFiles_normalized_copy = myProteinGroupsFiles_normalized.copy()
            ### Impute values from normal distribution as NbClust requires a complete dataset and RF imputation within each group cannot be done as no groups are defined yet
            myProteinGroupsFiles_imputedComplete = preprocessData.imputeFromNormal(myProteinGroupsFiles_log2_copy)
            #myProteinGroupsFiles_imputedComplete = preprocessData.imputeFromNormal(myProteinGroupsFiles_normalized_copy)
            ### if AssignedKFile is provided, use given k
            if (myAssignedKFiles is not None):
                myProteinGroupsFiles_groupsAssigned = preprocessData.predictGroupAssignmentFromGivenK(myProteinGroupsFiles_log2_copy, myProteinGroupsFiles_imputedComplete, myAssignedKFiles)
            #    myProteinGroupsFiles_groupsAssigned = preprocessData.predictGroupAssignmentFromGivenK(
            #        myProteinGroupsFiles_normalized_copy, myProteinGroupsFiles_imputedComplete, myAssignedKFiles)
                myProteinGroupsFiles_groupsAssigned.to_csv(os.path.join(myBasepath, myOutputFile, myId,
                                                                        'myProjects_groupsAssignedFromGivenK_' + myId + '.csv'),
                                                           index=False)
            ### if AssignedKFile is not provieded, predict k using NbClust
            else:
                myProteinGroupsFiles_groupsAssigned = preprocessData.predictGroupAssignment(myProteinGroupsFiles_log2_copy, myProteinGroupsFiles_imputedComplete, myId, myBasepath, myOutputFile)
                #myProteinGroupsFiles_groupsAssigned = preprocessData.predictGroupAssignment(
                #    myProteinGroupsFiles_normalized_copy, myProteinGroupsFiles_imputedComplete, myId, myBasepath,
                #    myOutputFile)
                myProteinGroupsFiles_groupsAssigned.to_csv(
                    os.path.join(myBasepath, myOutputFile, myId, 'myProjects_groupsAssignedPredicted_' + myId + '.csv'),
                    index=False)

        ### Keep only treatments with at least triplicate measurements (as this will be a pre-requisit for the t-test and meaningful outlier detection)
        myProteinGroupsFiles_triplicates = preprocessData.keepOnlyTriplicatesOrMore(myProteinGroupsFiles_groupsAssigned)
        myProteinGroupsFiles_triplicates.to_csv(
            os.path.join(myBasepath, myOutputFile, myId, 'myProjects_triplicates_' + myId + '.csv'), index=False)

        if myProteinGroupsFiles_triplicates.empty == False:

            ### Keep only proteins with at least 70% valid values in one group (to avoid too many false-positives)
            myProteinGroupsFiles_filteredValidValues = preprocessData.filterValid(myProteinGroupsFiles_triplicates)
            myProteinGroupsFiles_filteredValidValues.to_csv(
                os.path.join(myBasepath, myOutputFile, myId, 'myProjects_filteredValidValues_' + myId + '.csv'),
                index=False)

            if myProteinGroupsFiles_filteredValidValues.empty == False:
                generatePlots.plotHistogram(myProteinGroupsFiles_filteredValidValues, myId, 'filterValidValues',
                                            myBasepath, myOutputFile)

                myProteinGroupsFiles_filteredValidValues_copy = myProteinGroupsFiles_filteredValidValues.copy()

                ### impute values from Random Forest imputation (in case there are more than 70% valid values for that protein with this condition) and from downshifted, shrunk normal distribution (otherwise)
                myProteinGroupsFiles_imputed = preprocessData.imputeRFandNormal(myProteinGroupsFiles_filteredValidValues, myId)
                #myProteinGroupsFiles_imputed = preprocessData.imputeRF2(myProteinGroupsFiles_filteredValidValues)
                myProteinGroupsFiles_imputed.to_csv(os.path.join(myBasepath, myOutputFile, myId, 'myProjects_imputedRFAndNormal_' + myId + '.csv'), index=False)

                generatePlots.plotHistogram(myProteinGroupsFiles_imputed, myId, 'imputed', myBasepath, myOutputFile)
                generatePlots.plotBoxplot(myProteinGroupsFiles_imputed, myId, 'imputed', myBasepath, myOutputFile)
                generatePlots.plotPCA(myProteinGroupsFiles_imputed, myId, 'imputed', myBasepath, myOutputFile)
                generatePlots.plotHeatmap(myProteinGroupsFiles_imputed, myId, 'imputed', myBasepath, myOutputFile)

                # ### Detect outliers using PCAGrid
                # ### Remove outliers as long as there are at least still triplicate measurements left
                # myProteinGroupsFiles_outliersRemoved = preprocessData.detectOutliers2(myProteinGroupsFiles_imputed, myId, myBasepath, myOutputFile)
                #
                # myProteinGroupsFiles_outliersRemoved.to_csv(
                #     os.path.join(myBasepath, myOutputFile, myId, 'myProjects_outliersRemoved_' + myId + '.csv'),
                #     index=False)

#                if myProteinGroupsFiles_outliersRemoved.empty == False:
                if myProteinGroupsFiles_imputed.empty == False:

                    print('start linear modeling')

                    ### perform linear modeling to retrieve p-values
#                    myProteinGroupsFiles_lm = downstreamAnalysis.performLinearModeling(myProteinGroupsFiles_outliersRemoved, myProteinGroupsFiles_filteredValidValues_copy, myComparisonsFiles)
                    myProteinGroupsFiles_lm = downstreamAnalysis.performLinearModeling(myProteinGroupsFiles_imputed,
                                                                               myProteinGroupsFiles_filteredValidValues_copy,
                                                                               myComparisonsFiles)

            # myProteinGroupsFiles_tTest = downstreamAnalysis.performTTest(myProteinGroupsFiles_triplicatesAfterOutliers, myProteinGroupsFiles_filteredValidValues_copy, myComparisonsFiles)
                    for key, value in myProteinGroupsFiles_lm.items():
                        #    print(myId)
                        #    print(key)
                        # print(value.shape)
                        value.to_excel(os.path.join(myBasepath, myOutputFile, myId, 'myProjects_linearModelingResults_' + myId + '_' + key + '.xlsx'))

                        generatePlots.plotVolcano(value, myId, key, myBasepath, myOutputFile)

                    ### extract significant (FDR <= 0.05) and log2ratio >= log2(1.5)
                    myProteinGroupsFiles_significant = downstreamAnalysis.extractSignificant(myProteinGroupsFiles_lm)
                    for key, value in myProteinGroupsFiles_significant.items():
                        # print(myId)
                        # print(key)
                        # print(value.shape)
                        value.to_excel(os.path.join(myBasepath, myOutputFile, myId,
                                                    'myProjects_significantProteins_' + myId + '_' + key + '.xlsx'))

                    ### Perform GSEA

                    for db in dbs:
                        myEnrichedPathways = downstreamAnalysis.performGSEA(myProteinGroupsFiles_lm, myProteinGroupsFiles, db)
                        print("pathways:")
                        print(type(myEnrichedPathways))
                        print(type(value))
                        value.to_excel(os.path.join(myBasepath, myOutputFile, myId,'myProjects_GSEAResults_' + db + '_' + myId + '_' + key + '.xlsx'))

                        #for key, value in myEnrichedPathways.items():
                        #    #print(key)
                        #    print(value)
                        #    value.to_csv(os.path.join(myBasepath, myOutputFile, myId, 'myProjects_GSEAResults_'+db+'_'+myId+'_'+key+'.csv'))

# run main

if __name__ == "__main__":
    main()