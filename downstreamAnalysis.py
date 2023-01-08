import itertools
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest as multi
import pandas as pd
import requests
import json
import time
from datetime import datetime
import math
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.conversion import localconverter

class downstreamAnalysis:

    def performLinearModeling(value, myProteinGroupsFiles_filteredValidValues_copy, myComparisonsFiles):

        myListOfDfs_tTest = {}


        base = importr('base')
        stats = importr('stats')

        r = robjects.r

        from rpy2.robjects import r, pandas2ri
        pandas2ri.activate()

        # print(myProteinGroupsFiles_filteredValidValues_copy.iloc[0:3,440:442])

        def computePandFDR(myCondition1, myCondition2, value, myProteinGroupsFiles_filteredValidValues_copy):

            myUniprots = []
            myPs = []
            myLogRatios = []

            #print(value)

            for index, protein in value.iterrows():

                if 'sp|' in protein[0] or 'tr|' in protein[0]:
                    myUniprot = protein[0].split('|')[1]
                else:
                    myUniprot = protein[0].split(';')[0]

                myColsCondition1 = value.columns.get_loc(myCondition1)
                # print(myColsCondition1)
                myIndsCondition1 = ([i for i, e in enumerate(myColsCondition1) if e])
                # print(myIndsCondition1)
                myValuesCondition1 = protein[myIndsCondition1]
                # print(myValuesCondition1)

                myColsCondition2 = value.columns.get_loc(myCondition2)
                # print(myColsCondition2)
                myIndsCondition2 = ([i for i, e in enumerate(myColsCondition2) if e])
                # print(myIndsCondition2)
                myValuesCondition2 = protein[myIndsCondition2]
                # print(myValuesCondition2)

                myColsCondition1_NAs = myProteinGroupsFiles_filteredValidValues_copy.columns.get_loc(myCondition1)
                myIndsCondition1_NAs = ([i for i, e in enumerate(myColsCondition1_NAs) if e])
                # print(index)
                # print(myIndsCondition1_NAs)
                # print(myProteinGroupsFiles_filteredValidValues_copy.loc[index].iloc[myIndsCondition1_NAs])
                myValuesCondition1_NAs = myProteinGroupsFiles_filteredValidValues_copy.loc[index].iloc[
                    myIndsCondition1_NAs]  ### instead of chaining loc and iloc, one could use only iloc with Index.get_indexer for the row
                # print(myValuesCondition1_NAs)

                myColsCondition2_NAs = myProteinGroupsFiles_filteredValidValues_copy.columns.get_loc(myCondition2)
                myIndsCondition2_NAs = ([i for i, e in enumerate(myColsCondition2_NAs) if e])
                myValuesCondition2_NAs = myProteinGroupsFiles_filteredValidValues_copy.loc[index].iloc[
                    myIndsCondition2_NAs]
                # print(myValuesCondition2_NAs)

                myNApercentageCondition1 = myValuesCondition1_NAs.isnull().sum() / len(myValuesCondition1_NAs)
                myNApercentageCondition2 = myValuesCondition2_NAs.isnull().sum() / len(myValuesCondition2_NAs)
                # print(myNApercentageCondition1)
                # print(myNApercentageCondition2)

                if myNApercentageCondition1 <= 0.3 or myNApercentageCondition2 <= 0.3:

                    #print(pd.concat([myValuesCondition1,myValuesCondition2]))
                    #print(pd.concat([pd.Series([0]*len(myValuesCondition1)), pd.Series([1]*len(myValuesCondition2))]))

                    #myData = pd.DataFrame({'Abundance':pd.concat([myValuesCondition1,myValuesCondition2],ignore_index=True), 'Group':pd.concat([pd.Series([0]*len(myValuesCondition1)), pd.Series([1]*len(myValuesCondition2))],ignore_index=True)})
                    myData = pd.DataFrame({'Abundance': [round (elem,5) for elem in pd.concat([myValuesCondition1, myValuesCondition2], ignore_index=True)], 'Group': pd.concat([pd.Series([0] * len(myValuesCondition1)), pd.Series([1] * len(myValuesCondition2))], ignore_index=True)})

                    #print('start lm')

                    myLm = stats.lm('Abundance ~ Group', data=myData)

                    #print('Result of lm:')
                    #print(myLm)
                    #print(base.summary(myLm))
                    #print(myData)
                    #print(base.summary(myLm).rx2('coefficients'))

                    if (any(myValuesCondition1.notna()) and any(myValuesCondition2.notna())):
                        myP = base.summary(myLm).rx2('coefficients')[1][3]
                        import math
                        if (math.isnan(myP)):
                            print(myUniprot)
                            print(myValuesCondition1)
                            print(myValuesCondition2)
                            print(base.summary(myLm))
                        #myP = stats.ttest_ind(myValuesCondition1, myValuesCondition2)[1]
                        myLogRatio = np.mean(myValuesCondition1) - np.mean(myValuesCondition2)

                        if not (math.isnan(myP)):
                            myUniprots.append(myUniprot)
                            myPs.append(myP)
                            myLogRatios.append(myLogRatio)

            print(myPs)
            myFDRs = multi.fdrcorrection(myPs)[1]
            print(myFDRs)

            myResult_dict = {'ProteinID': myUniprots, 'Pvalue': myPs, 'FDR': myFDRs, 'LogRatio': myLogRatios}
            myResult = pd.DataFrame(myResult_dict)

            return (myResult)

        myLFQCols_unique = set([col for col in value.columns if 'LFQ intensity_' in col])

        if len(myLFQCols_unique) > 1:  ### may be deleted later

            if myComparisonsFiles is not None:

                print(myComparisonsFiles)

                for index, row in myComparisonsFiles.iterrows():

                    print(row)

                    myCondition1 = 'LFQ intensity_' + row[0]
                    myCondition2 = 'LFQ intensity_' + row[1]

                    #print(myCondition1)
                    #print(myCondition2)
                    #print(myLFQCols_unique)

                    if myCondition1 in myLFQCols_unique and myCondition2 in myLFQCols_unique:
                        myResult = computePandFDR(myCondition1, myCondition2, value, myProteinGroupsFiles_filteredValidValues_copy)

                        myListOfDfs_tTest[myCondition1 + '_' + myCondition2] = myResult

            elif 'LFQ intensity_control' in value.columns.values:

                for name in myLFQCols_unique:

                    if name != 'LFQ intensity_control':

                        myCondition1 = name
                        myCondition2 = 'LFQ intensity_control'

                        if myCondition1 in myLFQCols_unique:
                            # print(myCondition1)
                            myResult = computePandFDR(myCondition1, myCondition2, value,
                                                      myProteinGroupsFiles_filteredValidValues_copy)

                            myListOfDfs_tTest[myCondition1 + '_' + myCondition2] = myResult

            else:

                for pair in itertools.combinations(myLFQCols_unique, 2):
                    # print(pair)

                    myCondition1 = pair[0]
                    myCondition2 = pair[1]

                    myResult = computePandFDR(myCondition1, myCondition2, value,
                                              myProteinGroupsFiles_filteredValidValues_copy)

                    myListOfDfs_tTest[myCondition1 + '_' + myCondition2] = myResult

        return (myListOfDfs_tTest)

    def extractSignificant(myProteinGroupsFiles_tTest):

        myListOfDfs_significant = {}

        for key, value in myProteinGroupsFiles_tTest.items():
            mySignificantProteins_tmp = value.loc[value['FDR'] <= 0.05]
            mySignificantProteins = mySignificantProteins_tmp.loc[abs(mySignificantProteins_tmp['LogRatio']) >= math.log2(1.5)]

            myListOfDfs_significant[key] = mySignificantProteins

        return (myListOfDfs_significant)

    def performGSEA(myProteomicsData, myProteinGroupsFiles, db):

        #print(myProteinGroupsFiles_significant)

        r = robjects.r

        #r('install.packages("BiocManager", repos="http://cran.r-project.org")')
        #r('BiocManager::install("gage")')
        # robjects.r('BiocManager::install("MSstats")')

        fgsea = importr('fgsea')
        gage = importr('gage')
        base = importr('base')
        UniProt = importr('UniProt.ws')
        msigdbr = importr('msigdbr')
        AnnotationDbi = importr('AnnotationDbi')

        from rpy2.robjects import r, pandas2ri
        pandas2ri.activate()


        myGSEAResults = {}

        myKEGGsInflammation = pd.read_excel("C:/Users/1/Desktop/PROTEOMAS/LungInflammation_proteins.xlsx")
        print(myKEGGsInflammation)

        for label in myProteomicsData.keys():

            #print(myProteomicsData[label])

            myProteomicsData_Label = myProteomicsData[label]
            print(myProteomicsData_Label.head)


            # Map UniProt to gene names

            humanUp = UniProt.UniProt_ws(taxId=9606)

            # Map mouse or rat gene names to human ones

            mouse2human = msigdbr.msigdbr(species='mouse', category=db)
            rat2human = msigdbr.msigdbr(species='rat', category=db)

            myProteinGeneList = AnnotationDbi.select(humanUp, keys=myProteomicsData_Label['ProteinID'], columns = base.c('gene_primary'), keytype = 'UniProtKB')
            print(myProteinGeneList)

            myMergedFile = pd.merge(myProteomicsData_Label, myProteinGeneList, left_on='ProteinID', right_on='From', how='left')
            print(myMergedFile)

            myGeneNames_tmp = myMergedFile['Gene.Names..primary.']
            myGeneNames_single = myGeneNames_tmp.replace('-(\\d+)', '', regex=True)
            myGeneNames_withoutNumber = myGeneNames_single.replace('-(\\d+)', '', regex=True)

            myMergedFile["Gene.Names..primary."] = myGeneNames_withoutNumber

            print(myMergedFile.columns)
            print(myMergedFile.shape)

            myMergedFile = myMergedFile[myMergedFile['Gene.Names..primary.'].str.contains('None') == False]
                #myMergedFile[~myMergedFile['Gene.Names..primary.'].isin(['None','NA'])]
            print(myMergedFile)
            myMergedFile['colForSorting'] = robjects.vectors.FloatVector(base.as_numeric(-np.log10(myMergedFile['FDR']) * np.abs(myMergedFile['LogRatio'])))
            myMergedFile.sort_values(['colForSorting'],ascending=False)
            myMergedFile.drop_duplicates(subset='Gene.Names..primary.')

            print(myMergedFile)


            #myMergedFile_withoutDuplicates = myMergedFile.drop_duplicates(subset=['Gene.Names..primary.'])

            #myCombinedResult = robjects.ListVector({myMergedFile['Gene.Names..primary.']: -np.log10(myMergedFile['FDR']) * np.abs(myMergedFile['LogRatio'])})

            myCombinedResult = robjects.vectors.FloatVector(base.sort(base.as_numeric(myMergedFile['colForSorting']), decreasing=False))
            myCombinedResult.names = myMergedFile['Gene.Names..primary.']

            print(myCombinedResult)

            #myCombinedResult = myCombinedResult[-base.which(base.is_na(myCombinedResult.names))]
            #myCombinedResult_withoutDuplicates = myCombinedResult_sorted[set(myCombinedResult_sorted.names)]

            # Distinguish species

            myHumanLength = len([i for i in myMergedFile['Gene.Names..primary.'] if i in mouse2human['human_gene_symbol']])
            myMouseLength = len([i for i in myMergedFile['Gene.Names..primary.'] if i in mouse2human['gene_symbol']])
            myRatLength = len([i for i in myMergedFile['Gene.Names..primary.'] if i in rat2human['gene_symbol']])

            if (myHumanLength >= myMouseLength and myHumanLength >= myRatLength):

                msigdbr_list = robjects.ListVector(base.split(x = mouse2human['human_gene_symbol'], f = mouse2human['gs_name']))
                myKEGGsInflammation_names = AnnotationDbi.select(humanUp, keys=myKEGGsInflammation['UniProtID human'], columns=base.c("gene_primary"), keytype="UniProtKB")

            elif (myMouseLength > myHumanLength and myMouseLength > myRatLength):

                msigdbr_list = robjects.ListVector(base.split(x = mouse2human['gene_symbol'], f = mouse2human['gs_name']))
                myKEGGsInflammation_names = AnnotationDbi.select(humanUp, keys=myKEGGsInflammation['UniProtID mouse'], columns=base.c("gene_primary"), keytype="UniProtKB")

            elif (myRatLength > myHumanLength and myRatLength > myMouseLength):
                msigdbr_list = robjects.ListVector(base.split(x=rat2human['gene_symbol'], f=rat2human['gs_name']))
                myKEGGsInflammation_names = AnnotationDbi.select(humanUp, keys=myKEGGsInflammation['UniProtID rat'], columns=base.c("gene_primary"), keytype="UniProtKB")

                print(msigdbr_list)

            msigdbr_list_2 = robjects.ListVector({'Lung Inflammation Key Event': myKEGGsInflammation_names['Gene.Names..primary.']})
            print(msigdbr_list_2)

            msigdbr_list_added = robjects.ListVector(base.c(msigdbr_list,msigdbr_list_2))
            #msigdbr_list.append(myKEGGsInflammation_names['Gene.Names..primary.'],index='Lung Inflammation Key Event')
            #base.append(msigdbr_list.names,'Lung Inflammation Key Event')

            print(msigdbr_list_added)
            #print(len(myCombinedResult))


            ### test
            examplePathways = fgsea.__rdata__.fetch('examplePathways')
            print(examplePathways['examplePathways'])

            exampleRanks = fgsea.__rdata__.fetch('exampleRanks')
            print(exampleRanks['exampleRanks'])

            test = robjects.vectors.FloatVector([0.001,0.004,0.01,0.1])
            test.names = ['234695','19724','67610','12763']
            print(test)

            #fgseaRes = fgsea.fgsea(pathways=examplePathways['examplePathways'][1],stats=exampleRanks['exampleRanks'][1],minSize=15,maxSize=500)
            #print(fgseaRes)

            #msigdbr_list['Lung Inflammation Key Event'].append(myKEGGsInflammation_names['Gene.Names..primary.'])
            #res = fgsea.fgsea(pathways=msigdbr_list_added, stats=myCombinedResult, minSize=15, maxSize=600, scoreType='pos')

            pandas2ri.deactivate()
            res = fgsea.fgsea(pathways=msigdbr_list_added, stats=myCombinedResult, minSize=15, maxSize=600, scoreType='pos')
            print(res)

            res_df = pd.DataFrame(res).transpose()
            res_df.columns = base.colnames(res)

            myGSEAResults[label] = res_df

        return(myGSEAResults)