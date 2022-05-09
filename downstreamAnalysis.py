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

class downstreamAnalysis:

    def performTTest(value, myProteinGroupsFiles_filteredValidValues_copy, myComparisonsFiles):

        myListOfDfs_tTest = {}

        #print(myProteinGroupsFiles_filteredValidValues_copy.iloc[0:3,440:442])

        def computePandFDR(myCondition1, myCondition2, value, myProteinGroupsFiles_filteredValidValues_copy):

            ### Falsche Kontrolle bei PXD016148

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
                #print(myColsCondition1)
                myIndsCondition1 = ([i for i, e in enumerate(myColsCondition1) if e])
                #print(myIndsCondition1)
                myValuesCondition1 = protein[myIndsCondition1]
                #print(myValuesCondition1)

                myColsCondition2 = value.columns.get_loc(myCondition2)
                #print(myColsCondition2)
                myIndsCondition2 = ([i for i, e in enumerate(myColsCondition2) if e])
                #print(myIndsCondition2)
                myValuesCondition2 = protein[myIndsCondition2]
                #print(myValuesCondition2)

                myColsCondition1_NAs = myProteinGroupsFiles_filteredValidValues_copy.columns.get_loc(myCondition1)
                myIndsCondition1_NAs = ([i for i, e in enumerate(myColsCondition1_NAs) if e])
                #print(index)
                #print(myIndsCondition1_NAs)
                #print(myProteinGroupsFiles_filteredValidValues_copy.loc[index].iloc[myIndsCondition1_NAs])
                myValuesCondition1_NAs = myProteinGroupsFiles_filteredValidValues_copy.loc[index].iloc[myIndsCondition1_NAs] ### instead of chaining loc and iloc, one could use only iloc with Index.get_indexer for the row
                #print(myValuesCondition1_NAs)

                myColsCondition2_NAs = myProteinGroupsFiles_filteredValidValues_copy.columns.get_loc(myCondition2)
                myIndsCondition2_NAs = ([i for i, e in enumerate(myColsCondition2_NAs) if e])
                myValuesCondition2_NAs = myProteinGroupsFiles_filteredValidValues_copy.loc[index].iloc[myIndsCondition2_NAs]
                #print(myValuesCondition2_NAs)

                myNApercentageCondition1 = myValuesCondition1_NAs.isnull().sum() / len(myValuesCondition1_NAs)
                myNApercentageCondition2 = myValuesCondition2_NAs.isnull().sum() / len(myValuesCondition2_NAs)
                #print(myNApercentageCondition1)
                #print(myNApercentageCondition2)

                if myNApercentageCondition1 <= 0.3 or myNApercentageCondition2 <= 0.3:
                    myP = stats.ttest_ind(myValuesCondition1, myValuesCondition2)[1]
                    myLogRatio = np.mean(myValuesCondition1)-np.mean(myValuesCondition2)
                #else:
                    #print(myValuesCondition1)
                    #print(myValuesCondition2)
                    #print(myNApercentageCondition1)
                    #print(myNApercentageCondition2)
                #    myP = np.nan
                #    myLogRatio = np.nan

                    myUniprots.append(myUniprot)
                    myPs.append(myP)
                    myLogRatios.append(myLogRatio)


            #print(myPs)
            myFDRs = multi.fdrcorrection(myPs)[1]
            #print(myFDRs)

            myResult_dict = {'ProteinID': myUniprots, 'Pvalue': myPs, 'FDR': myFDRs, 'LogRatio': myLogRatios}
            myResult = pd.DataFrame(myResult_dict)

            return(myResult)


        myLFQCols_unique = set([col for col in value.columns if 'LFQ' in col])

        if len(myLFQCols_unique) > 1: ### may be deleted later

            if myComparisonsFiles is not None:

                #print(myComparisonsFiles)

                for index, row in myComparisonsFiles.iterrows():

                    myCondition1 = 'LFQ_' + row[0]
                    myCondition2 = 'LFQ_' + row[1]

                    if myCondition1 in myLFQCols_unique and myCondition2 in myLFQCols_unique:

                        myResult = computePandFDR(myCondition1, myCondition2, value, myProteinGroupsFiles_filteredValidValues_copy)

                        myListOfDfs_tTest[myCondition1 + '_' + myCondition2] = myResult

            elif 'LFQ_control' in value.columns.values:

                for name in myLFQCols_unique:

                    if name != 'LFQ_control':

                        myCondition1 = name
                        myCondition2 = 'LFQ_control'

                        if myCondition1 in myLFQCols_unique:
                            #print(myCondition1)
                            myResult = computePandFDR(myCondition1, myCondition2, value, myProteinGroupsFiles_filteredValidValues_copy)

                            myListOfDfs_tTest[myCondition1 + '_' + myCondition2] = myResult

            else:

                for pair in itertools.combinations(myLFQCols_unique,2):

                    #print(pair)

                    myCondition1 = pair[0]
                    myCondition2 = pair[1]

                    myResult = computePandFDR(myCondition1, myCondition2, value, myProteinGroupsFiles_filteredValidValues_copy)

                    myListOfDfs_tTest[myCondition1 + '_' + myCondition2] = myResult

        return(myListOfDfs_tTest)


    def extractSignificant(myProteinGroupsFiles_tTest):

        myListOfDfs_significant = {}

        for key, value in myProteinGroupsFiles_tTest.items():

            mySignificantProteins_tmp = value.loc[value['FDR'] <= 0.05]
            mySignificantProteins = mySignificantProteins_tmp.loc[abs(mySignificantProteins_tmp['LogRatio']) >= math.log2(1.5)]

            myListOfDfs_significant[key] = mySignificantProteins

        return(myListOfDfs_significant)


    def querySTRING(myProteinGroupsFiles_significant):

        myListOfDfs_STRING = {}

        for key, value in myProteinGroupsFiles_significant.items():

            if value.shape[0] > 0:

                ##################### 1) Mapping protein IDs to STRING IDs ####################
                ## For a given list of proteins the script resolves them
                ## (if possible) to the best matching STRING identifier
                ## and prints out the mapping on screen in the TSV format
                ###############################################################################

                string_api_url = "https://string-db.org/api"
                output_format = "tsv"
                method = "get_string_ids"

                print(datetime.now())

                stringIDs = []  # list with all String IDs, which will be input for next step

                # convert dataframe to list
                df2list = value['ProteinID'].tolist()
                #print(df2list)

                ## Set parameters
                for ids in df2list:

                    params = {

                        "identifiers": ids,  # your protein list
                        # "species": 9606,  # species NCBI identifier (HUMAN)
                        "limit": 5,  # only one (best) identifier per input protein
                        "echo_query": 1,  # see your input identifiers in the output
                        "caller_identity": "www.somewhere_app.org"  # your app name

                        }

                    ## Construct URL
                    request_url = "/".join([string_api_url, output_format, method])

                    ## Call STRING
                    results = requests.post(request_url, data=params)

                    ## Read and parse the results
                    for line in results.text.strip().split("\n"):
                        l = line.split("\t")
                        #print(l)
                        # time.sleep(2)
                        stringIDs.append(l[2])

                # delete the headers (stringIds) and save only the Human StringIDs
                sId = "stringId"

                # use a loop for deleting the headers
                try:
                    while True:
                        stringIDs.remove(sId)
                except ValueError:
                    pass

                print(stringIDs)

                if len(stringIDs) > 0:

                    ########### 2) Get Enrichments of the identified protein IDs ##################
                    ## The following script retrieves and prints out significantly enriched
                    ## (FDR < 1%) GO Processes for the given set of proteins.
                    ###############################################################################

                    string_api_url = "https://string-db.org/api"
                    output_format = "json"
                    method = "enrichment"

                    ## Construct the request
                    request_url = "/".join([string_api_url, output_format, method])

                    ## Set parameters
                    params = {

                        "identifiers": "%0d".join(stringIDs),  # your protein
                        # "species": 9606,  # species NCBI identifier
                        "caller_identity": "www.something_app.org"  # your app name

                    }

                    ## Call STRING
                    response = requests.post(request_url, data=params)
                    time.sleep(2)

                    ## Read and parse the results
                    # print(response.text)
                    data = json.loads(response.text)
                    time.sleep(2)

                    #print(data)

                    myCategories = []
                    myTerms = []
                    myFDRs = []
                    myPs = []
                    myDescriptions = []

                    #print(data)

                    for row in data:

                        #print(row)
                        #print(row['category'])
                        #print(float(row['fdr']))

                        myCategories.append(row['category'])
                        myTerms.append(row['term'])
                        myFDRs.append(float(row['fdr']))
                        myPs.append(float(row['p_value']))
                        myDescriptions.append(row['description'])

                    mySTRINGData_dict = {'Category': myCategories, 'Term': myTerms, 'FDR': myFDRs, 'Pvalue': myPs, 'Description': myDescriptions}
                    mySTRINGData = pd.DataFrame(mySTRINGData_dict)

                    myListOfDfs_STRING[key] = mySTRINGData

        return(myListOfDfs_STRING)