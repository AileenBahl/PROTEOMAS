import numpy as np
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
import pandas as pd
import re
import os
from generatePlots import generatePlots
import qnorm

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)
packnames = ('rrcov', 'ggplot2', 'missForest','rrcovHD','Uniprot.ws','msigdbr','AnnotationDbi')

names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))



class preprocessData:


    def filterUnnecessaryProteins(value,myColname):

        myColnameComplete = [s for s in value.columns if myColname.lower() in s.lower()]
        #print(myColnameComplete[0]) ### Can this be longer than 1? E.g Contaminant and Potential contaminant in the same proteinGroups file?
        myDf_new = value.loc[value[myColnameComplete[0]] != '+']

        return myDf_new


    def filterUncertainProteins(value):

        myDf_new = value.loc[(value['Peptides'] >= 2) & (value['Unique peptides'] >= 1)]
        print(myDf_new.shape)

        return myDf_new


    def log2Transform(value):

        myLFQCols = [col for col in value.columns if 'LFQ intensity ' in col]

        newValue = value.apply(lambda x: np.log2(x.replace(0, np.nan)) if x.name in myLFQCols else x)

        print("log2:")
        print(newValue)

        return newValue


    def VSNNormalize(value):

        #r = robjects.r

        #r('install.packages("BiocManager", repos="http://cran.r-project.org")')
        #r('BiocManager::install("vsn")')

        vsn = importr('vsn')
        base = importr('base')
        Biobase = importr('Biobase')

        from rpy2.robjects import r, pandas2ri
        pandas2ri.activate()

        set_seed = r('set.seed')
        set_seed(1)

        myNewValue_LFQ = pd.DataFrame(value.filter(regex='LFQ intensity '))
        print(myNewValue_LFQ)

        myNewValue_LFQ_matrix = base.data_matrix(myNewValue_LFQ/100000000)
        myNewValue_LFQ_matrix = myNewValue_LFQ_matrix * 100000000
        print(myNewValue_LFQ_matrix)

        fit = vsn.vsn2(myNewValue_LFQ_matrix)
        ynorm = vsn.predict(fit, myNewValue_LFQ_matrix)

        ynorm[ynorm < 0] = np.nan

        myLFQCols_logical = np.array([l.startswith('LFQ intensity ') for l in value.columns.values])
        myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]
        print("test here")
        print(myLFQCols)

        #value.iloc[:, myLFQCols] = myNewValue_LFQ
        value.iloc[:, myLFQCols] = ynorm
        print(value)

        return value


    def QNNormalize(value):

        base = importr('base')

        from rpy2.robjects import r, pandas2ri
        pandas2ri.activate()

        set_seed = r('set.seed')
        set_seed(1)

        myNewValue_LFQ = pd.DataFrame(value.filter(regex='LFQ intensity '))
        print(myNewValue_LFQ)

        myNewValue_LFQ_matrix = base.data_matrix(myNewValue_LFQ)

        ynorm = qnorm.quantile_normalize(myNewValue_LFQ_matrix, axis=1)

        myLFQCols_logical = np.array([l.startswith('LFQ intensity ') for l in value.columns.values])
        myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]
        print("test here")
        print(myLFQCols)

        #value.iloc[:, myLFQCols] = myNewValue_LFQ
        value.iloc[:, myLFQCols] = ynorm
        print(value)

        return value


    def assignGroups(value, myConditionAssignmentFiles):

        #if (myConditionAssignmentFiles is not None):
            # print(myConditionAssignmentFiles.columns)
            # print(value.columns)
        myLFQCols = [col for col in value.columns if col in myConditionAssignmentFiles.columns]
        myNewNames = ['LFQ intensity_' + str(s) for s in myConditionAssignmentFiles.iloc[0]]
            #myNewNames = ['LFQ intensity_' + str(s) for s in a]
            #print(myNewNames)
        #else:
        #    myLFQCols = [col for col in value.columns if 'LFQ intensity ' in col]
        #    myNewNames = ['LFQ intensity_' + str(s) for s in ['unknown'] * len(myLFQCols)]
            #print(myNewNames)

        #print("here:")
        #print(len(myLFQCols))
        #print(myLFQCols.shape[1])
        #print(myNewNames)

        value.rename(columns=dict(zip(myLFQCols, myNewNames)), inplace=True)

        #print(value)

        return value


    def predictGroupAssignment(value, myProteinGroupsFiles_imputedComplete, myId, myBasepath, myOutputFile):

        print('Start NbClust')

        ### Add option of obtaining k from an input file

        NbClust = importr('NbClust')
        base = importr('base')
        stats = importr('stats')

        r = robjects.r

        from rpy2.robjects import r, pandas2ri
        pandas2ri.activate()

        set_seed = r('set.seed')
        set_seed(1)

        print(myProteinGroupsFiles_imputedComplete.columns.values)
        myNewValue_LFQ = myProteinGroupsFiles_imputedComplete.filter(regex='LFQ intensity ')
        print(myNewValue_LFQ)

        # possible solution: imputation with one group (more than 70% overall = RF and less = normal distribution), then assign groups, then filtering and imputation again for the new groups
        #print(np.transpose(myNewValue_LFQ).shape)
        #print(base.nrow(np.transpose(myNewValue_LFQ)) - 1)
        #print(base.nrow(np.transpose(myNewValue_LFQ))/3+base.max(base.c(2,base.nrow(np.transpose(myNewValue_LFQ))/10)))

        print("start")
        print(np.transpose(myNewValue_LFQ))
        print(base.round(base.nrow(np.transpose(myNewValue_LFQ))/3+base.max(base.c(2,base.nrow(np.transpose(myNewValue_LFQ))/10))))
        print("NbClust")
        myNbClustResult = NbClust.NbClust(np.transpose(myNewValue_LFQ), min_nc=2, max_nc = base.round(base.nrow(np.transpose(myNewValue_LFQ))/3+base.max(base.c(2,base.nrow(np.transpose(myNewValue_LFQ))/10))), distance='euclidean', method='kmeans', index='cindex')
        print(myNbClustResult)
        #myNbClustResult = NbClust.NbClust(np.transpose(myNewValue_LFQ), min_nc=2, max_nc=7, distance='euclidean', method='kmeans', index='cindex')
        print("finish")

        #### Should we set the min_nc to 1 and max_nc to samples/3 as we'd expect at least triplicates or to 10 as with more conditions we anyways don't expect to get the clustering right?
        # myClusterNumber = myNbClustResult[[2]][1]

        ##### instead look into best partition argument
        # km = stats.kmeans(np.transpose(myNewValue_LFQ),myClusterNumber, nstart=30)
        # fviz_cluster(km, np.transpose(myNewValue_LFQ))
        #print(myNbClustResult.rx2('Best.partition'))
        print(myNbClustResult)
        k = myNbClustResult.rx2('Best.nc')[0]
        print("k:")
        print(k)
        km = myNbClustResult.rx2('Best.partition')
        #print(km)

        myNewValue_LFQ.to_csv("C:/Users/1/Desktop/PROTEOMAS/test_km.csv",index=False)

        generatePlots.plotClusterplotForKNN(myNewValue_LFQ, k, myId, myBasepath, myOutputFile)

        myNewNames = ['LFQ intensity_' + str(s) for s in km]
        #print(myNewNames)

        myLFQCols_logical = np.array([l.startswith('LFQ intensity ') for l in value.columns.values])
        myLFQColsInOriginal = np.arange(0, value.shape[1])[myLFQCols_logical]
        #print(myLFQColsInOriginal)
        #print(dict(zip(myProteinGroupsFiles_groupsAssigned.columns[myLFQColsInOriginal], myNewNames)))

        value.rename(columns=dict(zip(value.columns[myLFQColsInOriginal], myNewNames)), inplace=True)
        #print(myProteinGroupsFiles_groupsAssigned.columns)

        return value


    def predictGroupAssignmentFromGivenK(value, myProteinGroupsFiles_imputedComplete, myK):

        print('Assign from given k')

        base = importr('base')
        stats = importr('stats')

        r = robjects.r

        from rpy2.robjects import r, pandas2ri
        pandas2ri.activate()

        set_seed = r('set.seed')
        set_seed(1)

        #print(value.shape)
        #print(value.columns.values)
        myNewValue_LFQ = myProteinGroupsFiles_imputedComplete.filter(regex='LFQ intensity ')
        print(myNewValue_LFQ)
        #print(myK.squeeze())
        #print(type(np.transpose(myNewValue_LFQ)))
        #print(type(myK.squeeze()))

        km = stats.kmeans(np.transpose(myNewValue_LFQ),myK.squeeze(), nstart=30)
        km_clusters = km.rx2('cluster')
        #print(km)

        myNewNames = ['LFQ intensity_' + str(s) for s in km_clusters]
        #print(myNewNames)

        myLFQCols_logical = np.array([l.startswith('LFQ intensity ') for l in value.columns.values])
        myLFQColsInOriginal = np.arange(0, value.shape[1])[myLFQCols_logical]
        #print(myLFQColsInOriginal)
        #print(dict(zip(myProteinGroupsFiles_groupsAssigned.columns[myLFQColsInOriginal], myNewNames)))

        value.rename(columns=dict(zip(value.columns[myLFQColsInOriginal], myNewNames)), inplace=True)
        #print(myProteinGroupsFiles_groupsAssigned.columns)

        return value


    def keepOnlyTriplicatesOrMore(value):

        myLFQCols_all = [col for col in value.columns if 'LFQ intensity_' in col]
        #print(myLFQCols_all)

        myNonTriplicateNames = set([i for i in myLFQCols_all if myLFQCols_all.count(i)<3])
        #print(myNonTriplicateNames)

        newValue = value.drop(myNonTriplicateNames, axis=1)
        test=[col for col in newValue.columns if 'LFQ intensity_' in col]
        #print(test)

        return newValue


    def filterValid(value):

        myLFQCols_unique = set([col for col in value.columns if 'LFQ intensity_' in col])

        myNAList = []
        for name in myLFQCols_unique:
            myCurrentCols = value[name]
            myNAsum = myCurrentCols.isnull().sum(axis=1)
            myNApercentage = myNAsum/myCurrentCols.shape[1]
            myNAList.append(myNApercentage)

        myNADf = pd.DataFrame(myNAList)
        myNADf_min = myNADf.min(axis=0)
        #print(type(myNADf_max))

        myNADf_valid = [i for i, v in enumerate(myNADf_min) if v <= 0.3]
        #myNADf_valid = [i for i, v in enumerate(myNADf_min) if v <= 0.4]

        newValue = value.iloc[myNADf_valid, :]

        return newValue


    def imputeFromNormal(value):

        myLFQCols_logical = np.array([l.startswith('LFQ intensity ') for l in value.columns.values])
        myLFQCols = np.arange(0,value.shape[1])[myLFQCols_logical]
        #print(myLFQCols)

        for i in myLFQCols:

            #print(i)
            myValues = value.iloc[:, i]
            myNAValues = myValues.isnull().values
            myMean = myValues.mean()
            mySd = myValues.std()

            myMean_down = myMean - 1.8 * mySd
            mySd_down = mySd * 0.3

            myRandomValues = np.random.normal(loc=myMean_down, scale=mySd_down, size=value.shape[0])

            value.iloc[myNAValues, i] = myRandomValues[myNAValues]

        return value


    def imputeRF(value):

        missForest = importr('missForest')
        base = importr('base')

        r = robjects.r

        from rpy2.robjects import r, pandas2ri
        pandas2ri.activate()

        set_seed = r('set.seed')
        set_seed(1)

        myNewValue_LFQ = value.filter(regex='LFQ intensity_')
        print(myNewValue_LFQ.shape)

        print("RF imputation")
        #print(missForest.missForest(myNewValue_LFQ))
        #print(missForest.missForest(myNewValue_LFQ)[0])

        newValue = missForest.missForest(myNewValue_LFQ)[0]
        print(newValue.shape)

        return newValue


    def imputeRF2(value):

        myLFQCols_unique = set([col for col in value.columns if 'LFQ intensity_' in col])
        print(myLFQCols_unique)
        #print(value.shape)
        value_copy = value.copy()

        for name in myLFQCols_unique:

            print(name)

            myNewValue_LFQ = value_copy.loc[:,name]
            #print('myNewValue_LFQ')
            #print(myNewValue_LFQ.shape)

            if myNewValue_LFQ.shape[0] > 0 and any(myNewValue_LFQ.isnull()):

                print("now:")
                #print(myNewValue_LFQ.shape)
                #print(any(myNewValue_LFQ.isnull()))

                missForest = importr('missForest')
                base = importr('base')

                r = robjects.r

                from rpy2.robjects import r, pandas2ri
                pandas2ri.activate()

                set_seed = r('set.seed')
                set_seed(1)

                #test=missForest.missForest(np.transpose(myNewValue_LFQ), maxiter = 3, ntree=30)
                #print(test)
                print(np.transpose(myNewValue_LFQ))
                myImputedRows = missForest.missForest(np.transpose(myNewValue_LFQ), maxiter=3, ntree=30)[0]
                print(myImputedRows)

                test = np.transpose(myImputedRows)
                #print(test.shape)

                value.loc[:, name] = test
                #print(value.shape)

        return value


    def imputeRFandNormal(value, myId):

        ### RF imputation on NAs with more than 70% valid values in that group and row

        myLFQCols_unique = set([col for col in value.columns if 'LFQ intensity_' in col])
        #print(myLFQCols_unique)
        #print(value.shape)
        value_copy = value.copy()

        for name in myLFQCols_unique:

            print(name)
            myCurrentCols = value[name]
            print(myCurrentCols)
            myNAsum = myCurrentCols.isnull().sum(axis=1)
            print(myNAsum)
            myNApercentage = myNAsum / myCurrentCols.shape[1]
            print(myNApercentage[0:10])
            myNADf = pd.DataFrame(myNApercentage)
            print(myNADf)
            #print(type(myNADf))

            myNADf_RF = [i for i, v in enumerate(myNADf.iloc[:,0]) if 0 < v <= 0.3]
            #myNADf_RF = [i for i, v in enumerate(myNADf.iloc[:,0]) if v <= 0.4]
            print(myNADf_RF)

            #newValue = value.iloc[myNADf_RF, :]

            #print(name)

            myNewValue_LFQ_tmp = value_copy.loc[:,name]
            myNewValue_LFQ = myNewValue_LFQ_tmp.iloc[myNADf_RF, :]
            print('myNewValue_LFQ')
            print(myNewValue_LFQ.shape)

            if myNewValue_LFQ.shape[0] > 0 and any(myNewValue_LFQ.isnull()):

                print("now:")
                print(myNewValue_LFQ.shape[0])
                print(any(myNewValue_LFQ.isnull()))

                missForest = importr('missForest')
                base = importr('base')

                r = robjects.r

                from rpy2.robjects import r, pandas2ri
                pandas2ri.activate()

                set_seed = r('set.seed')
                set_seed(1)

                #test=missForest.missForest(np.transpose(myNewValue_LFQ), maxiter = 3, ntree=30)
                #print(test)
                #print(np.transpose(myNewValue_LFQ))
                myImputedRows = missForest.missForest(np.transpose(myNewValue_LFQ), maxiter=3, ntree=30)[0]

                test = np.transpose(myImputedRows)
                myNewValue_LFQ_tmp.columns = test.columns
                test.index = test.index.astype(int)

                myNewValue_LFQ_tmp.loc[test.index, :] = test[:]

                print("after RF imputation:")
                print(myNewValue_LFQ_tmp)

                myLFQCols_logical = value.columns.get_loc(name)
                #    print(myLFQCols_logical)
                myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]

                #print(myLFQCols)

                value.iloc[:, myLFQCols] = myNewValue_LFQ_tmp
                #print(value)


        # imputation from normal distribution for all others

        myLFQCols_logical = np.array([l.startswith('LFQ intensity_') for l in value.columns.values])
        myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]

        for i in myLFQCols:
            #print(i)
            myValues = value.iloc[:, i]
            myNAValues = myValues.isnull().values
            myMean = myValues.mean()
            mySd = myValues.std()

            myMean_down = myMean - 1.8 * mySd
            mySd_down = mySd * 0.3

            myRandomValues = np.random.normal(loc=myMean_down, scale=mySd_down, size=value.shape[0])

            print("after norm imputation:")
            print(myRandomValues)

            value.iloc[myNAValues, i] = myRandomValues[myNAValues]

            #print(value)


        return value


    def detectOutliers(value, myId, myBasepath, myOutputFile):

        myNewValue_LFQ = value.filter(regex='LFQ intensity_')

        rrcov = importr('rrcov')
        base = importr('base')
        stats = importr('stats')

        r = robjects.r

        from rpy2.robjects import r, pandas2ri
        pandas2ri.activate()

        set_seed = r('set.seed')
        set_seed(1)

        ### overall

        myNewValue_LFQ_transposed = np.transpose(myNewValue_LFQ)
        #myPCAGrid = rrcov.PcaGrid(myNewValue_LFQ_transposed, 0.975, 2)

        #cutoffSD = myPCAGrid.slots['cutoff.sd']
        #cutoffOD = myPCAGrid.slots['cutoff.od']
        #myFlag = myPCAGrid.slots['flag']


        #### per group

        #print(myNewValue_LFQ_transposed.index.values)
        myLFQRows_unique = set(myNewValue_LFQ_transposed.index.values)

        myFlagsAll = []
        for name in myLFQRows_unique:

        #    print(name)
            myCurrentRows = myNewValue_LFQ_transposed.loc[name, :]

            myLFQCols_logical = value.columns.get_loc(name)
        #    print(myLFQCols_logical)
            myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]

            #print(myNewValue_LFQ_transposed.loc[name, :].index.values)
            #print(myNewValue_LFQ_transposed.loc[name, :].index.names)
            #print(myCurrentRows.shape)

            #myMaxOutliers = myCurrentRows.shape[0] - 2
            myMaxOutliers = myCurrentRows.shape[0] - 3

        #    print(myMaxOutliers)

            if myMaxOutliers > 0:

                myPCAGrid = rrcov.PcaGrid(myCurrentRows, 0.975, 2)

                #import plotly.express as px
                #fig = px.scatter(myPCAGrid, x=0, y=1)
                #fig.show()
                # , color=target, labels={'0': 'PC1 (' + str(round(varianceExplainedPC1 * 100, 2)) + '%)', '1': 'PC2 (' + str(round(varianceExplainedPC2 * 100, 2)) + '%)'})
                #fig.write_image(os.path.join(myBasepath, myOutputFile, myId, 'Plots', 'OutlierPlot_' + myId + '.jpg'))

                generatePlots.plotOutlierPlot(myPCAGrid, myId, myBasepath, myOutputFile)

                cutoffSD = myPCAGrid.slots['cutoff.sd']
                cutoffOD = myPCAGrid.slots['cutoff.od']
                myFlag = myPCAGrid.slots['flag']

                myOds = pd.Series(myPCAGrid.slots['od'],index=myLFQCols)
                mySds = pd.Series(myPCAGrid.slots['sd'], index=myLFQCols)

        #        print(myOds)

                myOdsDiff = myOds - myPCAGrid.slots['cutoff.od']
                mySdsDiff = mySds - myPCAGrid.slots['cutoff.sd']

        #        print(myOdsDiff)
        #        print(myOdsDiff[myOdsDiff > 0])

                myValues_all_max = myOdsDiff.combine(mySdsDiff, max)

                myValues_Outliers = myValues_all_max[myValues_all_max > 0]

        #        print(myValues_Outliers)

                #myOutliers_maxN = myValues_Outliers.nsmallest(n=2)

                #print(myOutliers_maxN)

                #print((myOds-myPCAGrid.slots['cutoff.od'])[base.which(base.as_vector(myOds) > myPCAGrid.slots['cutoff.od'])])

                #myDistancesToCutoffForOutliers = base.sort(base.c((myOds-myPCAGrid.slots['cutoff.od'])[base.which(base.as_vector(myOds) > myPCAGrid.slots['cutoff.od'])], (mySds-myPCAGrid.slots['cutoff.sd'])[base.which(base.as_vector(mySds) > myPCAGrid.slots['cutoff.sd'])]), decreasing = True)

                #print(myDistancesToCutoffForOutliers)
                #print(myDistancesToCutoffForOutliers.index)

                if myMaxOutliers >= len(myValues_Outliers):
                    myOutliersForRemoval = myValues_Outliers.index
                else:
                    myOutliersForRemoval = myValues_Outliers.nlargest(n=myMaxOutliers).index
        #        print(myOutliersForRemoval)

                myFlagsAll.extend(myOutliersForRemoval)

        #newValue = value.drop(columns=myFlagsAll, axis=1)

        #myLFQCols_logical = np.array([l.startswith('LFQ') for l in value.columns.values])
        #myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]

        #myZeroInds = [i for i, e in enumerate(myFlagsAll) if e == 0]
        #myZeroInds = [i for i, e in enumerate(myFlag) if e == 0]
        #myOutlierCols = myLFQCols[myZeroInds]

        #print(myFlagsAll)

        all_cols = set(range(0, len(value.columns)))
        keep_cols = all_cols - set(myFlagsAll)
        newValue = value.iloc[:, list(keep_cols)]

        return(newValue)


    def detectOutliers2(value, myId, myBasepath, myOutputFile):

        myNewValue_LFQ = value.filter(regex='LFQ intensity_')

        rrcovHD = importr('rrcovHD')
        base = importr('base')
        stats = importr('stats')
        rrcov = importr('rrcov')

        r = robjects.r

        from rpy2.robjects import r, pandas2ri
        pandas2ri.activate()

        set_seed = r('set.seed')
        set_seed(1)

        ### overall

        myNewValue_LFQ_transposed = np.transpose(myNewValue_LFQ)
        #myNewValue_LFQ_transposed['group'] = myNewValue_LFQ.columns
        myGroup = myNewValue_LFQ.columns

        print(myNewValue_LFQ_transposed)
        print(myGroup)

        obj = rrcovHD.OutlierPCOut(myNewValue_LFQ_transposed,myGroup)

        outlierDist = rrcov.getDistance(obj)
        outlierLabels = rrcov.getFlag(obj)

        ### for group in groups:
        ###     maxOutlier = len(Group) -3
        ###     plotOutlierPlot
        ###     if OutliersInGroup > 0 and < maxOutlier:
        ###         remove all outlier
        ###     else if OutliersInGroup > maxOutliers:
        ###         sort by Distance
        ###     remove largest maxOutlier outliers

        myLFQRows_unique = set(myNewValue_LFQ_transposed.index.values)

        myFlagsAll = []

        for name in myLFQRows_unique:

#            generatePlots.plotOutlierPlot2(obj, myId, myBasepath, myOutputFile)

            myCurrentRows = myNewValue_LFQ_transposed.loc[name, :]

            myLFQCols_logical = value.columns.get_loc(name)
            #    print(myLFQCols_logical)

            myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]

            myMaxOutliers = myCurrentRows.shape[0] - 3


            if myMaxOutliers > 0:

                print(outlierLabels)
                myValues_Outliers = np.where(np.array(outlierLabels)==0)
                print(myValues_Outliers)

                if myMaxOutliers >= len(myValues_Outliers):
                    #print(myCurrentRows)
                    #print(myLFQCols_logical)
                    #print(myLFQCols)
                    print("current:")
                    print(rrcovHD.getClassLabels(obj, 1))

#        myValues_Outliers = myValues_all_max[myValues_all_max > 0]

#        myOutliersForRemoval = myOutlierLabels.index
#                else:
#                    myOutliersForRemoval = myValues_Outliers.nlargest(n=myMaxOutliers).index
        #        print(myOutliersForRemoval)

        myFlagsAll.extend(myOutliersForRemoval)

        #newValue = value.drop(columns=myFlagsAll, axis=1)

        #myLFQCols_logical = np.array([l.startswith('LFQ') for l in value.columns.values])
        #myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]

        #myZeroInds = [i for i, e in enumerate(myFlagsAll) if e == 0]
        #myZeroInds = [i for i, e in enumerate(myFlag) if e == 0]
        #myOutlierCols = myLFQCols[myZeroInds]

        #print(myFlagsAll)

        all_cols = set(range(0, len(value.columns)))
        keep_cols = all_cols - set(myFlagsAll)
        newValue = value.iloc[:, list(keep_cols)]


        return(newValue)