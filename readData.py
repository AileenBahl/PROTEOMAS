import os
import pandas as pd
import re


class readData:


    def getAllFilenames(myBasepath):

        myFilePaths = []

        for dirpath, dirnames, files in os.walk(os.path.join(myBasepath,'Input_MWCNT_tmp')):
            #print(dirpath)
            for file_name in files:
                myFilePath = os.path.join(dirpath,file_name)
                myFilePaths.append(myFilePath)

        return(myFilePaths)


    def extractProjectIds(myFilenames):

        myProjectIds = []

        for file in myFilenames:

            if 'proteinGroup' in file:

                myProjectId = re.search('proteinGroups' + '_' + '(.+?)' + '.txt', file).group(1)
                #myProjectId = re.search('proteinGroups' + '_' + '(.+?)' + '.xlsx', file).group(1)
                myProjectIds.append(myProjectId)

        return(myProjectIds)


    def readFiles(file,myFileType,myFileExtension):
        # read all ConditionAssignment files
        if myFileExtension == '.txt':
            myFile = pd.read_csv(file,delimiter="\t")
        elif myFileExtension == '.xlsx':
            myFile = pd.read_excel(file)

        return myFile
