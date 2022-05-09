import matplotlib.pyplot as plt
import numpy as np
import os
from bioinfokit import analys, visuz

class generatePlots:

    def plotHistogram(value, myId, myStep, myBasepath):

        myLFQCols_logical = np.array([l.startswith('LFQ') for l in value.columns.values])
        #print(myLFQCols_logical)
        myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]
        #print(myLFQCols)

        for i in myLFQCols:

            myValues = value.iloc[:, i]
            print(value.columns[i])
            myPlotName = value.columns[i]+"_"+str(i)

            #plt.style.use('_mpl-gallery')

            # plot:
            fig, ax = plt.subplots()

            ax.hist(myValues, bins=20)
            plt.xlabel('LFQ values')
            plt.ylabel('Count')
            plt.title('Histogram_' + myPlotName)
            plt.xlim(myValues.min(),myValues.max())
            plt.grid(False)
            plt.tight_layout()

            plt.savefig(os.path.join(myBasepath, 'Output', myId, 'Plots', 'Histogram_'+myStep+myId+myPlotName+'.jpg'))
            plt.close()

            #plt.show


    def plotBoxplot(value, myId, myStep, myBasepath):

        myValues = value.filter(regex='LFQ')

        myLFQCols_logical = np.array([l.startswith('LFQ') for l in value.columns.values])
        myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]
        myPlotNames = [m + '_' + str(n) for m, n in zip(value.columns[myLFQCols], list(range(len(myLFQCols))))]


        #plt.figure(figsize=(1, 1))
        plt.boxplot(myValues, positions=np.arange(2, len(myLFQCols)*2+1, 2).tolist(), widths=1.5)
        # Why doesn't it work for filtered data

        plt.xticks(ticks=np.arange(2, len(myLFQCols)*2+1, 2).tolist(), labels=myPlotNames, rotation='vertical')

        plt.xlabel('Sample')
        plt.ylabel('LFQ values')
        plt.title('Boxplot')
        plt.grid(False)
        plt.tight_layout()

        plt.savefig(os.path.join(myBasepath, 'Output', myId, 'Plots', 'Boxplot_' + myStep + myId + '.jpg'))
        plt.close()


    def plotVolcano(value, myId, key, myBasepath):

        print(value)
        print(myId)
        print(key)

        # only if significant
        #if abs(value['LogRatio'].any()) >= 0.5849625 and abs(value['Pvalue'] <= 0.05)):
        visuz.GeneExpression.volcano(df=value, lfc='LogRatio', pv='Pvalue', sign_line=True, lfc_thr=[0.5849625,0.5849625], figtype='jpg')
        # change path


    def plotPCA(value, myId, myStep, myBasepath):

        myValues = value.filter(regex='LFQ')

        myLFQCols_logical = np.array([l.startswith('LFQ') for l in value.columns.values])
        myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]
        myPlotNames = [m + '_' + str(n) for m, n in zip(value.columns[myLFQCols], list(range(len(myLFQCols))))]

        sklearn.decomposition.PCA().fit_transfrom()

        bioinfokit.visuz.cluster.biplot



    # def plotHeatmap

    # def plotOutlierPlot

    # def plotClusterplotForKNN