import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from bioinfokit import analys, visuz
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import plotly.express as px
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
#import dash_bio
import plotly.graph_objects as go
import gc

matplotlib.use("Agg")

class generatePlots:

    def plotHistogram(value, myId, myStep, myBasepath, myOutputFile):

        if myStep == "beforeNormalization" or myStep == "vsn" or myStep == "quantileNorm":
            myLFQCols_logical = np.array([l.startswith('LFQ intensity ') for l in value.columns.values])
        else:
            myLFQCols_logical = np.array([l.startswith('LFQ intensity_') for l in value.columns.values])

        #print(myLFQCols_logical)
        myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]
        #print(myLFQCols)

        for i in myLFQCols:

            myValues = value.iloc[:, i]
            print(value.columns[i])
            print(myValues)
            myPlotName = value.columns[i]+"_"+str(i)

            #plt.style.use('_mpl-gallery')

            # plot:
            if not all(myValues.isnull()):

                fig, ax = plt.subplots()

                ax.hist(myValues, bins=20)
                plt.xlabel('LFQ intensity values')
                plt.ylabel('Count')
                plt.title('Histogram_' + myPlotName)
                plt.xlim(myValues.min(),myValues.max())
                plt.grid(False)
                plt.tight_layout()

                plt.savefig(os.path.join(myBasepath, myOutputFile, myId, 'Plots', 'Histogram_'+myStep+myId+myPlotName+'.jpg'))
                plt.close()


    def plotBoxplot(value, myId, myStep, myBasepath, myOutputFile):

        print(value.shape)
        #print(value)

        if myStep == "beforeNormalization" or myStep == "vsn" or myStep == "quantileNorm":
            print("Creating boxplot")
            myValues = value.filter(regex='LFQ intensity ')
            myLFQCols_logical = np.array([l.startswith('LFQ intensity ') for l in value.columns.values])
        else:
            myValues = value.filter(regex='LFQ intensity_')
            myLFQCols_logical = np.array([l.startswith('LFQ intensity_') for l in value.columns.values])

        myLFQCols = np.arange(0, value.shape[1])[myLFQCols_logical]
        myPlotNames = [m + '_' + str(n) for m, n in zip(value.columns[myLFQCols], list(range(len(myLFQCols))))]

        print("Boxplot")
        print(myValues)
        print(np.arange(2, len(myLFQCols)*2+1, 2).tolist())

        #plt.figure(figsize=(1, 1))

        plt.boxplot(myValues, positions=np.arange(2, len(myLFQCols)*2+1, 2).tolist(), widths=1.5)
        plt.xticks(ticks=np.arange(2, len(myLFQCols)*2+1, 2).tolist(), labels=myPlotNames, rotation='vertical')

        plt.xlabel('Sample')
        plt.ylabel('LFQ intensity values')
        plt.title('Boxplot')
        plt.grid(False)
        plt.tight_layout()

        plt.savefig(os.path.join(myBasepath, myOutputFile, myId, 'Plots', 'Boxplot_' + myStep + '_' + myId + '.jpg'))
        plt.close()


    def plotVolcano(value, myId, key, myBasepath, myOutputFile):

        print(value)
        print(myId)
        print(key)

        # only if significant
        #if abs(value['LogRatio'].any()) >= 0.5849625 and abs(value['Pvalue'] <= 0.05)):
        #visuz.GeneExpression.volcano(df=value, lfc='LogRatio', pv='Pvalue', sign_line=True, lfc_thr=[0.5849625,0.5849625], figname= os.path.join(myBasepath, myOutputFile, myId, 'Plots', 'VolcanoPlot_' + myId+'_'+key),figtype='jpg')

        #fig = dash_bio.VolcanoPlot(dataframe=value)
        #fig.write_image(os.path.join(myBasepath, myOutputFile, myId, 'Plots', 'VolcanoPot_' + myId + '.jpg'))

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=value['LogRatio'],y=-np.log10(value['Pvalue']),mode='markers'))
        fig.add_vline(x=np.log2(1.5), line_color="red")
        fig.add_vline(x=-np.log2(1.5), line_color="red")
        fig.add_hline(y=-np.log10(0.05), line_color="red")
        fig.update_layout(xaxis_title="log2 ratio", yaxis_title="-log10(p-value)")
        fig.write_image(os.path.join(myBasepath, myOutputFile, myId, 'Plots', 'VolcanoPlot_' + myId + '.jpg'))

        ### separate names for multiple conditions

    def plotPCA(value, myId, myStep, myBasepath, myOutputFile):

        if myStep == "beforeLog2" or myStep == "vsn" or myStep == "quantileNorm":
            myValues = value.filter(regex='LFQ intensity')
            print(myValues)
        else:
            myValues = value.filter(regex='LFQ intensity_')
            print(myValues)

        myValues_transposed = myValues.transpose()

        target = myValues_transposed.index.to_numpy()

        myValues_transposed_st = StandardScaler().fit_transform(myValues_transposed)

        pca=PCA(n_components=2)
        components = pca.fit_transform(myValues_transposed_st)

        varianceExplainedPC1 = pca.explained_variance_ratio_[0]
        varianceExplainedPC2 = pca.explained_variance_ratio_[1]

        fig = px.scatter(components, x=0, y=1, color=target, labels={'0': 'PC1 ('+str(round(varianceExplainedPC1*100,2))+'%)', '1': 'PC2 ('+str(round(varianceExplainedPC2*100,2))+'%)'})
        fig.write_image(os.path.join(myBasepath, myOutputFile, myId, 'Plots', 'PCABiplot_' + myStep + myId + '.jpg'))


    # def plotHeatmap

    def plotHeatmap(value, myId, myStep, myBasepath, myOutputFile):

        import seaborn as sns

        myValues = value.filter(regex='LFQ intensity_')
        myValues.columns = myValues.columns

        sns.clustermap(myValues, cmap="coolwarm", yticklabels=False)
        plt.savefig(os.path.join(myBasepath, myOutputFile, myId, 'Plots', 'Heatmap_' + myStep + myId + '.jpg'))
        plt.close()


    def plotOutlierPlot(myPCAGrid, myId, myBasepath, myOutputFile):

        r = robjects.r

        r.jpeg(os.path.join(myBasepath, myOutputFile, myId, 'Plots', 'OutlierPlot_' + myId + '.jpg'))
        #     res=500,
        #     width=2000, height=2000)
        # par(cex=0.6)
        r.plot(myPCAGrid)
        r('dev.off()')


    def plotClusterplotForKNN(myNewValue_LFQ, k, myId, myBasepath, myOutputFile):

        factoextra = importr('factoextra')
        stats = importr('stats')

        km = stats.kmeans(myNewValue_LFQ, k, nstart=30)

        print("km:")
        print(km)

        r = robjects.r

        r.jpeg(os.path.join(myBasepath, myOutputFile, myId, 'Plots', 'ClusterPlotKNN_' + myId + '.jpg'))
        #     res=500,
        #     width=2000, height=2000)
        # par(cex=0.6)
        #print(km)
        #print(myNewValue_LFQ)
        factoextra.fviz_cluster(km, myNewValue_LFQ)
        r('dev.off()')