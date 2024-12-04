# -*- coding: utf-8 -*-
"""

A module to perform Disease Specific Genomic Analysis (DSGA). 

A healthy state model is built from normal tissue data. 

Diseased tissue data is transformed into two components - the normal component and the diseased component. 

"""

import numpy as np
from sklearn.decomposition import PCA
import scipy.signal as ss
from scipy.stats.stats import pearsonr
import pandas as pd




def flat_construction(df_N):
    """Perform FLAT construction by constructing a linear model fit of the normal tumour vector genes as rows and samples as columns"""

    df_normal = df_N.to_numpy()
    normal_T = df_normal.transpose()
    normal_FLAT = np.empty(normal_T.shape)

    #for each gene in the normal vector
    for i in range(0,normal_T.shape[1]):
        n_i = normal_T[:,i] #take each gene seperately
        normal_i = np.delete(normal_T, i, 1) #remove it from the dataframe
        x = np.linalg.lstsq(normal_i, n_i, rcond=None) #construct linear model
        b = normal_i@x[0] #find fit of model
        normal_FLAT[:,i] = b #this is the new normal vector value

    return normal_FLAT


def wold_invariant(normal_FLAT):
    """Compute and plot the Wold invariant of PCA. Identify the number of principal components for the healthy state model"""

    features_no = normal_FLAT.shape[1] #number of samples
    pca = PCA(n_components=features_no)
    principalComponents = pca.fit_transform(normal_FLAT)

    sin_val = pca.singular_values_ #singular values
    sin_val_sqr = np.square(sin_val) #squared singular values

    R = normal_FLAT.shape[1]
    n = normal_FLAT.shape[0]

    #calculate list holding Wold invariant values
    wold = [(sin_val_sqr[l]/(np.sum(sin_val_sqr[l+1:])))*(((n-l-1)*(R-l))/(n+R-2*l)) for l in range(0,R)]

    #find the most prominent peak in Wold value
    peaks, _ = ss.find_peaks(wold)
    prominences = ss.peak_prominences(wold, peaks)[0]

    #find the index (peak) where the prominence spikes
    spike = np.amax(prominences)
    spike_index = peaks[np.where(prominences == spike)]

    return (principalComponents, spike_index)


def HSM(df_T, principalComponents, spike_index):
    """Choose the number of Wold components to build the Healthy State Model"""

    df_tumour = df_T.to_numpy()
    tumour_T = df_tumour.transpose()
    HSM = principalComponents[:,:int(spike_index)]

    #fit tumour vector to healthy state model
    x = np.linalg.lstsq(HSM, tumour_T, rcond=None)

    NcT = HSM@x[0] #healthy component
    DcT = tumour_T - NcT #diseased component

    return DcT


def threshold_coord(DcT):
    """Threshold data coordinates (genes, proteins, etc.) so that only the genes that show a significant deviation from the healthy state are retained """

    #find the 5th and 95th quantile of each gene 
    #take the absolute value
    q_df = pd.DataFrame({"Q5":DcT.quantile(q=0.05, axis=1),
                         "Q95":DcT.quantile(q=0.95, axis=1)})

    q_abs = q_df[["Q5", "Q95"]].max(axis=1)

    #list of genes that pass the 85th and 98th percentile
    q85 = q_abs.quantile(q=0.85) #relaxed
    q98 = q_abs.quantile(q=0.98) #stringent

    relaxed = DcT[q_abs > q85]
    stringent = DcT[q_abs > q98]


    #empty dataframe of patients
    Dc_mat = pd.DataFrame(columns = DcT.columns)

    #if the gene passess the 85th threshold
    for i in range(0 , len(relaxed)):
        correlation_QR = []
        gene_r = relaxed.iloc[i,:]

        #find the correlation between that gene and the stringent genes
        for j in range(0 , len(stringent)):
            gene_s = stringent.iloc[j,:]

            corr, _ = pearsonr(gene_r, gene_s)
            correlation_QR.append(corr)

        #if r>0.6 with over 3 stringent genes
        sig_sum = sum(val > 0.6 for val in correlation_QR)

        #retain that gene
        if sig_sum > 0.6:
            Dc_mat = Dc_mat.append(gene_r)

    #return matrix so patients are rows and genes are columns
    Dc_mat_T = Dc_mat.T
    return Dc_mat_T


def DSGA(df_normal, df_tumour, threshold = True):
    """Perform disease-specific genomic analysis on a tumour dataset, referencing against a healthy tissue dataset"""

    print(str(df_tumour.shape[1]) + " co-ordinates as input")
    #obtain flat construction of normal genes
    df_normal_flat = flat_construction(df_normal)

    #calculate the number Wold principal components
    #and the value at which they spike
    principalComponents, spike_index = wold_invariant(df_normal_flat)

    #construct a healthy state model using the tumour data
    DcT = HSM(df_tumour, principalComponents, spike_index)

    #convert DcT to a pandas dataframe with the genes and columns kept
    DcT_df =  pd.DataFrame(data=DcT,
                           index=df_tumour.columns,
                           columns=df_tumour.index)
    
    Dc_mat = DcT_df

    print("\n")

    if threshold == True:
        #threshold the data coordinates
        Dc_mat = threshold_coord(DcT_df)
        print(str(Dc_mat.shape[1]) + " co-ordinates are retained")
    else:
        print(str(Dc_mat.shape[0]) + " co-ordinates are retained")

    return(Dc_mat)
