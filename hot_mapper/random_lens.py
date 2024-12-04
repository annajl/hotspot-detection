# -*- coding: utf-8 -*-
import numpy as np


def Lens(data, nonzero_features = None, weights = None, feature_list = None, weight_range = [-1,1]):
    """Return a linear combination of a subset of features for each vector"""
    
    #define the number of features samples used in the data
    #do not perform feature selection unless specified
    total_samples, total_features = data.shape

    #randomly select the subset of features and their corresponding weights
    if weights is None and feature_list is None:

        #create an index for the subset of features sampled in the lens function
        rng = np.random.default_rng()
        feature_list = rng.choice(total_features, nonzero_features, replace=False)
        
        #create a list of corresponding weights for each feature
        weights = np.random.uniform(low=weight_range[0], high=weight_range[1], size=len(feature_list))
    
    if nonzero_features is None:
        nonzero_features = total_features
        
    
    #for each sample in the dataset calculate the lens function 
    lens = np.zeros(total_samples)
    
    #for each vector in the dataset
    for v in range(total_samples):
        x = data[v,:]
    
        for k in range(nonzero_features):
            feature_i = feature_list[k]
            #find the sum of all features x random
            lens[v] += x[feature_i] * weights[k]
            
    lens_settings = {"lens": lens,
                    "weights": weights,
                    "feature_list": feature_list}
    
    return lens_settings