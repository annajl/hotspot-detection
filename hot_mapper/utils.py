# -*- coding: utf-8 -*-

import networkx as nx
from sklearn.metrics import  silhouette_score
import numpy as np
from statistics import mean
from itertools import chain
import pandas as pd
import matplotlib as mpl

def sample_index_in_nodes(node_index_dataframe, node_list):
    return node_index_dataframe[node_index_dataframe[node_list].max(axis=1) == 1].index

def colour_nodes_by_attribute(node_index_dataframe, attribute, norm = False):
    """for each node, create a list of the y values for each index"""
    if norm == True:
        #normalise y values to range
        attribute = range01(np.array(attribute))


    #the node values are averaged over all patients contained in each node
    node_values = []
    for node in node_index_dataframe:
        index = list(node_index_dataframe[node].loc[node_index_dataframe[node] == 1].index)
        node_y = [attribute[i] for i in index]
        node_values.append(np.mean(node_y))

    return node_values

#

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r


def range01(x):
    """Return the scaled transofmration of an array between 0 and 1"""
    r01 = np.array([(x- min(x))/ (max(x)- min(x))])
    return r01


def generate_artificial_hotspot(X, radius, random_seed = 42):
    #random sample from X
    np.random.seed(random_seed)
    core_sample = X[np.random.randint(X.shape[0], size=1), :]

    #generate y labels
    y = [0] * X.shape[0]
    for i in X:
        dist = np.linalg.norm(i - core_sample)
        if dist < radius:
            l = np.where(X==i)
            index = l[0][1]
            y[index] = 1
    return y





def norm(power,mag,matrix):
    """Computing euclidean norm of the matrix"""
    m = np.abs(matrix)
    m = np.power(matrix,power)
    m = np.sum(m,axis=1)
    m = np.power(m,mag/power)
    return m


def flatten(listOfLists):
    "Convert a list of lists to a single list"
    return chain.from_iterable(listOfLists)



def check_hotspot_nodes_in_intervals(mapper, hotspot_nodes):
    # map  nodes to intervals
    nodes_to_intervals = {}
    count = 0
    while count < mapper.samples_in_nodes.shape[1]:
        for k,v in mapper.node_count_in_intervals.items():
            for i in range(0,v):
                if k in nodes_to_intervals:
                    nodes_to_intervals[k].append(count)
                else:
                    nodes_to_intervals[k] = [count]
                count =  count + 1


    # count the number of hotspot nodes and non-hotspot nodes in an interval
    class_count_df = []
    for k,v in nodes_to_intervals.items():
        class_count_h = 0
        class_count_nh = 0
        for i in v:
            if i in hotspot_nodes:
                class_count_h += 1
            else:
                class_count_nh += 1
        class_count_df.append([class_count_nh, class_count_h])

    #plot hotspot nodes in intervals
    class_df = pd.DataFrame(class_count_df, columns = ["Non-Hotspot", "Hotspot"])
    mpl.style.use('ggplot')
    ax = class_df.plot.bar(figsize=(10,5),
                            xlabel='Intervals',
                            ylabel='Node count',
                            title  = 'Number of nodes in intervals',
                            color = ["indigo","yellow"],
                            stacked = True)
    ax.legend(loc=2)
