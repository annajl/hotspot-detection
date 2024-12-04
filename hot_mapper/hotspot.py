# -*- coding: utf-8 -*-

import hot_mapper.utils as hmu
import hot_mapper.visualisation as hmv
import hot_mapper.utils as utils

import numpy as np
import networkx as nx
import pandas as pd
from statsmodels import robust
from itertools import compress
from scipy.spatial import distance
from scipy.cluster import hierarchy



import matplotlib as mpl
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
plt.rcParams.update({'font.size': 22})




#----------------------hotspot class--------------------------------#
class HotspotSearch:
    """ The hotspot class searches a collection of interconnected nodes
    obtained from the Mapper grah, identifying any clusters of nodes that
    present anomalous attribute filter_values
        """


    def __init__(self, mapper_graph, attribute_function, samples_in_nodes, text = True):

        #each subgraph needs to be assigned new labels as an new graph, previously contain
        #labels for original mapper. labels retained in nx attributes'_node'
        self.graph = mapper_graph
        self.samples_in_nodes = samples_in_nodes

        #if attribute function provided as values for each sample, average per node
        if len(attribute_function) == samples_in_nodes.shape[0]:
            self.node_attribute = utils.colour_nodes_by_attribute(samples_in_nodes, attribute_function)

        elif len(attribute_function) == samples_in_nodes.shape[1]: 
            self.node_attribute = attribute_function

        else:
            print("attribute size wrong length: must be value for each sample or value for each node")

    def _calculate_edge_weights(self):
        #assigns edge weights for every pair of edges in the original graph
        edge_weights = {}
        graph_copy = self.graph.copy()
        for e in graph_copy.edges():
            #find the absolute difference in nodes
            node_difference = abs(self.node_attribute[e[0]] - self.node_attribute[e[1]])
            graph_copy[e[0]][e[1]]['weight'] = node_difference #assign as edge weight
            edge_weights[(e[0],e[1])] = node_difference
        return edge_weights



    def _seperate_graph_connected_components(self, graph):
        """Hotspot detection is performed for each seperate component of the graph. These
        are converted to seperate networkx subgraphs"""

        connected_components = {i:graph.subgraph(c).copy() for i,c in enumerate(nx.connected_components(graph))}
        return connected_components



    def _sort_subgraph_edge_weights(self, subgraph, edge_weights):

        # #subset edges and weights to those in the subgraph
        subgraph_edge_weights = {}
        for k in subgraph.edges:
            try:
                subgraph_edge_weights[k] = edge_weights[k]

            #if edge order has swapped try reversing key
            except KeyError:
                k_rev = (k[1],k[0])
                subgraph_edge_weights[k_rev] = edge_weights[k_rev]


        #sort labels in order of ascending weight
        subgraph_edge_weights_sorted = {k: v for k, v in sorted(subgraph_edge_weights.items(), key=lambda item: item[1])}

        return subgraph_edge_weights_sorted




    def _identify_edge_cut_off(self, subgraph, edge_weights, plot_dendrogram = False):
        """ Define edge weights according to the difference in attribute value between nodes.
            Sort the edge weights, then identify the cut-off point to build the clusters - between the
            edge weights with the largest difference in values"""

        #only define cutoff if more than one node is present, as otherwise no edges exist

        if len(subgraph.nodes) == 1:
            cutoff = 1

        else:
            #return the edge weights found in the new subgraph
            edge_weights_sorted = self._sort_subgraph_edge_weights(subgraph, edge_weights)
            #empty matrix of node length x needed for linkage tree
            a = pd.DataFrame(1.0, index= subgraph.nodes, columns=subgraph.nodes)

            #construct matrix of nodes and edges
            for i,j in edge_weights_sorted.items():
                u = i[0]
                v = i[1]
                a.loc[u,v] = j
                a.loc[v,u] = j
                a.loc[u,u] = 0.0
                a.loc[v,v] = 0.0

            #construct a single linkage matrix of connected node
            dists = distance.squareform(a)
            Z = hierarchy.linkage(dists)
            
            #retain the threshold for edge cutoff and get the edge weight values from the dict
            #find the differences between all the values. We ignore the very last value
            edge_differences = list(edge_weights_sorted.values())
            edge_difference_distances = [j-i for i, j in zip(edge_differences[:-2], edge_differences[1:])]

            #find the maximum difference and the index
            try:
                m = max(edge_difference_distances)
                max_index = edge_difference_distances.index(m)
                cut_index = max_index

            except ValueError:
                m = 0
                cut_index = 0

            #index at the cutoff value + half way between the difference
            cutoff = edge_differences[cut_index]+ m*0.5





            if plot_dendrogram == True:
                #specify the nodes contained in this subgraph
                lab = list(subgraph.nodes())
                print(Z)
                print(lab)
                print(cutoff)
                graph_colours = {i:v for i,v in enumerate(self.node_attribute)}
                subgraph_colours = [graph_colours[i] for i in lab]

                #assign matching node colour values to dendrogram labels
                cmap = mpl.cm.viridis
                norm = mpl.colors.Normalize(vmin=min(self.node_attribute), vmax=max(self.node_attribute))
                colouring = cmap(norm((subgraph_colours)))
                vir_col = {str(v): colouring[i] for i, v in enumerate(lab) }

                dd = hierarchy.dendrogram(Z, 
                                        labels = list(lab),
                                        color_threshold = cutoff,
                                        above_threshold_color = "grey",
                                        leaf_font_size = 14)

                #add the matching colours to the tick labels
                ax = plt.gca()
                xlbls = ax.get_xmajorticklabels()
                for lbl in xlbls:
                    lbl.set_color(vir_col[lbl.get_text()])
                plt.show()

        return cutoff


    def _identify_attribute_clusters_below_cutoff(self, subgraph, edge_weights, cutoff):
        """Function seperates the graph into clusters according to the edge
        weight cutoff point. Any edges above this threshold are removed, leaving
        a graph seperated into groups of similar node values"""

        #return the edge weights found in the new subgraph
        edge_weights_sorted = self._sort_subgraph_edge_weights(subgraph, edge_weights)

        #make a copy of the original subgraph
        graph_copy = subgraph.copy()

        #identify the any edges above the cut-off
        remove_edges = [i for i in edge_weights_sorted if edge_weights_sorted[i] > cutoff]

        #remove these values from the graph
        for e1,e2 in remove_edges:
            graph_copy.remove_edge(e1,e2)

        #map between the old and new labels
        communities = self._seperate_graph_connected_components(graph_copy)
        cluster_nodes = [list(cluster.nodes) for cluster in list(communities.values())]
        
        return(cluster_nodes)


    def _find_attribute_value_of_node_cluster(self, nodes):
        """Function obtains the mean filter value for each community in the Mapper graph"""

        #for each node obtain the list of indexes
        node_attribute_values = [self.node_attribute[i] for n,i in enumerate(nodes)]
        return round(np.mean(node_attribute_values),3)


    def _find_no_samples_in_nodes(self, nodes):
        """Function finds the sample size for specified nodes in the Mapper graph"""
        return np.sum(self.samples_in_nodes[nodes].max(axis=1))


    def _cluster_classification(self, subgraph, community_clusters, attribute_threshold, min_sample_size, attribute_extreme = "either"):
        ## Set up the hotspot dictionary containing information for each subgraph ##
        hotspot = {}
        hotspot_class = [True]*len(community_clusters)

        #find all the nodes in this graph community
        component = subgraph.nodes

        for i,cluster in enumerate(community_clusters):

            #find neighbour nodes - the other remaining nodes in the component
        #    print(f"cluster {list(cluster)}")
            neighbour = list(np.setdiff1d(component, cluster))

            #find mean attribute values of cluster and neighbour
            cluster_mean_attribute = self._find_attribute_value_of_node_cluster(cluster)
            neighbour_mean_attribute = self._find_attribute_value_of_node_cluster(neighbour)
        #    print(f"mean attribute: \nhotspot {cluster_mean_attribute} neighbours {neighbour_mean_attribute}\n ")
            #find sample size of cluster and neighbour
            cluster_sample_size = self._find_no_samples_in_nodes(cluster)
            neighbour_sample_size = self._find_no_samples_in_nodes(neighbour)
        #    print(f"mean size: \nhotspot {cluster_sample_size} neighbours {neighbour_sample_size}\n ")

            # CHECK 1 - Size of samples in the cluster is sufficiently larger than threshold
            if cluster_sample_size < min_sample_size:
            #    print("Size of samples in the cluster is sufficiently lower than threshold")
                hotspot_class[i] = False


            # CHECK 2 - Size of samples in the cluster is smaller than neighbourhood
            mad_threshold = robust.mad([cluster_sample_size,neighbour_sample_size], c=1)
            if (neighbour_sample_size - cluster_sample_size) < mad_threshold:
            #    print("# Size of samples in the cluster is larger than neighbourhood")
                hotspot_class[i] = False


            #CHECK 3 - Check the attribute difference between cluster and neighbourhood is large enough
            #conditional on parameter describing extreme of attribute function to be investigated
            low_check = (neighbour_mean_attribute < cluster_mean_attribute)
            high_check = (neighbour_mean_attribute > cluster_mean_attribute)
            att_check = abs(neighbour_mean_attribute - cluster_mean_attribute) < attribute_threshold

            extreme_options = {"lower": any([att_check, low_check]),
                                "higher": any([att_check, high_check]),
                                "either": att_check}

            if extreme_options[attribute_extreme] == True:
            #    print(" attribute difference between cluster and neighbourhood is not large enough")
                hotspot_class[i] = False


        return list(compress(community_clusters, hotspot_class))




    def search_graph(self, attribute_threshold, min_sample_size, attribute_extreme = "either", plot_dendrogram = False):
        #calculate the weight of the edges as the difference in attribute between the nodes
        edge_weights = self._calculate_edge_weights()
        
        #identify subgraphs of connected components in networkx graph
        subgraphs = self._seperate_graph_connected_components(self.graph)
        

        #for each subgraph identify the cut-off point between edges
        subgraph_cutoffs = [self._identify_edge_cut_off(subgraphs[i], edge_weights, plot_dendrogram = plot_dendrogram) for i in subgraphs]

        #identify the community clusters in the graph that lie below the attribute cut-off
        community_cluster_nodes = [self._identify_attribute_clusters_below_cutoff(subgraphs[i], edge_weights, subgraph_cutoffs[i]) for i in subgraphs ]

        #classify each commmunity cluster in the graph as a hotspot or non-hotspot
        hotspot_clusters = [self._cluster_classification(subgraphs[i], community_cluster_nodes[i], attribute_threshold, min_sample_size, attribute_extreme) for i in subgraphs]

        #return flat list of hotspot nodes from all clusters in all components
        hotspot_nodes = [nodes for component in hotspot_clusters for nodes in component ]
        self.hotspots = hotspot_nodes
        return hotspot_nodes



    def visualise_hotspots_in_graph(self, size = 10, style = 1, labels = False):
        #draw graph highlighting all hotspot nodes that may be present in each components
        #draw as seperate graphs
        for hotspot_nodes in self.hotspots:
            #for hotspot_nodes in component:
            print(hotspot_nodes)
            hmv.draw_graph(mapper_graph = self.graph,
                          attribute_function = self.node_attribute,
                          samples_in_nodes = self.samples_in_nodes,
                          hotspot_nodes = hotspot_nodes,
                          size = size,
                          style = style,
                          labels = labels)
