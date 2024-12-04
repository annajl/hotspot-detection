
import numpy as np
import pandas as pd
import sklearn.cluster as sklc
from sklearn import manifold, decomposition
import networkx as nx

#supporting python scripts
import hot_mapper.utils as utils
#single linkage dendogram
from scipy.spatial import distance







def _build_cover_on_lens_function(data, lens_function, intervals, overlap):
    """Build a cover by dividing the lens into overlapping intervals and retrieve
    the samples contained in each interval"""

    #find the range of the lens function
    lens_min = np.amin(lens_function)
    lens_max = np.amax(lens_function)

    #calculate the size of each interval
    interval_length = (lens_max - lens_min) / (((intervals-1)  *  (1 - overlap)) + 1)

    #the number of samples in the dataset
    no_samples_in_data = np.array(range(data.shape[0]))

    #these are the fundamental properties of the intervals
    interval_samples = [False for i in range(intervals)]
    samples_in_interval = {}
    interval_sets = []

    for i in range (0, intervals):
        #the starting point of each interval is the start of the lens function
        ai = lens_min + (i * interval_length) * (1 - overlap)
        bi = ai + interval_length
        interval_sets.append([ai,bi])
        #for each point, assign it to an interval if the function value
        #of this point lies between interval start and end points
        points_in_set = ((ai <= lens_function) & (lens_function <= bi))
        samples = no_samples_in_data[points_in_set]
        samples_in_interval[i] = np.unique(samples)


    return samples_in_interval, interval_sets



def _minimum_samples_for_clustering_algorithm(clustering_algorithm):
    """ clustering algorithms can set the minimum number of clusters
    #and will throw up an error if there are too little samples in cluster
    # n_clusters = agglo,  min_cluster_size = hdbscan. Can add more if necessary"""
    algorithm_parameters = clustering_algorithm.get_params()
    algorithm_minimum = 2
    options = ["n_clusters", "min_cluster_size"]
    for setting in algorithm_parameters:
        if setting in options:
            algorithm_minimum = algorithm_parameters[setting]

    return algorithm_minimum



def _cluster_data_in_intervals(data, intervals, clustering_algorithm, samples_in_interval):
    """Perform clustering within each interval on the data in the original space.
    These clusters form nodes in the graph, and overlapping clusters are reperesented
    by edges. """

    cluster_samples_in_interval = {}
    n = 0 #n is the number of clusters at the start
    n_i = 0 # n_i is the number of clusters generated in that interval set
    min_samples_in_cluster = _minimum_samples_for_clustering_algorithm(clustering_algorithm)
    #for each interval, if there is more than two data points in that interval
    #Fit a clustering algorithm to those datapoints
    for i in range(0, intervals):
        #access the unique samples within each interval
        samples = samples_in_interval[i]
        points =  data[samples]

        #THIS SETS THE MIMINMUM NUMBER OF SAMPLES IN A NODE - we lose intervals if we do not incorporate nodes
        if len(samples) > min_samples_in_cluster:
            cl = clustering_algorithm.fit(points)
            c_i = [samples[np.where(cl.labels_ == label)] for label in set(cl.labels_)]
            cluster_samples_in_interval[i] = c_i
        n =  n + n_i
    return cluster_samples_in_interval


#
def _convert_sampleID_dict_to_matrix(data, ID_dictionary):
    mtx = pd.DataFrame(0, index=np.arange(len(data)), columns=ID_dictionary.keys())
    for node in ID_dictionary:
        node_columns = mtx[node]
        samples = ID_dictionary[node]
        node_columns[samples] = 1
    return mtx


def _build_cluster_index_labels(samples_in_clusters):
    #samples in clusters is composed of = [interval 0: [[samples in cluster 0][samples in cluster 1]]] etc..
    #return clusters_dict, composed of [interval 0 : [0,1], interval 1 : [2, 3] etc...]
    count = 0
    clusters_dict = {}

    for k,v in samples_in_clusters.items():
        start = count
        clusters_list = []
        for i, samples in enumerate(v):
            clusters_list.append(count)
            count += 1
        clusters_dict[k] = clusters_list
    return clusters_dict







class MapperGraph():
    """ The Mapper class builds a network graph from data. The workflow allows you to:
            1. Transform the dataset
            2. Construct a lens function to project the data to a lower dimension
            4. Build a cover of the lens from overlapping intervals
            3. Perform clustering within each interval
            5. Construct a graph from the clustering
            6. Visualise the graph using networkx

            """


    def __init__(self, data, lens_function, intervals, overlap, clustering_algorithm, text = True):

        self.data = data
        self.lens_function = lens_function
        self.intervals = intervals
        self.overlap = overlap
        self.clustering_algorithm = clustering_algorithm
        self.text = text

        if self.text == True:
            print("Initializing Mapper class...")



    def build_graph(self):
        """This involves two steps - build a cover on the lens function to divide it into overlapping intervals, then
        clustering in each interval on the original point cloud. The networkx graph is built by converting the clusters
        to nodes and edges are formed when two clusters have overlapping samples.  """


        if self.text == True:
            print("Build cover...")
        samples_in_intervals, interval_sets = _build_cover_on_lens_function(self.data, self.lens_function, self.intervals, self.overlap)

        if self.text == True:
            print("Build clusters...")
        samples_in_clusters = _cluster_data_in_intervals(self.data, self.intervals, self.clustering_algorithm, samples_in_intervals)

        #networkx graph class
        G = nx.Graph()
        node = 0
        node_dict = {}
        samples_in_node = {}

        if self.text == True:
            print("Build graph...")
        #construct a graph with a vertex for each cluster
        #access the samples in each cluster in each interval to build a node
        #samples in clusters is composed of = [interval 0: [[samples in cluster 0][samples in cluster 1]]] etc..
        node_count_in_intervals = {}
        nodes_in_intervals = _build_cluster_index_labels(samples_in_clusters)
        for i in samples_in_clusters:
            #create dictionary describing nodes in interval_sets
            node_count_in_intervals[i] = len(samples_in_clusters[i])
            #for each cluster in that interval
            for j in range(0, len(samples_in_clusters[i])):
                #define the node index values simultaneously
                samples_in_node[node] = samples_in_clusters[i][j]
                #samplecount_node[node] = len(np.unique(samples_in_clusters[i][j]))
                #add the node of samples to the graph
                G.add_node(node)
                #build the node dictionary to understand which intervals and clusters the nodes correspond to
                node_dict[(i,j)] = node
                node = node + 1



        #to build an edge check for overlapping samples between neighbouring nodes
        #define the key list to iterate through dict and check neighbouring intervals
        keyList=sorted(samples_in_clusters.keys())
        for i, l in enumerate(keyList[:-1]):
            l_neighbour = keyList[i+1]
            cluster_points = samples_in_clusters[l]
            cluster_points_neighbour = samples_in_clusters[l_neighbour]

            for cluster_label, main_cluster in enumerate(cluster_points):
                for cluster_label_neighbour, neighbour_cluster in enumerate(cluster_points_neighbour):
                    #if any points in the main cluster being check are in the neighbouring cluster
                    if any(point in main_cluster for point in neighbour_cluster):
                        G.add_edge(node_dict[(l,cluster_label)], node_dict[(l_neighbour,cluster_label_neighbour)])


        #convert samples_in_node from dict with node in keys and samples in values,
        #to df with nodes in columns and samples in rows and binary values indicating the presence of sample in node
        #def convert_index_dict_to_matrix():
        interval_clusters = _convert_sampleID_dict_to_matrix(self.data, samples_in_intervals)
        node_clusters = _convert_sampleID_dict_to_matrix(self.data, samples_in_node)




        self.graph = G
        self.samples_in_nodes = node_clusters
        self.node_count_in_intervals = node_count_in_intervals
        self.nodes_in_intervals = nodes_in_intervals
        self.samples_in_intervals = interval_clusters
        self.interval_sets = interval_sets

        return G
