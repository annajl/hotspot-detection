import hot_mapper.mapper as mapper_algorithm
import hot_mapper.utils as utils
import hot_mapper.hotspot as hotspot_algorithm
import hot_mapper.random_lens as linear_lens_combination
import hot_mapper.visualisation as mapper_plot
from itertools import product

class Search():
    """This class searches across the Mapper parameters to identify the parameters that build a graph which contains a hotspot

    Parameters
    ----------

    X : pandas dataframe
        Dataset from which to obtain lens and build mapper graphs

    runs: into
        How many times to run the search for a lens
            """

    def __init__(self, X, runs = 1):
        self.X = X
        self.runs = runs
        self.parameters = {}
        self.parameter_lens = []
        self.parameter_samples = {}

    def build_graphs(self, parameters, visualise = False):
        """Search through the parameter options and build mapper graphs

        Parameters
        ----------

        parameters : dictionary
            Dictionary of parameter options.
            Must include lens / clustering algorithm / interval list / overlap list / epsilon / minimum sample size / attribute function / hotspot extremity
        """

        #Runs = lens space
        count = 0
        print('Testing Testing')
        signficance = False
        print("Building parameters and searching for hotspots")
        while count < self.runs:
            if signficance == False:
                #generate a random lens from features
                if parameters["predefined_lens"] is None:
                    random_lens = linear_lens_combination.Lens(self.X, nonzero_features = parameters["non_zero_lens_features"])
                else: 
                    random_lens = parameters["predefined_lens"]
                    
                #build a grid of interval and overlap combinations
                io_list = list(product(parameters["interval_list"], parameters["overlap_list"]))

                #for each grid combination of the interval & overlap, search lens for hotspot
                for i in io_list:
                    #define intervals and overlap
                    i_param = i[0]
                    o_param = i[1]

                    # Run a clustering algorithm and build the graph
                    mapper = mapper_algorithm.MapperGraph(data = self.X,
                                                            lens_function = random_lens["lens"],
                                                            intervals = i_param,
                                                            overlap = o_param,
                                                            clustering_algorithm = parameters["clustering_algorithm"],
                                                            text = False)

                    #build the graph with edges and nodes
                    mapper.build_graph()

                    #visualise graph
                    if visualise == True:
                        mapper_plot.draw_graph(mapper_graph = mapper.graph,
                                        attribute_function = parameters["attribute_function"],
                                        samples_in_nodes = mapper.samples_in_nodes,
                                        size = 5,
                                        style = 2,
                                        labels = False)

                    #run hotspot detection
                    hotspot_search = hotspot_algorithm.HotspotSearch(mapper_graph = mapper.graph,
                                                 attribute_function = parameters["attribute_function"],
                                                 samples_in_nodes = mapper.samples_in_nodes)

                    hotspots = hotspot_search.search_graph(attribute_threshold = parameters["epsilon"],
                                                                min_sample_size = parameters["min_samples"],
                                                                attribute_extreme = parameters["extreme"])

                    #return list of samples in each hotspot found
                    sample_list = []
                    for n in hotspots:
                        sample_list.append(utils.sample_index_in_nodes(mapper.samples_in_nodes, n))


                    #if hotspot present, save properties
                    if any(hotspots):
                        
                        self.parameters[(i_param,o_param)] = hotspots # list of hotspots
                        self.parameter_lens = {"weights": random_lens["weights"],
                                                "feature_list": random_lens["feature_list"]}
                        self.parameter_samples[(i_param,o_param)] = sample_list
                        significance = True

                #if hotspots exist in the filter function search
                if self.parameters:
                    print("\nHotspots search successful")
                    return
                else:
                    count += 1
                    print(F"\nCompleted {count} searches")
