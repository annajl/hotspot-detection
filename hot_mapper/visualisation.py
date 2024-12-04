
import networkx as nx
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import hot_mapper.utils as utils



plt.rcParams.update({'font.size': 36})

def draw_graph(mapper_graph, attribute_function, samples_in_nodes, hotspot_nodes = None,  style = 1, size = 1, labels = False, tick_labels = False, col_legend_title = "Legend", file_name = None, file_format = "png"):
    """Visualise the networkx graph.

    Parameters
    ----------

    style : [1],[2], default: ``1``
        Selects either fruchterman_reingold_layout from networkx (1) or kamada_kawai_layout (2) to structure the graph

    size : int, default: ``10``
        Size of node legends specifying the number of samples per nodes

    labels = boolean, default: ``False``
        Specifies whether to label nodes with number ids

    col_legend_title = str, default: ```Legend```
        Labels the attribute legend

        """

    graph = mapper_graph.copy()
    network_styles = {1: nx.fruchterman_reingold_layout(graph, seed=300),
                      2: nx.kamada_kawai_layout(graph),
                      3: nx.nx_pydot.graphviz_layout(graph)}

    #Plot figure with legend, specifying style
    cmap = mpl.cm.viridis
    fig = plt.figure(figsize=(12, 12), constrained_layout=True)
    pos =  network_styles[style]


    #colour the nodes by the attribute of choice
    #if attribute function provided as values for each sample, average per node
    if len(attribute_function) == samples_in_nodes.shape[0]:
        attribute_by_node = utils.colour_nodes_by_attribute(samples_in_nodes, attribute_function)

    elif len(attribute_function) == samples_in_nodes.shape[1]: 
        attribute_by_node = attribute_function

    else:
        print("attribute size wrong length: must be value for each sample or value for each node")
    #if attribute function provided as values per node, keep values 
    norm = mpl.colors.Normalize(vmin=min(attribute_by_node), vmax=max(attribute_by_node))
    colouring = cmap(norm((attribute_by_node)))

    #specify the number of samples in each node according to size attribute
    nsize = [np.sum(samples_in_nodes[node]) for node in samples_in_nodes]
    nodes = nx.draw_networkx_nodes(graph,
                              pos= pos,
                              node_color=colouring,
                              alpha=1,
                              node_size = [size * n for n in nsize])

    nodes.set_edgecolor('grey')
    nodes.set_linewidth(2)
    nx.draw_networkx_edges(graph,
                           pos = pos,
                           width= 3,
                           alpha = 0.5,
                           edge_color='dimgray')



    #if hotspot nodes are specified colour the outline of nodes and edges red
    if hotspot_nodes:
        H = graph.subgraph(hotspot_nodes)
        h_nsize = [nsize[i] for i in H.nodes]
        h_colouring = [colouring[i] for i in H.nodes]
        hnodes = nx.draw_networkx_nodes(H,
                                      pos= pos,
                                      node_color =  h_colouring,
                                      alpha=1,
                                      node_size = [size * n for n in h_nsize])
        hnodes.set_edgecolor('red')
        hnodes.set_linewidth(2.5)
        nx.draw_networkx_edges(H, pos, width= 8, alpha = 0.3, edge_color= "red")


    #position labels slighter offset to nodes
    if labels == True:
        pos_higher = {}
        off = 0.05  # offset on the y axis

        for k, v in pos.items():
            pos_higher[k] = (v[0], v[1])

        nx.draw_networkx_labels(graph,
                                pos_higher,
                                font_size=26,
                                font_color="black",
                                bbox = {"ec": "k", "fc": "white", "alpha": 0.6})


    #legend
    ticks = [min(attribute_by_node), max(attribute_by_node)]
    ax = plt.gca()
    ax.axis('off')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    cax = ax.inset_axes([1.1, 0.05, 0.04, 0.2])
    cbar = plt.colorbar(sm,cax=cax, ticks = ticks)
    if tick_labels == True:
        cbar.ax.set_yticklabels(['Low', 'High'])
    cbar.ax.set_title(col_legend_title, fontsize = 30, pad = 30)
    cbar.ax.tick_params(labelsize=26, pad =10, bottom = True)

    #define the size of the node labels and plot as non-existent dots
    #obtain the 10th, 50th, and 100th percentile
    legend_n = {int(round(np.percentile(list(nsize), 10),0)),
                int(round(np.percentile(list(nsize), 50),-1)),
                int(round(np.percentile(list(nsize), 100),-1))}

    for v in legend_n:
        plt.scatter([],[], s= (size * v), label='{}'.format(v))

    #get the legend handles and arrange to correct order
    handles,labels = ax.get_legend_handles_labels()
    handles, labels = zip(*[ (handles[i], labels[i]) for i in sorted(range(len(handles)), key=lambda k: list(map(int,labels))[k])] )
    ax.legend(handles, labels,title_fontsize = 30, fontsize=26, loc=2, bbox_to_anchor=(1, 1),fancybox=True, title = "no. of samples", labelspacing  = 1.6,frameon=False)

    #set colour of points
    leg = ax.get_legend()
    for i in [0,1,2]:
        try:
            leg.legendHandles[i].set_color('grey')
        except Exception:
            pass


    if file_name:
        plt.savefig(file_name, format = file_format)

    plt.show()
