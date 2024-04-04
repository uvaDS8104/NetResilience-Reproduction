import numpy as np
import graph_tool.all as gt
import networkx as nx
import glob
from csv import reader as csvreader
from csv import list_dialects


def _from_csv(csvdir, create_using=nx.Graph):
    ## Build g from edgelist

    # add edges
    edgelist = []
    with open(f'{csvdir}edges.csv', 'r') as csvfile:
        edgereader = csvreader(csvfile)
        header = next(edgereader)
        for row in edgereader:
            source = int(row[0])
            target = int(row[1])
            edgelist.append((source, target))

    g = nx.from_edgelist(edgelist, create_using=create_using)
    # print('original')
    # print(nx.number_of_nodes(g))
    # print(nx.number_of_edges(g))

    # convert to undirected graph and remove self loops
    # print('after to undirected and remove self loops')
    g = g.to_undirected()
    g.remove_edges_from(nx.selfloop_edges(g))
    # print(nx.number_of_nodes(g))
    # print(nx.number_of_edges(g))

    return g


def nx_gen_save_synth_graphs(n_graphs):
    N = 10000
    M = 20000
    m = 2

    for i in range(n_graphs):
        g_er = nx.gnm_random_graph(N, M)
        g_pa = nx.barabasi_albert_graph(N, m)

        # remove self loops, if any
        g_er.remove_edges_from(nx.selfloop_edges(g_er))
        g_pa.remove_edges_from(nx.selfloop_edges(g_pa))

        nx.write_gml(g_er, f'data/ER/{i:02d}.gml')
        nx.write_gml(g_pa, f'data/PA/{i:02d}.gml')


def nx_save_real_graphs():
    internet = 'data/internet/'
    www = 'data/worldwideweb/'

    g1 = _from_csv(internet)
    g2 = _from_csv(www, create_using=nx.DiGraph)

    nx.write_gml(g1, 'data/internet/internet.gml')
    nx.write_gml(g2, 'data/worldwideweb/worldwideweb.gml')


def nx_load_synth_graphs():
    root_dir_er = 'data/ER/'
    root_dir_pa = 'data/PA/'
    ers = glob.glob('*.gml', root_dir=root_dir_er)
    pas = glob.glob('*.gml', root_dir=root_dir_pa)
    ers.sort()
    pas.sort()

    er_list = []
    pa_list = []

    # print(ers)
    # print(pas)

    for i in range(len(ers)):
        g_er = nx.read_gml(f'data/ER/{i:02d}.gml')
        g_pa = nx.read_gml(f'data/PA/{i:02d}.gml')

        er_list.append(g_er)
        pa_list.append(g_pa)

    # build gdict
    gdict = {}
    gdict['N'] = 10000
    gdict['n_graphs'] = len(ers)
    gdict['ER'] = er_list
    gdict['PA'] = pa_list

    return gdict


def nx_load_real_graphs():
    g1 = nx.read_gml('data/internet/internet.gml')
    g2 = nx.read_gml('data/worldwideweb/worldwideweb.gml')

    ## build gdict
    gdict = {}
    gdict['gs'] = [g1, g2]
    gdict['names'] = ['Internet', 'World Wide Web']
    gdict['descs'] = ['Topology map of the Internet', 'Sample of web page net based on University of Notre Dame pages']

    return gdict


def generate_save_synth_graphs(n_graphs):
    '''
    Creates n_graphs of Erdos-Renyi (ER) and Preferential Attachment (PA or Barabasi-Albert).
    All graphs have N = 10000 nodes and |E| = 20000 edges for <k> = 4.
    Saves graphs to data/ER/{g_index}.gt or data/PA/{g_index}.gt
    
    
    params:
        n_graphs : int
    '''
    
    N = 10000
    m = 2  # for PA : each new node is attached to 2 other existing nodes
           #          in order to get |E| = 20000 and therefore <k> = 4
    
    for g_index in range(n_graphs):
        ## generate ER graph 
        g_er = gt.random_graph(N, lambda : 4, directed=False)
        gt.random_rewire(g_er, model='erdos')


        ## generate PA graph
        # try seeding with 
        g_pa = gt.price_network(N, m=m, directed=False)
        
        ## save graphs
        g_er.save(f'data/ER/{g_index:02d}.gt')
        g_pa.save(f'data/PA/{g_index:02d}.gt')


def load_synth_graphs():
    '''
    Loads saved graphs and returns dictionary containing graph info
    
    returns:
        gdict : {'N'        : number of nodes;
                 'n_graphs' : number of graphs for ER, PA;
                 'ER'       : list of ER graphs;
                 'PA'       : list of PA graphs}
                 
    '''
    root_dir_er = 'data/ER/'
    root_dir_pa = 'data/PA/'
    ers = glob.glob('*.gt', root_dir=root_dir_er)
    pas = glob.glob('*.gt', root_dir=root_dir_pa)
    ers.sort()
    pas.sort()
    
    er_list = []
    pa_list = []
    
    for i in range(len(ers)):
        g_er = gt.Graph()
        g_pa = gt.Graph()
        
        g_er.load(f'{root_dir_er}{ers[i]}')
        g_pa.load(f'{root_dir_pa}{pas[i]}')
        
        er_list.append(g_er)
        pa_list.append(g_pa)
        
    ## build gdict
    gdict = {}
    gdict['N'] = 10000
    gdict['n_graphs'] = len(ers)
    gdict['ER'] = er_list
    gdict['PA'] = pa_list
    
    return gdict


def load_real_graphs():
    ## topological map of the internet
    ## https://networks.skewed.de/net/route_views
    g1 = gt.collection.ns["route_views/19971108"]

    ## Sample of World-Wide Web based on University of Notre Dame pages
    ## https://networks.skewed.de/net/notre_dame_web
    g2 = gt.collection.ns["notre_dame_web"]
    
    ## build gdict
    gdict = {}
    gdict['gs'] = [g1, g2]
    gdict['names'] = ['Internet', 'World Wide Web']
    gdict['descs'] = ['Topology map of the Internet', 'Sample of web page net based on University of Notre Dame pages']
    
    return gdict


def main():
    n_graphs = 10
    nx_gen_save_synth_graphs(n_graphs)
    nx_save_real_graphs()


if __name__ == '__main__':
    main()
