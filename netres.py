import numpy as np
import matplotlib.pyplot as plt
import graph_tool.all as gt
import networkx as nx
import graphblas_algorithms as ga
from tqdm import tqdm
import sys
import argparse
from graph_saveload import load_synth_graphs, nx_load_synth_graphs
from graph_saveload import load_real_graphs, nx_load_real_graphs


def eprint(*args, **kwargs):
    # print to stderr to be consistent with how tqdm prints to stderr
    print(*args, file=sys.stderr, **kwargs)
    

def _nx2gt(nxg):
    N = nx.number_of_nodes(nxg)
    gtg = gt.Graph(directed=False)
    gtg.add_vertex(N)  ## adds N vertices
    vmap = {}
    for i, v in enumerate(nxg.nodes()):
        vmap[v] = i
        
    for vi, vj in nxg.edges():
        i = vmap[vi]
        j = vmap[vj]
        gtg.add_edge(i, j)
        
    return gtg


def _diameter(g):
    ## TODO
    ## convert to graphblas?
    lengths = dict(nx.all_pairs_shortest_path_length(g))
    diam = 0
    for i in range(len(V) - 1):
        for j in range(i, len(V)):
            vi = V[i]
            vj = V[j]
            d = lengths[vi][vj]
            if d > diam:
                diam = d
                
    return diam


def _pseudo_diameter(g):
    ## create gt graph and run gt.pseudo_diameter()
    g2 = _nx2gt(g)
    pdiam, _ = gt.pseudo_diameter(g2)
    return pdiam


## wrapper for core loop to show/not show the tqdm bar
def nx_synth_removals(gdict, fracs, nperms, statistic, seed=None, progress=False):
    if progress:
        total_iters = 2 * gdict['n_graphs'] * nperms * fracs.shape[0]
        with tqdm(total=total_iters) as pbar:
            res = _nx_synth_removals(gdict, fracs, nperms, statistic, seed, progress, pbar)
    else:
        pbar = None
        res = _nx_synth_removals(gdict, fracs, nperms, statistic, seed, progress, pbar)
    
    rand_diams, targ_diams, rand_Ss, targ_Ss, rand_s_means, targ_s_means = res
    
    return rand_diams, targ_diams, rand_Ss, targ_Ss, rand_s_means, targ_s_means


def _nx_synth_removals(gdict, fracs, nperms, statistic, seed, progress, pbar):
    if statistic == 'all':
        _stat = 0
    elif statistic == 'diam':
        _stat = 1
    elif statistic == 'compsizes':
        _stat = 2
    else:
        raise Exception('Bad statistic inputted')

    N = gdict['N']
    n_graphs = gdict['n_graphs']

    prng = np.random.default_rng(seed=seed)
    
    V = np.arange(N)  ## list of nodes to successively mask over for node removals
    
    ## create seq of nodes to permute then remove
    rand_remove_seq = np.arange(N)
    
    # stores all resilience statistic info for all combs of ER/PA, graph, removal permutation, frac removed
    # denote the ndarray shape as (t, j, p, f)
    # ER index is t=0, PA index is t=1
    # graph is the generated synth graph
    rand_diams = np.zeros((2, n_graphs, nperms, fracs.shape[0]))
    targ_diams = np.zeros((2, n_graphs, nperms, fracs.shape[0]))
    
    rand_Ss = np.zeros((2, n_graphs, nperms, fracs.shape[0]))
    targ_Ss = np.zeros((2, n_graphs, nperms, fracs.shape[0]))
    
    rand_s_means = np.zeros((2, n_graphs, nperms, fracs.shape[0]))
    targ_s_means = np.zeros((2, n_graphs, nperms, fracs.shape[0]))
    
    N_f = np.round(N * fracs).astype(int)  # num of total nodes removed based on nfrac
    N_diff = np.zeros(N_f.shape, dtype=int)  # num of nodes to remove at each step to range from f_start to f_end
    for i in range(1, N_diff.shape[0]):
        N_diff[i] = N_f[i] - N_f[i-1]
        
    for t, g_type in enumerate(['ER', 'PA']):
        for j in range(n_graphs):
            g = gdict[g_type][j]

            ## for some reason, loaded graph nodes are labeled as strings
            ## so relabel everything using integers
            g = nx.convert_node_labels_to_integers(g)

            targ_remove_seq_blocks = [[] for _ in range(N)]
            for i, k in nx.degree(g):
                targ_remove_seq_blocks[k].append(i)

            targ_remove_seq_blocks = targ_remove_seq_blocks[::-1]

            targ_remove_seq_blocks = [np.array(targ_remove_seq_block, dtype=int) \
                                      for targ_remove_seq_block in targ_remove_seq_blocks]

            for p in range(nperms):
                ## shuffle node removal order
                prng.shuffle(rand_remove_seq)
                for targ_remove_seq_block in targ_remove_seq_blocks:
                    prng.shuffle(targ_remove_seq_block)
                targ_remove_seq = np.concatenate(targ_remove_seq_blocks)

                ckpt = 0  # checkpoint into remove_seqs; used to keep track of which nodes have been removed

                for f in range(N_diff.shape[0]):
                    r = N_diff[f]  # number of nodes to remove in this iteration

                    rand_V_keep = V[rand_remove_seq[ckpt+r:]]
                    targ_V_keep = V[targ_remove_seq[ckpt+r:]]

                    rand_subg = nx.subgraph(g, rand_V_keep)
                    targ_subg = nx.subgraph(g, targ_V_keep)       
                    
                    if _stat == 0:  ## get all stats
                        ## convert to graphblas for faster diameter computations
#                         rand_ga_ccVs = [(ga.Graph.from_networkx(rand_subg.subgraph(c)), list(c)) \
#                                        for c in nx.connected_components(rand_subg)]
#                         targ_ga_ccVs = [(ga.Graph.from_networkx(targ_subg.subgraph(c)), list(c)) \
#                                        for c in nx.connected_components(targ_subg)]

#                         rand_diam = max([_diameter(ga_cc, ga_cc_V) \
#                                          for (ga_cc, ga_cc_V) in rand_ga_ccVs])
#                         targ_diam = max([_diameter(ga_cc, ga_cc_V) \
#                                          for (ga_cc, ga_cc_V) in targ_ga_ccVs])

                        rand_diam = max([_diameter(rand_subg.subgraph(c)) \
                                         for c in nx.connected_components(rand_subg)])
                        targ_diam = max([_diameter(targ_subg.subgraph(c)) \
                                         for c in nx.connected_components(targ_subg)])
                        
                        rand_S, rand_s_mean = nx_calc_component_size_stats(rand_subg, N)
                        targ_S, targ_s_mean = nx_calc_component_size_stats(targ_subg, N)
                        
                    elif _stat == 1:  ## get diams
                        ## convert to graphblas for faster diameter computations
#                         rand_ga_ccVs = [(ga.Graph.from_networkx(rand_subg.subgraph(c)), list(c)) \
#                                        for c in nx.connected_components(rand_subg)]
#                         targ_ga_ccVs = [(ga.Graph.from_networkx(targ_subg.subgraph(c)), list(c)) \
#                                        for c in nx.connected_components(targ_subg)]

#                         rand_diam = max([_diameter(ga_cc, ga_cc_V) \
#                                          for (ga_cc, ga_cc_V) in rand_ga_ccVs])
#                         targ_diam = max([_diameter(ga_cc, ga_cc_V) \
#                                          for (ga_cc, ga_cc_V) in targ_ga_ccVs])
                        
                        rand_diam = max([_diameter(rand_subg.subgraph(c)) \
                                         for c in nx.connected_components(rand_subg)])
                        targ_diam = max([_diameter(targ_subg.subgraph(c)) \
                                         for c in nx.connected_components(targ_subg)])
                        
                        rand_S, rand_s_mean = (0, 0)
                        targ_S, targ_s_mean = (0, 0)
                    
                    elif _stat == 2:  ## get comp sizes
                        rand_diam = 0
                        targ_diam = 0
                        
                        rand_S, rand_s_mean = nx_calc_component_size_stats(rand_subg, N)
                        targ_S, targ_s_mean = nx_calc_component_size_stats(targ_subg, N)

                    rand_diams[t, j, p, f] = rand_diam
                    targ_diams[t, j, p, f] = targ_diam

                    rand_Ss[t, j, p, f] = rand_S
                    targ_Ss[t, j, p, f] = targ_S

                    rand_s_means[t, j, p, f] = rand_s_mean
                    targ_s_means[t, j, p, f] = targ_s_mean

                    ckpt += r
                    
                    if progress:
                        pbar.update(1)

    ## reshape to consider each graph + permutation paring as independent samples
    rand_diams = rand_diams.reshape((2, -1, fracs.shape[0]))  
    targ_diams = targ_diams.reshape((2, -1, fracs.shape[0]))
    
    rand_Ss = rand_Ss.reshape((2, -1, fracs.shape[0]))
    targ_Ss = targ_Ss.reshape((2, -1, fracs.shape[0]))
    
    rand_s_means = rand_s_means.reshape((2, -1, fracs.shape[0]))
    targ_s_means = targ_s_means.reshape((2, -1, fracs.shape[0]))
    
    return rand_diams, targ_diams, rand_Ss, targ_Ss, rand_s_means, targ_s_means


## wrapper for core loop to show/not show the tqdm bar
def nx_real_removals(gdict, fracs, nperms, statistic, seed=None, progress=False):
    if progress:
        total_iters = 2 * nperms * fracs.shape[0]
        with tqdm(total=total_iters) as pbar:
            res = _nx_real_removals(gdict, fracs, nperms, statistic, seed, progress, pbar)
    else:
        pbar = None
        res = _nx_synth_removals(gdict, fracs, nperms, statistic, seed, progress, pbar)
    
    rand_diams, targ_diams, rand_Ss, targ_Ss, rand_s_means, targ_s_means = res
    
    return rand_diams, targ_diams, rand_Ss, targ_Ss, rand_s_means, targ_s_means


def _nx_real_removals(gdict, fracs, nperms, statistic, seed, progress, pbar):
    if statistic == 'all':
        _stat = 0
    elif statistic == 'diam':
        _stat = 1
    elif statistic == 'compsizes':
        _stat = 2
    else:
        raise Exception('Bad statistic inputted')

    prng = np.random.default_rng(seed=seed)
    
    Ns = [nx.number_of_nodes(g) for g in gdict['gs']]
    rand_remove_idxs_list = [np.arange(N) for N in Ns]
    
    # stores all resilience statistic info for all combs of internet/www, removal permutation, frac removed
    # denote the ndarray shape as (t, p, f)
    # internet index is t=0, www index is t=1
    rand_diams = np.zeros((2, nperms, fracs.shape[0]))
    targ_diams = np.zeros((2, nperms, fracs.shape[0]))
    
    rand_Ss = np.zeros((2, nperms, fracs.shape[0]))
    targ_Ss = np.zeros((2, nperms, fracs.shape[0]))
    
    rand_s_means = np.zeros((2, nperms, fracs.shape[0]))
    targ_s_means = np.zeros((2, nperms, fracs.shape[0]))
    
    Ns_f = [np.round(N * fracs).astype(int) for N in Ns]
    Ns_diff = [np.zeros(N_f.shape, dtype=int) for N_f in Ns_f]
    for i, N_diff in enumerate(Ns_diff):
        N_f = Ns_f[i]
        for j in range(1, N_diff.shape[0]):
            N_diff[j] = N_f[j] - N_f[j-1]
    
    for t, (g, name) in enumerate(zip(gdict['gs'], gdict['names'])):
        N = Ns[t]
        N_diff = Ns_diff[t]
        
        ## for some reason, loaded graph nodes are labeled as strings
        ## so relabel everything using integers
        g = nx.convert_node_labels_to_integers(g)
        
        V = np.arange(N)  ## list of nodes to successively mask over for node removals
    
        ## create seq of nodes to permute then remove
        rand_remove_seq = np.arange(N)

        targ_remove_seq_blocks = [[] for _ in range(N)]
        for i, k in nx.degree(g):
            targ_remove_seq_blocks[k].append(i)

        targ_remove_seq_blocks = targ_remove_seq_blocks[::-1]

        targ_remove_seq_blocks = [np.array(targ_remove_seq_block, dtype=int) \
                                  for targ_remove_seq_block in targ_remove_seq_blocks]
            
        for p in range(nperms):
            prng.shuffle(rand_remove_seq)
            
            for targ_remove_seq_block in targ_remove_seq_blocks:
                prng.shuffle(targ_remove_seq_block)
                
            targ_remove_seq = np.concatenate(targ_remove_seq_blocks)
            
            ckpt = 0
            
            for f in range(N_diff.shape[0]):
                r = N_diff[f]

                rand_V_keep = V[rand_remove_seq[ckpt+r:]]
                targ_V_keep = V[targ_remove_seq[ckpt+r:]]

                rand_subg = nx.subgraph(g, rand_V_keep)
                targ_subg = nx.subgraph(g, targ_V_keep)
                
                if _stat == 0:  ## get all stats
                    ## convert to graphblas for faster diameter computations
#                     rand_ga_ccVs = [(ga.Graph.from_networkx(rand_subg.subgraph(c)), list(c)) \
#                                    for c in nx.connected_components(rand_subg)]
#                     targ_ga_ccVs = [(ga.Graph.from_networkx(targ_subg.subgraph(c)), list(c)) \
#                                    for c in nx.connected_components(targ_subg)]

#                     rand_diam = max([_diameter(ga_cc, ga_cc_V) \
#                                      for (ga_cc, ga_cc_V) in rand_ga_ccVs])
#                     targ_diam = max([_diameter(ga_cc, ga_cc_V) \
#                                      for (ga_cc, ga_cc_V) in targ_ga_ccVs])
                    
                    rand_diam = max([_diameter(rand_subg.subgraph(c)) \
                                         for c in nx.connected_components(rand_subg)])
                    targ_diam = max([_diameter(targ_subg.subgraph(c)) \
                                     for c in nx.connected_components(targ_subg)])

                    rand_S, rand_s_mean = nx_calc_component_size_stats(rand_subg, N)
                    targ_S, targ_s_mean = nx_calc_component_size_stats(targ_subg, N)

                elif _stat == 1:  ## get diams
                    ## convert to graphblas for faster diameter computations
#                     rand_ga_ccVs = [(ga.Graph.from_networkx(rand_subg.subgraph(c)), list(c)) \
#                                    for c in nx.connected_components(rand_subg)]
#                     targ_ga_ccVs = [(ga.Graph.from_networkx(targ_subg.subgraph(c)), list(c)) \
#                                    for c in nx.connected_components(targ_subg)]

#                     rand_diam = max([_diameter(ga_cc, ga_cc_V) \
#                                      for (ga_cc, ga_cc_V) in rand_ga_ccVs])
#                     targ_diam = max([_diameter(ga_cc, ga_cc_V) \
#                                      for (ga_cc, ga_cc_V) in targ_ga_ccVs])
                    
                    rand_diam = max([_diameter(rand_subg.subgraph(c)) \
                                         for c in nx.connected_components(rand_subg)])
                    targ_diam = max([_diameter(targ_subg.subgraph(c)) \
                                     for c in nx.connected_components(targ_subg)])

                    rand_S, rand_s_mean = (0, 0)
                    targ_S, targ_s_mean = (0, 0)

                elif _stat == 2:  ## get comp sizes
                    rand_diam = 0
                    targ_diam = 0

                    rand_S, rand_s_mean = nx_calc_component_size_stats(rand_subg, N)
                    targ_S, targ_s_mean = nx_calc_component_size_stats(targ_subg, N)
                
                rand_diams[t, p, f] = rand_diam
                targ_diams[t, p, f] = targ_diam

                rand_Ss[t, p, f] = rand_S
                targ_Ss[t, p, f] = targ_S

                rand_s_means[t, p, f] = rand_s_mean
                targ_s_means[t, p, f] = targ_s_mean
                
                ckpt += r
                
                if progress:
                    pbar.update(1)
                
    return rand_diams, targ_diams, rand_Ss, targ_Ss, rand_s_means, targ_s_means


def nx_calc_component_size_stats(g, N_total):
    N = nx.number_of_nodes(g)

    Cs = np.array([len(C) for C in nx.connected_components(g)])
    i = np.argmax(Cs)
    Cmax = Cs[i]
    s_total = Cs.sum() - Cmax
    k = Cs.shape[0] - 1
    
    S = Cmax / N_total
    s_mean =  s_total / k if k > 0 else np.nan
    
    return S, s_mean

    
def nx_plot_real_deg_dist(gdict):
    fig, axs = plt.subplots(ncols=2, figsize=(12, 6))
    fig.suptitle('Real Network Degree Distributions (log-log scale)')
    for i, (g, name) in enumerate(zip(gdict['gs'], gdict['names'])):
        ax = axs[i]
        counts = nx.degree_histogram(g)
        kmean = np.mean(counts)
        
        ax.plot(counts, 'b.', linestyle='')
        ax.text(0.8, 0.8, f'<k> = {kmean:.02f}', fontsize=10)
        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)
        ax.set_xlabel('Degree')
        ax.set_ylabel('Counts')
        ax.set_title(name)
        
    plt.savefig('nx_real-deg-dist.png')
    plt.clf()


def nx_plot_synth_deg_dist(gdict):
    n_graphs = gdict['n_graphs']
    for t, g_type in enumerate(['ER', 'PA']):
        g_list = gdict[g_type]
        
        fig, axs = plt.subplots(nrows=4, ncols=3, figsize=(16, 15))
        fig.suptitle(f'{g_type} Degree Distributions')
        
        for k in range(n_graphs):
            g = g_list[k]
            
            i, j = divmod(k, 3)
            ax = axs[i, j]
            
            # counts, bins = gt.vertex_hist(g, 'out')
            counts = nx.degree_histogram(g)
            ax.plot(counts, 'b.', linestyle='')
            ax.set_xlabel('Degree')
            ax.set_ylabel('Counts')
            if g_type == 'PA':
                ax.set_xscale('log', base=10)
                ax.set_yscale('log', base=10)
                
        plt.savefig(f'nx-{g_type}-deg-dist.png')
        plt.clf()
    
    
def nx_plot_synth_diam_results(rand_diams, targ_diams, fracs):
    q = [0.025, 0.975]  # quantiles to compute
    
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), sharey=True)
    offset = 0.0005  # for better plotting so symbols don't overlap too much
    
    _diam_max = np.max((rand_diams.max(), targ_diams.max()))
    rand_diam_means = np.nanmean(rand_diams, axis=1)
    rand_diam_quantiles = np.nanquantile(rand_diams, q=q, axis=1)  ## returns shape (q.shape[0], 2, fracs.shape[0])
    
    targ_diam_means = np.nanmean(targ_diams, axis=1)
    targ_diam_quantiles = np.nanquantile(targ_diams, q=q, axis=1)
    
    ## plot results from single graph + idx permutation combination
    axs[0].plot(fracs, rand_diams[0, 0, :], 'b^',
                linestyle='', markerfacecolor='none',
                label='ER Random Removal')  # ER random removal
    
    axs[0].plot(fracs+offset, targ_diams[0, 0, :], 'rD',
                linestyle='', markerfacecolor='none',
                label='ER Targeted Removal')  # ER targeted removal
    
    axs[0].plot(fracs, rand_diams[1, 0, :], 'bs',
                linestyle='', markerfacecolor='none',
                label='PA Random Removal')  # PA random removal
    
    axs[0].plot(fracs+offset, targ_diams[1, 0, :], 'ro',
                linestyle='', markerfacecolor='none',
                label='PA Targeted Removal')  # PA targeted removal
    
    axs[0].set_title('Single Graph + Permutation Combination')
    axs[0].set_xlabel('Fraction of total nodes removed')
    axs[0].set_ylabel('Diameter')
    axs[0].set_yticks(ticks=np.arange(_diam_max))
    axs[0].legend()
    
    
    ## plot mean + quantile results using all graph + idx permutation combinations
    axs[1].plot(fracs, rand_diam_means[0], 'b^', linestyle='', markerfacecolor='none')  # ER random removal means
    axs[1].fill_between(fracs, rand_diam_quantiles[0, 0], rand_diam_quantiles[1, 0],
                        color='blue', alpha=0.15)  # ER random removal quantiles
    
    axs[1].plot(fracs+offset, targ_diam_means[0], 'rD', linestyle='', markerfacecolor='none')  # ER random removal means
    axs[1].fill_between(fracs+offset, targ_diam_quantiles[0, 0], targ_diam_quantiles[1, 0],
                        color='red', alpha=0.15)  # ER targeted removal quantiles
    
    axs[1].plot(fracs, rand_diam_means[1], 'bs', linestyle='', markerfacecolor='none')  # PA random removal means
    axs[1].fill_between(fracs, rand_diam_quantiles[0, 1], rand_diam_quantiles[1, 1],
                        color='blue', alpha=0.15)  # ER random removal quantiles
    
    axs[1].plot(fracs+offset, targ_diam_means[1], 'ro', linestyle='', markerfacecolor='none')  # ER random removal means
    axs[1].fill_between(fracs+offset, targ_diam_quantiles[0, 1], targ_diam_quantiles[1, 1],
                        color='red', alpha=0.15)  # ER targeted removal quantiles
    
    axs[1].set_title(f'Mean and [{q[0]}, {q[1]}] Quantiles of\nMultiple Graph + Permutation Combinations')
    axs[1].set_xlabel('Fraction of total nodes removed')
    axs[1].set_ylabel('Diameter')
    axs[1].set_yticks(ticks=np.arange(_diam_max))
    
    fig.savefig('nx_synth_diams.png')
    plt.clf()
    
    
def nx_plot_synth_frag_results(rand_Ss, targ_Ss, rand_s_means, targ_s_means, fracs):
    q = [0.025, 0.975]
    
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))
    offset = 0.0005  # for better plotting so symbols don't overlap too much
    
    rand_Ss_means = np.nanmean(rand_Ss, axis=1)
    rand_Ss_quantiles = np.nanquantile(rand_Ss, q=q, axis=1)  ## returns shape (q.shape[0], 2, fracs.shape[0])
    
    targ_Ss_means = np.nanmean(targ_Ss, axis=1)
    targ_Ss_quantiles = np.nanquantile(targ_Ss, q=q, axis=1)
    
    rand_s_means_means = np.nanmean(rand_s_means, axis=1)
    rand_s_means_quantiles = np.nanquantile(rand_s_means, q=q, axis=1)  ## returns shape (q.shape[0], 2, fracs.shape[0])
    
    targ_s_means_means = np.nanmean(targ_s_means, axis=1)
    targ_s_means_quantiles = np.nanquantile(targ_s_means, q=q, axis=1)
    
    
    gtypes = ['ER', 'PA']
    datanames = ['Single Graph + Permutation Combination',
                 f'Mean and [{q[0]}, {q[1]}] Quantiles of\nMultiple Graph + Permutation Combinations']

    for i, gtype in enumerate(gtypes):
        ## plot single graph + permutation combination
        j = 0
        ax = axs[i, j]
        

        ax.plot(fracs, rand_Ss[i, 0, :], 'bs',
                linestyle='', markerfacecolor='none',
                label='S Random Removal')

        ax.plot(fracs, rand_s_means[i, 0, :], 'bs',
                linestyle='',
                label='<s> Random Removal')

        ax.plot(fracs+offset, targ_Ss[i, 0, :], 'ro',
                linestyle='', markerfacecolor='none',
                label='S Targeted Removal')

        ax.plot(fracs+offset, targ_s_means[i, 0, :], 'ro',
                linestyle='',
                label='<s> Targeted Removal')
        
        ax.set_title(f'{gtype} {datanames[j]}')
        ax.set_xlabel('Fraction of total nodes removed')
        ax.set_ylabel('S and <s>')
        ax.legend()

        ## plot mean + quantile of multiple graph + permutation combinations
        j = 1
        ax = axs[i, j]
        
        
        ax.plot(fracs, rand_Ss_means[i], 'bs',
                linestyle='', markerfacecolor='none')
        ax.fill_between(fracs, rand_Ss_quantiles[0, i], rand_Ss_quantiles[1, i],
                        color='blue', alpha=0.15)
        
        ax.plot(fracs, rand_s_means_means[i], 'bs',
                linestyle='')
        ax.fill_between(fracs, rand_s_means_quantiles[0, i], rand_s_means_quantiles[1, i],
                        color='blue', alpha=0.15)
        
        ax.plot(fracs+offset, targ_Ss_means[i], 'ro',
                linestyle='', markerfacecolor='none')
        ax.fill_between(fracs+offset, targ_Ss_quantiles[0, i], targ_Ss_quantiles[1, i],
                        color='red', alpha=0.15)
        
        ax.plot(fracs+offset, targ_s_means_means[i], 'ro',
                linestyle='')
        ax.fill_between(fracs+offset, targ_s_means_quantiles[0, i], targ_s_means_quantiles[1, i],
                        color='red', alpha=0.15)
        
        ax.set_title(f'{gtype} {datanames[j]}')
        ax.set_xlabel('Fraction of total nodes removed')
        ax.set_ylabel('S and <s>')
        
    fig.savefig('nx_synth_frags.png')
    plt.clf()

    
def nx_plot_real_diam_results(rand_diams, targ_diams, fracs):
    q = [0.025, 0.975]  # quantiles to compute
    
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 10), sharey=True)
    offset = 0.0005  # for better plotting so symbols don't overlap too much
    
    names = ['Internet', 'World Wide Web']
    psingle = 'Single Permutation'
    pall = f'Mean + [{q[0]}, {q[1]}] Quantiles of\nMultiple Permutations'
    plottypes = [psingle, pall]
    
    _diam_maxs = [np.max((rand_diams[i, :, :].max(), targ_diams[i, :, :].max())) for i in range(2)]
    rand_diam_means = np.nanmean(rand_diams, axis=1)
    rand_diam_quantiles = np.nanquantile(rand_diams, q=q, axis=1)  ## returns shape (q.shape[0], 2, fracs.shape[0])
    
    targ_diam_means = np.nanmean(targ_diams, axis=1)
    targ_diam_quantiles = np.nanquantile(targ_diams, q=q, axis=1)
    
    for i, name in enumerate(names):
        ax_psingle = axs[i, 0]
        ax_pall = axs[i, 1]
        
        
        ax_psingle.plot(fracs, rand_diams[i, 0, :], 'bs',
                        linestyle='', markerfacecolor='none',
                        label=f'{name} Random Removal')
            
        ax_psingle.plot(fracs+offset, targ_diams[i, 0, :], 'ro',
                        linestyle='', markerfacecolor='none',
                        label=f'{name} Targeted Removal')
        
        ax_pall.plot(fracs, rand_diam_means[i], 'bs',
                     linestyle='', markerfacecolor='none')
        ax_pall.fill_between(fracs, rand_diam_quantiles[0, i], rand_diam_quantiles[1, i],
                             color='blue', alpha=0.15)
        
        ax_pall.plot(fracs+offset, targ_diam_means[i], 'ro',
                     linestyle='', markerfacecolor='none')
        ax_pall.fill_between(fracs, targ_diam_quantiles[0, i], targ_diam_quantiles[1, i],
                             color='red', alpha=0.15)
        
        _diam_max_psingle = np.max((rand_diams[i, 0, :].max(), targ_diams[i, 0, :].max()))
        _diam_max_pall = np.round(np.max((rand_diam_means[i].max(), targ_diam_means[i].max())))

        ax_psingle.set_title(f'{name} {psingle}')
        ax_psingle.set_xlabel('Fraction of total nodes removed')
        ax_psingle.set_ylabel('Diameter')
        ax_psingle.set_ylim((-5, _diam_max_psingle+5))
        ax_psingle.legend()
        
        ax_pall.set_title(f'{name} {pall}')
        ax_pall.set_xlabel('Fraction of total nodes removed')
        ax_pall.set_ylabel('Diameter')
        ax_pall.set_ylim((-5, _diam_max_pall+5))
        
    fig.savefig('nx_real_diams.png')
    plt.clf()


def nx_plot_real_frag_results(rand_Ss, targ_Ss, rand_s_means, targ_s_means, fracs):
    q = [0.025, 0.975]
    
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))
    offset = 0.0005  # for better plotting so symbols don't overlap too much
    
    names = ['Internet', 'World Wide Web']
    psingle = 'Single Permutation'
    pall = f'Mean + [{q[0]}, {q[1]}] Quantiles of\nMultiple Permutations'
    plottypes = [psingle, pall]
    
    rand_Ss_means = np.nanmean(rand_Ss, axis=1)
    rand_Ss_quantiles = np.nanquantile(rand_Ss, q=q, axis=1)  ## returns shape (q.shape[0], 2, fracs.shape[0])
    
    targ_Ss_means = np.nanmean(targ_Ss, axis=1)
    targ_Ss_quantiles = np.nanquantile(targ_Ss, q=q, axis=1)
    
    rand_s_means_means = np.nanmean(rand_s_means, axis=1)
    rand_s_means_quantiles = np.nanquantile(rand_s_means, q=q, axis=1)  ## returns shape (q.shape[0], 2, fracs.shape[0])
    
    targ_s_means_means = np.nanmean(targ_s_means, axis=1)
    targ_s_means_quantiles = np.nanquantile(targ_s_means, q=q, axis=1)
    
    for i, name in enumerate(names):
        ## plot single graph + permutation combination
        j = 0
        ax = axs[i, j]
        

        ax.plot(fracs, rand_Ss[i, 0, :], 'bs',
                linestyle='', markerfacecolor='none',
                label='S Random Removal')

        ax.plot(fracs, rand_s_means[i, 0, :], 'bs',
                linestyle='',
                label='<s> Random Removal')

        ax.plot(fracs+offset, targ_Ss[i, 0, :], 'ro',
                linestyle='', markerfacecolor='none',
                label='S Targeted Removal')

        ax.plot(fracs+offset, targ_s_means[i, 0, :], 'ro',
                linestyle='',
                label='<s> Targeted Removal')
        
        ax.set_title(f'{name} {plottypes[j]}')
        ax.set_xlabel('Fraction of total nodes removed')
        ax.set_ylabel('S and <s>')
        ax.legend()

        ## plot mean + quantile of multiple graph + permutation combinations
        j = 1
        ax = axs[i, j]
        
        
        ax.plot(fracs, rand_Ss_means[i], 'bs',
                linestyle='', markerfacecolor='none')
        ax.fill_between(fracs, rand_Ss_quantiles[0, i], rand_Ss_quantiles[1, i],
                        color='blue', alpha=0.15)
        
        ax.plot(fracs, rand_s_means_means[i], 'bs',
                linestyle='')
        ax.fill_between(fracs, rand_s_means_quantiles[0, i], rand_s_means_quantiles[1, i],
                        color='blue', alpha=0.15)
        
        ax.plot(fracs+offset, targ_Ss_means[i], 'ro',
                linestyle='', markerfacecolor='none')
        ax.fill_between(fracs+offset, targ_Ss_quantiles[0, i], targ_Ss_quantiles[1, i],
                        color='red', alpha=0.15)
        
        ax.plot(fracs+offset, targ_s_means_means[i], 'ro',
                linestyle='')
        ax.fill_between(fracs+offset, targ_s_means_quantiles[0, i], targ_s_means_quantiles[1, i],
                        color='red', alpha=0.15)
        
        ax.set_title(f'{name} {plottypes[j]}')
        ax.set_xlabel('Fraction of total nodes removed')
        ax.set_ylabel('S and <s>')
        
    fig.savefig('nx_real_frags.png')
    plt.clf()


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--progress', '-p', action='store_true')
    parser.add_argument('--nperms', type=int, default=1)
    parser.add_argument('--seed', '-s', type=int, default=1000)
    args = parser.parse_args()
    
    # progress = True

    ## set some hyperparams
    f_synth_diam_start = 0
    f_synth_diam_end = 0.05
    num_synth_diam = 21
    # num_synth_diam = 6
    ## fraction of nodes to remove
    fracs_synth_diam = np.linspace(f_synth_diam_start, f_synth_diam_end, num=num_synth_diam)  # includes endpoints

    f_real_diam_start = 0
    f_real_diam_end = 0.025
    num_real_diam = 21
    # num_real_diam = 3
    fracs_real_diam = np.linspace(f_real_diam_start, f_real_diam_end, num=num_real_diam)

    f_synth_frag_start = 0
    f_synth_frag_end = 0.5
    num_synth_frag = 51
    # num_synth_frag = 6
    fracs_synth_frag = np.linspace(f_synth_frag_start, f_synth_frag_end, num=num_synth_frag)
    
    f_real_frag_start = 0
    f_real_frag_end = 0.15
    num_real_frag = 21
    # num_real_frag = 3
    fracs_real_frag = np.linspace(f_real_frag_start, f_real_frag_end, num=num_real_frag)
    
    # nperms = 10
    # nperms = 1
    # seed = 1000
    
    nperms = args.nperms
    progress = args.progress
    seed = args.seed
    
    print('nperms:', nperms)
    print('progress:', progress)
    print('seed:', seed)
    print('num_synth_diam:', num_synth_diam)
    print('num_real_diam:', num_real_diam)
    print('num_synth_frag:', num_synth_frag)
    print('num_real_frag:', num_real_frag)
    
    if progress:
        eprint('Loading Synth Graphs...')
    gdict_synth = nx_load_synth_graphs()
    
    if progress:
        eprint('Plotting Degree Distribution of Synth Graphs...')
    nx_plot_synth_deg_dist(gdict_synth)
    
    if progress:
        eprint('Computing Diam Stats on Synth Graphs')
    synth_res = nx_synth_removals(gdict_synth, fracs_synth_diam, nperms, 'diam', seed=seed, progress=progress)
    rand_diams, targ_diams, _, _, _, _ = synth_res
    
    if progress:
        eprint('Computing Component Size Stats on Synth Graphs')
    synth_res = nx_synth_removals(gdict_synth, fracs_synth_frag, nperms, 'compsizes', seed=seed, progress=progress)
    _, _, rand_Ss, targ_Ss, rand_s_means, targ_s_means = synth_res
    
    if progress:
        eprint('Plotting Results on Synth Graphs...')
    nx_plot_synth_diam_results(rand_diams, targ_diams, fracs_synth_diam)
    nx_plot_synth_frag_results(rand_Ss, targ_Ss, rand_s_means, targ_s_means, fracs_synth_frag)

    if progress:
        eprint('Loading Real Graphs...')
    gdict_real = nx_load_real_graphs()
    
    if progress:
        eprint('Plotting Degree Distribution of Real Graphs...')
    nx_plot_real_deg_dist(gdict_real)

    if progress:
        eprint('Computing Diam Stats on Real Graphs')
    real_res = nx_real_removals(gdict_real, fracs_real_diam, nperms, 'diam', seed=seed, progress=progress)
    rand_diams, targ_diams, _, _, _, _ = real_res
    
    if progress:
        eprint('Computing Component Size Stats on Real Graphs')
    real_res = nx_real_removals(gdict_real, fracs_real_frag, nperms, 'compsizes', seed=seed, progress=progress)
    _, _, rand_Ss, targ_Ss, rand_s_means, targ_s_means = real_res
    
    if progress:
        eprint('Plotting Results on Real Graphs...')
    nx_plot_real_diam_results(rand_diams, targ_diams, fracs_real_diam)
    nx_plot_real_frag_results(rand_Ss, targ_Ss, rand_s_means, targ_s_means, fracs_real_frag)


if __name__ == '__main__':
    main()

