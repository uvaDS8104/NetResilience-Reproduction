import numpy as np
import matplotlib.pyplot as plt
import graph_tool.all as gt
import networkx as nx
import graphblas_algorithms as ga
from tqdm import tqdm
import sys
import argparse
import os
import pickle
import functools
import operator
from graph_saveload import nx_load_synth_graphs
from graph_saveload import nx_load_synth_multilayer_graphs
from graph_saveload import nx_load_real_graphs


def eprint(*args, **kwargs):
    # print to stderr to be consistent with how tqdm prints to stderr
    print(*args, file=sys.stderr, **kwargs)


def flatten(a):
    # flatten a nd list into a n-1 d list
    return functools.reduce(operator.iconcat, a, [])


def _calc_diameter(subg_ccs):
    # compute diameter in graphblas by taking max_{d_c} where d_c is diameter of connected component subgraph
    subgs = [ga.Graph.from_networkx(subg) for subg in subg_ccs]  ## use graphblas
    Vs = [list(subg.nodes()) for subg in subg_ccs]
    lengths = [dict(nx.all_pairs_shortest_path_length(subg)) for subg in subgs]
    
    subg_diams = np.zeros(len(subgs), dtype=int)
    for k, subg in enumerate(subgs):
        V = Vs[k]
        lendict = lengths[k]
        sub_d = 0
        for i in range(len(V) - 1):
            for j in range(i+1, len(V)):
                vi = V[i]
                vj = V[j]
                d = lendict[vi][vj]
                if d > sub_d:
                    sub_d = d
        subg_diams[k] = sub_d
        
    return subg_diams.max()


def _calc_fragmentation(subg_ccs, N_total):
    # compute S and <s>
    Cs = np.array([subg.number_of_nodes() for subg in subg_ccs])
    i = np.argmax(Cs)
    Cmax = Cs[i]
    
    Crest = np.delete(Cs, i)
    
    S = Cmax / N_total
    s_mean = Crest.mean()
    
    return np.array([S, s_mean])


def get_new_cc_subgs(cc_subgs, V):
    ## cc_subgs is a list of subgraphs, all connected components
    ## V list of nodes to create new induced subgraphs
    
    # get list of induced subgraphs which may not be ccs anymore
    new_subgs = [nx.subgraph(cc_subg, V) for cc_subg in cc_subgs]
    # get ccs for subgraphs :: list of generators of ccs
    new_subg_cc_gens = [nx.connected_components(subg) for subg in new_subgs]
    # list of lists of new induced subgraphs, which are all connected components
    a = [[nx.subgraph(subg, cc) for cc in subg_cc_gen] \
         for subg, subg_cc_gen in zip(new_subgs, new_subg_cc_gens)]
    
    new_cc_subgs = flatten(a)

    return new_cc_subgs


def get_rem_strat_shuffler(rem_strat):
    # returns the shuffling function
    if rem_strat == 'rand':
        shuffle_fn = rand_shuffler
    elif rem_strat == 'targ':
        shuffle_fn = targ_shuffler
        
    return shuffle_fn


def rand_shuffler(curr_rem_seq, cc_subgs, Nrem, prng):
    # shuffles the last Nrem elements of curr_rem_seq in-place
    prng.shuffle(curr_rem_seq[-Nrem:])


def targ_shuffler(curr_rem_seq, cc_subgs, Nrem, prng):
    # shuffles the last Nrem elements of curr_rem_seq in-place
    # bucket nodes by degree
    remove_seq_blocks = [[] for _ in range(Nrem)]
    for cc in cc_subgs:
        for i, k in nx.degree(cc):
            remove_seq_blocks[k].append(i)

    # reverse order of buckets
    remove_seq_blocks = remove_seq_blocks[::-1]

    # convert to np.ndarray
    remove_seq_blocks = [np.array(remove_seq_block, dtype=int) \
                               for remove_seq_block in remove_seq_blocks]

    # shuffle within buckets
    for remove_seq_block in remove_seq_blocks:
        prng.shuffle(remove_seq_block)
        
    _rem_shuffled = np.concatenate(remove_seq_blocks)
    
    # update curr_rem_seq
    curr_rem_seq[-Nrem:] = _rem_shuffled
    
    
def calc_stats(gs, rem_strat, fracs, nperms, stat, progress=True, snapshots=False, seed=1000):
    ## Wrapper function which sets up the output array, shuffle_fn, and tqdm progressbar
    prng = np.random.default_rng(seed=seed)
    
    # create output array
    if stat == 'diam':
        outshape = (len(gs), nperms, fracs.shape[0])
    elif stat == 'frag':
        outshape = (len(gs), nperms, fracs.shape[0], 2)  # 2 is for S and s_mean

    all_stats = np.zeros(outshape)
    all_snapshots = [[] for _ in range(len(gs))]
    
    # get node removal order shuffle function
    shuffle_fn = get_rem_strat_shuffler(rem_strat)
    
    # handle pbar
    if progress:
        total_iters = len(gs) * nperms * fracs.shape[0]
        with tqdm(total=total_iters) as pbar:
            for i, g in enumerate(gs):
                for j in range(nperms):
                    res = _calc_stats(g, shuffle_fn, fracs, stat, progress, pbar, snapshots, prng)
                    if snapshots:
                        stats_arr, snapshots_list = res
                    else:
                        stats_arr = res
                    all_stats[i, j] = stats_arr
                    
    else:
        pbar = None
        for i, g in enumerate(gs):
            for j in range(nperms):
                res = _calc_stats(g, shuffle_fn, fracs, stat, progress, pbar, snapshots, prng)
                if snapshots:
                    stats_arr, snapshots_list = res
                else:
                    stats_arr = res
                all_stats[i, j] = stats_arr
                            
    ## flatten gs and nperms dims
    if stat == 'diam':
        all_stats = all_stats.reshape((-1, fracs.shape[0]))
    elif stat == 'frag':
        all_stats = all_stats.reshape((-1, fracs.shape[0], 2))
    
    if snapshots:
        return all_stats, all_snapshots
    else:
        return all_stats

    
def _calc_stats(g, shuffle_fn, fracs, stat, progress, pbar, snapshots, prng):
    ## Calculate statistics for single graph + permutation
    N = g.number_of_nodes()
    Nrem = N  ## remaining number of nodes after removals
    V = np.arange(N)
    
    N_f = np.round(N * fracs).astype(int)  # num of total nodes removed based on fracs
    N_diff = np.zeros(N_f.shape, dtype=int)  # num of nodes to remove at each step
    
    for i in range(1, N_diff.shape[0]):
        N_diff[i] = N_f[i] - N_f[i-1]
        
    if stat == 'diam':
        stat_fn = _calc_diameter
        outshape = (fracs.shape[0],)
    elif stat == 'frag':
        stat_fn = lambda x: _calc_fragmentation(x, N)  # hack so stat_fn has same input signature
        outshape = (fracs.shape[0], 2)

    stats_arr = np.zeros(outshape)
    snapshots_list = []
    
    ## keep list of connected components rather than full graph
    cc_subgs = [nx.subgraph(g, c) for c in nx.connected_components(g)]
    
    remove_seq = np.arange(N)
    ckpt = 0  # checkpoint into remove_seqs; used to keep track of which nodes have been removed thus far

    for f in range(N_diff.shape[0]):
        r = N_diff[f]  # number of nodes to remove in this iteration
        V_keep = V[remove_seq[ckpt+r:]]
        cc_subgs = get_new_cc_subgs(cc_subgs, V_keep)
        
        ## NEED TO RECOMPUTE removal order for targeted removals
        ## because induced subgraph may have different node ordering
        Nrem -= r
        shuffle_fn(remove_seq, cc_subgs, Nrem, prng)

        ## store snapshots
        if snapshots:
            g_snap = nx.Graph()
            for subg in cc_subgs:
                g_snap.update(subg)
            g_snap_gml_string = '\n'.join(nx.generate_gml(g_snap))  ## CURRENTLY NOT WORKING
            snapshots_list.append(g_snap_gml_string)

        _stat = stat_fn(cc_subgs)
        stats_arr[f] = _stat

        ckpt += r

        if progress:
            pbar.update(1)

    if snapshots:
        return stats_arr, snapshots_list
    else:
        return stats_arr


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--graphtype', '-g', choices=['er', 'pa', 'er2', 'pa2', 'erpa', 'int', 'www'])
    parser.add_argument('--ngraphs', '-n', type=int, choices=range(1, 11))
    parser.add_argument('--stat', '-s', choices=['diam', 'frag'])
    parser.add_argument('--remstrat', '-r', choices=['rand', 'targ'])
    parser.add_argument('--progress', '-p', action='store_true')
    parser.add_argument('--snapshots', action='store_true')
    parser.add_argument('--nperms', type=int, default=1)
    parser.add_argument('--seed', type=int, default=1000)
    args = parser.parse_args()
    
    graphtype = args.graphtype
    n_graphs = args.ngraphs
    stat = args.stat
    rem_strat = args.remstrat
    nperms = args.nperms
    snapshots = args.snapshots
    progress = args.progress
    seed = args.seed
    
    if progress:
        print('Supplied Arguments:')
        print('--graphtype:', graphtype)
        print('--ngraphs:', n_graphs)
        print('--stat:', stat)
        print('--remstrat:', rem_strat)
        print('--nperms:', nperms)
        print('--snapshots:', snapshots)
        print('--progress:', progress)
        print('--seed:', seed)
    
    outdir = 'outfiles/results'
    os.makedirs(outdir, exist_ok=True)
    
    if progress:
        print('Loading graphs...')
    if graphtype == 'er' or graphtype == 'pa':
        gdict = nx_load_synth_graphs()
        gs = gdict[graphtype.upper()]
        if stat == 'diam':
            f_start = 0
            f_end = 0.05
            f_num = 21
        elif stat == 'frag':
            f_start = 0
            f_end = 0.5
            f_num = 26
            
    elif graphtype == 'er2' or graphtype == 'pa2' or graphtype == 'erpa':
        gdict = nx_load_synth_multilayer_graphs()
        gs = gdict[graphtype.upper()]
        if stat == 'diam':
            f_start = 0
            f_end = 0.05
            f_num = 21
        elif stat == 'frag':
            f_start = 0
            f_end = 0.5
            f_num = 26
        
    elif graphtype == 'int' or graphtype == 'www':
        gdict = nx_load_real_graphs()
        gs, gname, gdesc = gdict[graphtype.upper()]
        if stat == 'diam':
            f_start = 0
            f_end = 0.025
            f_num = 21
        elif stat == 'frag':
            f_start = 0
            f_end = 0.15
            f_num = 26

    if progress:
        print('Using linspace arguments:')
        print('f_start:', f_start)
        print('f_end:', f_end)
        print('f_num:', f_num)
            
    fracs = np.linspace(f_start, f_end, f_num)
    
    if progress:
        print(f'Calculating {stat} on {n_graphs} {graphtype} graph(s) using {nperms} permutations of {rem_strat} removal strategy...')
    res = calc_stats(gs[:n_graphs], rem_strat, fracs, nperms, stat, 
                     progress=progress, snapshots=snapshots, seed=seed)
    if snapshots:
        all_stats, all_snapshots = res
        snapshot_outname = f'{outdir}/{graphtype}_{rem_strat}_{stat}_snapshots.pkl'
        if progress:
            print(f"Saving snapshots to '{snapshot_outname}'...")
            pickle.dump(all_snapshots, snapshot_outname)
    else:
        all_stats = res
    
    outfilename = f'{outdir}/{graphtype}_{rem_strat}_{stat}.npy'
    if progress:
        print(f"Saving statistics to '{outfilename}'...")
    np.save(outfilename, all_stats)


if __name__ == '__main__':
    main()

