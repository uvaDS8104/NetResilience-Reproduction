import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def plot_all_deg_dist(gs, graphtypes, outdir):
    for i, graphtype in enumerate(graphtypes):
        g = gs[i]
        if graphtype == 'er' or graphtype == 'er2':
            loglog = False
        else:
            loglog = True
        nx_plot_deg_dist(g, graphtype, outdir, loglog=loglog)


def nx_plot_deg_dist(g, graphtype, outdir, loglog=False):
    fig, ax = plt.subplots(figsize=(8, 8))
    fig.suptitle(f'{graphtype.upper()} Degree Distribution')
    
    counts = nx.degree_histogram(g)
    ax.plot(counts, 'b.', linestyle='')
    ax.set_xlabel('Degree')
    ax.set_ylabel('Counts')
    
    if loglog:
        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)

    plt.savefig(f'{outdir}/plots/nx_{graphtype}_degdist.png')
    plt.clf()


def plot_snapshots(snaps, fracs, ncols, title, graphtype, rem_strat, outdir,
                   step=1, sidesize=6, fontsize=32, node_size=20, width=0.3):
    _fracs = fracs[::step]
    _snaps = snaps[::step]
    
    nrows, r = divmod(_fracs.shape[0], ncols)
    if r > 0:
        nrows += 1
    
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*sidesize, nrows*sidesize))
    fig.suptitle(title, fontsize=fontsize)
    
    for i, snap in enumerate(_snaps):
        row, col = divmod(i, ncols)
        ax = axs[row, col] if nrows > 1 else axs[col]
        ax.set_title(f'f = {_fracs[i]:.2f}')
        nx.draw(snap, pos=nx.spring_layout(snap), ax=ax, node_size=node_size, width=width)
        
    fig.savefig(f'{outdir}/plots/nx_{graphtype}_{rem_strat}_snapshots.png')
    plt.clf()

    
def plot_synth_diam_results(er_rand_diams, er_targ_diams, pa_rand_diams, pa_targ_diams, fracs, outdir):
    q = [0.025, 0.975]  # quantiles to compute
    
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), sharey=True)
    offset = 0.0005  # for better plotting so symbols don't overlap too much
    
    _diam_max = np.max((er_rand_diams.max(), er_targ_diams.max(), pa_rand_diams.max(), pa_targ_diams.max()))
    
    diams_list = [er_rand_diams, er_targ_diams, pa_rand_diams, pa_targ_diams]
    labels = ['ER Random Removal', 'ER Targeted Removal', 'PA Random Removal', 'PA Targeted Removal']
    markers = ['b^', 'rD', 'bs', 'ro']
    
    for i, (diams, label, marker) in enumerate(zip(diams_list, labels, markers)):
        if i % 2 == 0:  # rands
            _fracs = fracs
            color = 'blue'
        else:
            _fracs = fracs + offset
            color = 'red'
        
        diam_means = np.nanmean(diams, axis=0)
        diam_quantiles = np.nanquantile(diams, q=q, axis=0)
        
        ## plot single graph + permutation combination
        axs[0].plot(_fracs, diams[0, :], marker,
                    linestyle='', markerfacecolor='none',
                    label=label)
        axs[0].set_title('Single Graph + Permutation Combination')
        axs[0].set_xlabel('Fraction of total nodes removed')
        axs[0].set_ylabel('Diameter')
        axs[0].set_yticks(ticks=np.arange(_diam_max))
        axs[0].legend()
        
        ## plot mean + quantile results using all graph + idx permutation combinations
        axs[1].plot(_fracs, diam_means, marker, linestyle='', markerfacecolor='none')
        axs[1].fill_between(fracs, diam_quantiles[0], diam_quantiles[1],
                            color=color, alpha=0.15)
        axs[1].set_title(f'Mean and [{q[0]}, {q[1]}] Quantiles of\nMultiple Graph + Permutation Combinations')
        axs[1].set_xlabel('Fraction of total nodes removed')
        axs[1].set_ylabel('Diameter')
        axs[1].set_yticks(ticks=np.arange(_diam_max))
        
    fig.savefig(f'{outdir}/plots/nx_synth_diams.png')
    plt.clf()
    

def plot_synth_frag_results(er_rand_frags, er_targ_frags, pa_rand_frags, pa_targ_frags, fracs, outdir):
    q = [0.025, 0.975]
    
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))
    offset = 0.0005  # for better plotting so symbols don't overlap too much
    
    gtypes = ['ER', 'PA']
    datanames = ['Single Graph + Permutation Combination',
                 f'Mean and [{q[0]}, {q[1]}] Quantiles of\nMultiple Graph + Permutation Combinations']
    rem_strats = ['Random Removal', 'Targeted Removal']
    stat_names = ['S', '<s>']
    markers = ['bs', 'ro']
    colors = ['blue', 'red']
    
    frag_stats = [[er_rand_frags, er_targ_frags],
                  [pa_rand_frags, pa_targ_frags]]
    
    for i, gtype in enumerate(gtypes):
        _fracs = fracs if i % 2 == 0 else fracs + offset
        for k, (rem_strat, color, marker) in enumerate(zip(rem_strats, colors, markers)):
            frags = frag_stats[i][k]

            for l, stat_name in enumerate(stat_names):
                label = f'{stat_name} {rem_strat}'
                mfc = 'none' if l == 0 else marker[0]
                
                frag_stat_means = np.nanmean(frags[:, :, l], axis=0)
                frag_stat_quantiles = np.nanquantile(frags[:, :, l], q=q, axis=0)
                
                ## plot single graph + permutation combination
                ax = axs[i, 0]
                ax.plot(_fracs, frags[0, :, l], marker, linestyle='',
                        markerfacecolor=mfc, label=label)
                    
                ## plot mean + quantile results using all graph + idx permutation combinations
                ax = axs[i, 1]
                ax.plot(_fracs, frag_stat_means, marker,
                        linestyle='', markerfacecolor=mfc)
                ax.fill_between(_fracs, frag_stat_quantiles[0], frag_stat_quantiles[1],
                                color=color, alpha=0.15)
                
        axs[i, 0].set_title(f'{gtype} {datanames[0]}')
        axs[i, 0].set_xlabel('Fraction of total nodes removed')
        axs[i, 0].set_ylabel('S and <s>')
        axs[i, 0].legend()
        
        axs[i, 1].set_title(f'{gtype} {datanames[1]}')
        axs[i, 1].set_xlabel('Fraction of total nodes removed')
        axs[i, 1].set_ylabel('S and <s>')
        
    fig.savefig(f'{outdir}/plots/nx_synth_frags.png')
    plt.clf()

        
def plot_real_diam_results(rand_diams, targ_diams, fracs, gtype, outdir):
    q = [0.025, 0.975]  # quantiles to compute
    
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), sharey=True)
    fig.suptitle(f'{gtype} Diameter')
    
    offset = 0.0005  # for better plotting so symbols don't overlap too much
    
    _diam_max = np.max((rand_diams.max(), targ_diams.max()))
    
    rand_diam_means = np.nanmean(rand_diams, axis=0)
    rand_diam_quantiles = np.nanquantile(rand_diams, q=q, axis=0)  ## returns shape (2, fracs.shape[0])
    
    targ_diam_means = np.nanmean(targ_diams, axis=0)
    targ_diam_quantiles = np.nanquantile(targ_diams, q=q, axis=0)
    
    ## plot results from single graph + idx permutation combination
    axs[0].plot(fracs, rand_diams[0, :], 'bs',
                linestyle='', markerfacecolor='none',
                label='Random Removal')  # PA random removal
    
    axs[0].plot(fracs+offset, targ_diams[0, :], 'ro',
                linestyle='', markerfacecolor='none',
                label='Targeted Removal')  # PA targeted removal
    
    axs[0].set_title('Single Graph + Permutation Combination')
    axs[0].set_xlabel('Fraction of total nodes removed')
    axs[0].set_ylabel('Diameter')
    axs[0].set_yticks(ticks=np.arange(_diam_max))
    axs[0].legend()
    
    ## plot mean + quantile results using all graph + idx permutation combinations
    axs[1].plot(fracs, rand_diam_means, 'bs', linestyle='', markerfacecolor='none')  # PA random removal means
    axs[1].fill_between(fracs, rand_diam_quantiles[0], rand_diam_quantiles[1],
                        color='blue', alpha=0.15)  # ER random removal quantiles
    
    axs[1].plot(fracs+offset, targ_diam_means, 'ro', linestyle='', markerfacecolor='none')  # ER random removal means
    axs[1].fill_between(fracs+offset, targ_diam_quantiles[0], targ_diam_quantiles[1],
                        color='red', alpha=0.15)  # ER targeted removal quantiles
    
    axs[1].set_title(f'Mean and [{q[0]}, {q[1]}] Quantiles of\nMultiple Graph + Permutation Combinations')
    axs[1].set_xlabel('Fraction of total nodes removed')
    axs[1].set_ylabel('Diameter')
    axs[1].set_yticks(ticks=np.arange(_diam_max))
    
    fig.savefig(f'{outdir}/plots/nx_{gtype.lower()}_diams.png')
    plt.clf()
    

def plot_real_frag_results(rand_frags, targ_frags, fracs, gtype, outdir):
    q = [0.025, 0.975]
    
    rand_Ss = rand_frags[0, :, 0]
    rand_s_means = rand_frags[0, :, 1]
    
    targ_Ss = targ_frags[0, :, 0]
    targ_s_means = targ_frags[0, :, 1]
    
    rand_Ss_means = np.nanmean(rand_frags[:, :, 0], axis=0)
    rand_Ss_quantiles = np.nanquantile(rand_frags[:, :, 0], q=q, axis=0)
    
    targ_Ss_means = np.nanmean(targ_frags[:, :, 0], axis=0)
    targ_Ss_quantiles = np.nanquantile(targ_frags[:, :, 0], q=q, axis=0)
    
    rand_s_means_means = np.nanmean(rand_frags[:, :, 1], axis=0)
    rand_s_means_quantiles = np.nanquantile(rand_frags[:, :, 1], q=q, axis=0)
    
    targ_s_means_means = np.nanmean(targ_frags[:, :, 1], axis=0)
    targ_s_means_quantiles = np.nanquantile(targ_frags[:, :, 1], q=q, axis=0)
    
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 8))
    fig.suptitle(f'{gtype} Fragmentation Plots')
    
    offset = 0.0005  # for better plotting so symbols don't overlap too much

    ## plot single graph + permutation combination
    ax = axs[0]
    
    ax.plot(fracs, rand_Ss, 'bs',
            linestyle='', markerfacecolor='none',
            label='S Random Removal')

    ax.plot(fracs, rand_s_means, 'bs',
            linestyle='',
            label='<s> Random Removal')

    ax.plot(fracs+offset, targ_Ss, 'ro',
            linestyle='', markerfacecolor='none',
            label='S Targeted Removal')

    ax.plot(fracs+offset, targ_s_means, 'ro',
            linestyle='',
            label='<s> Targeted Removal')

    ax.set_title(f'{gtype} Single Graph + Permutation Combination')
    ax.set_xlabel('Fraction of total nodes removed')
    ax.set_ylabel('S and <s>')
    ax.legend()

    ## plot mean + quantile of multiple graph + permutation combinations
    ax = axs[1]


    ax.plot(fracs, rand_Ss_means, 'bs',
            linestyle='', markerfacecolor='none')
    ax.fill_between(fracs, rand_Ss_quantiles[0], rand_Ss_quantiles[1],
                    color='blue', alpha=0.15)

    ax.plot(fracs, rand_s_means_means, 'bs',
            linestyle='')
    ax.fill_between(fracs, rand_s_means_quantiles[0], rand_s_means_quantiles[1],
                    color='blue', alpha=0.15)

    ax.plot(fracs+offset, targ_Ss_means, 'ro',
            linestyle='', markerfacecolor='none')
    ax.fill_between(fracs+offset, targ_Ss_quantiles[0], targ_Ss_quantiles[1],
                    color='red', alpha=0.15)

    ax.plot(fracs+offset, targ_s_means_means, 'ro',
            linestyle='')
    ax.fill_between(fracs+offset, targ_s_means_quantiles[0], targ_s_means_quantiles[1],
                    color='red', alpha=0.15)

    ax.set_title(f'{gtype} Mean and [{q[0]}, {q[1]}] Quantiles of\nMultiple Graph + Permutation Combinations')
    ax.set_xlabel('Fraction of total nodes removed')
    ax.set_ylabel('S and <s>')
    
    fig.savefig(f'{outdir}/plots/nx_{gtype.lower()}_frags.png')
    plt.clf()
    
    
def nx_plot_multilayer_synth_diam_results(rand_diams, targ_diams, fracs, graphtype, outdir):
    q = [0.025, 0.975]  # quantiles to compute
    
    fig, axs = plt.subplots(ncols=2, figsize=(12, 6), sharey=True)
    offset = 0.0005  # for better plotting so symbols don't overlap too much
    
    _diam_max = np.max((rand_diams.max(), targ_diams.max()))
    
    diams_list = [rand_diams, targ_diams]
    labels = ['Random Removal', 'Targeted Removal']
    markers = ['b*', 'rX']
    
    for i, (diams, label, marker) in enumerate(zip(diams_list, labels, markers)):
        if i % 2 == 0:  # rands
            _fracs = fracs
            color = 'blue'
        else:
            _fracs = fracs + offset
            color = 'red'
        
        diam_means = np.nanmean(diams, axis=0)
        diam_quantiles = np.nanquantile(diams, q=q, axis=0)
        
        ## plot single graph + permutation combination
        axs[0].plot(_fracs, diams[0, :], marker,
                    linestyle='', markerfacecolor='none',
                    label=f'{graphtype.upper()} {label}')
        axs[0].set_title('Single Graph + Permutation Combination')
        axs[0].set_xlabel('Fraction of total nodes removed')
        axs[0].set_ylabel('Diameter')
        axs[0].set_yticks(ticks=np.arange(_diam_max))
        axs[0].legend()
        
        ## plot mean + quantile results using all graph + idx permutation combinations
        axs[1].plot(_fracs, diam_means, marker, linestyle='', markerfacecolor='none')
        axs[1].fill_between(fracs, diam_quantiles[0], diam_quantiles[1],
                            color=color, alpha=0.15)
        axs[1].set_title(f'Mean and [{q[0]}, {q[1]}] Quantiles of\nMultiple Graph + Permutation Combinations')
        axs[1].set_xlabel('Fraction of total nodes removed')
        axs[1].set_ylabel('Diameter')
        axs[1].set_yticks(ticks=np.arange(_diam_max))
        
    fig.savefig(f'{outdir}/plots/nx_{graphtype}_diams.png')
    plt.clf()
    
    
def nx_plot_multilayer_synth_frag_results(rand_frags, targ_frags, fracs, graphtype, outdir):
    q = [0.025, 0.975]
    
    fig, axs = plt.subplots(ncols=2, figsize=(12, 6), sharey=True)
    offset = 0.0005  # for better plotting so symbols don't overlap too much
    
    datanames = ['Single Graph + Permutation Combination',
                 f'Mean and [{q[0]}, {q[1]}] Quantiles of\nMultiple Graph + Permutation Combinations']
    rem_strats = ['Random Removal', 'Targeted Removal']
    stat_names = ['S', '<s>']
    markers = ['b*', 'rX']
    colors = ['blue', 'red']
    
    frag_stats = [rand_frags, targ_frags]
    
    for i, (frags, rem_strat, marker) in enumerate(zip(frag_stats, rem_strats, markers)):
        if i % 2 == 0:  # rands
            _fracs = fracs
            color = 'blue'
        else:
            _fracs = fracs + offset
            color = 'red'
            
        for l, stat_name in enumerate(stat_names):
            label = f'{stat_name} {rem_strat}'
            mfc = 'none' if l == 0 else marker[0]
            
            frag_stat_means = np.nanmean(frags[:, :, l], axis=0)
            frag_stat_quantiles = np.nanquantile(frags[:, :, l], q=q, axis=0)
            
            axs[0].plot(_fracs, frags[0, :, l], marker, linestyle='',
                        markerfacecolor=mfc, label=label)
            
            axs[1].plot(_fracs, frag_stat_means, marker,
                        linestyle='', markerfacecolor=mfc)
            axs[1].fill_between(_fracs, frag_stat_quantiles[0], frag_stat_quantiles[1],
                                color=color, alpha=0.15)
                        
    axs[0].set_title(f'{graphtype.upper()} {datanames[0]}')
    axs[0].set_xlabel('Fraction of total nodes removed')
    axs[0].set_ylabel('S and <s>')
    axs[0].legend()

    axs[1].set_title(f'{graphtype.upper()} {datanames[1]}')
    axs[1].set_xlabel('Fraction of total nodes removed')
    axs[1].set_ylabel('S and <s>')
    
    fig.savefig(f'{outdir}/plots/nx_{graphtype}_frags.png')
    plt.clf()
    
    
def main():
    import argparse
    import os
    
    outdir = 'outfiles'
    os.makedirs(f'{outdir}/plots', exist_ok=True)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--graphtype', '-g', choices=['synth1', 'synth2', 'int', 'www'])
    parser.add_argument('--stat', '-s', choices=['diam', 'frag', 'degdist'])
    args = parser.parse_args()
    
    graphtype = args.graphtype
    stat = args.stat
    
    if stat == 'degdist':
        if graphtype == 'synth1':
            from graph_saveload import nx_load_synth_graphs
            gdict = nx_load_synth_graphs()
            gs = [gdict['ER'][0], gdict['PA'][0]]
            gtypes = ['er', 'pa']
            
        elif graphtype == 'synth2':
            from graph_saveload import nx_load_synth_multilayer_graphs
            gdict = nx_load_synth_multilayer_graphs()
            gs = [gdict['ER2'][0], gdict['PA2'][0], gdict['ERPA'][0]]
            gtypes = ['er2', 'pa2', 'erpa']
            
        else:
            from graph_saveload import nx_load_real_graphs
            gdict = nx_load_real_graphs()
            gs = gdict[graphtype.upper()][0]
            gtypes = [graphtype]
        
        plot_all_deg_dist(gs, gtypes, outdir)
    
    if graphtype == 'synth1' and stat == 'diam':
        # load calculated stats
        er_rand_diams = np.load(f'{outdir}/results/er_rand_diam.npy')
        er_targ_diams = np.load(f'{outdir}/results/er_targ_diam.npy')
        pa_rand_diams = np.load(f'{outdir}/results/pa_rand_diam.npy')
        pa_targ_diams = np.load(f'{outdir}/results/pa_targ_diam.npy')
        
        f_start = 0
        f_end = 0.05
        f_num = 21
        # f_num = 6
        fracs = np.linspace(f_start, f_end, f_num)
        
        plot_synth_diam_results(er_rand_diams, er_targ_diams, pa_rand_diams, pa_targ_diams, fracs, outdir)
    
    elif graphtype == 'synth1' and stat == 'frag':
        # load calculated stats
        er_rand_frags = np.load(f'{outdir}/results/er_rand_frag.npy')
        er_targ_frags = np.load(f'{outdir}/results/er_targ_frag.npy')
        pa_rand_frags = np.load(f'{outdir}/results/pa_rand_frag.npy')
        pa_targ_frags = np.load(f'{outdir}/results/pa_targ_frag.npy')
        
        f_start = 0
        f_end = 0.5
        f_num = 26
        # f_num = 6  ## ONLY DURING TESTING
        fracs = np.linspace(f_start, f_end, f_num)
        
        plot_synth_frag_results(er_rand_frags, er_targ_frags, pa_rand_frags, pa_targ_frags, fracs, outdir)
        
    elif graphtype == 'synth2' and stat == 'diam':
        f_start = 0
        f_end = 0.05
        f_num = 21
        # f_num = 6
        fracs = np.linspace(f_start, f_end, f_num)
        
        gtypes = ['er2', 'pa2', 'erpa']
        for gtype in gtypes:
            rand_diams = np.load(f'{outdir}/results/{gtype}_rand_diam.npy')
            targ_diams = np.load(f'{outdir}/results/{gtype}_targ_diam.npy')
            nx_plot_multilayer_synth_diam_results(rand_diams, targ_diams, fracs, gtype, outdir)
    
    elif graphtype == 'synth2' and stat == 'frag':
        f_start = 0
        f_end = 0.5
        f_num = 26
        # f_num = 6
        fracs = np.linspace(f_start, f_end, f_num)
        
        gtypes = ['er2', 'pa2', 'erpa']
        for gtype in gtypes:
            rand_frags = np.load(f'{outdir}/results/{gtype}_rand_frag.npy')
            targ_frags = np.load(f'{outdir}/results/{gtype}_targ_frag.npy')
            nx_plot_multilayer_synth_frag_results(rand_frags, targ_frags, fracs, gtype, outdir)
    
    elif (graphtype == 'int' or graphtype == 'www') and stat == 'diam':
        f_start = 0
        f_end = 0.025
        f_num = 21
        # f_num = 6
        fracs = np.linspace(f_start, f_end, f_num)
        
        rand_diams = np.load(f'{outdir}/results/{graphtype}_rand_diam.npy')
        targ_diams = np.load(f'{outdir}/results/{graphtype}_targ_diam.npy')
        
        plot_real_diam_results(rand_diams, targ_diams, fracs, graphtype, outdir)
    
    elif (graphtype == 'int' or graphtype == 'www') and stat == 'frag':
        f_start = 0
        f_end = 0.15
        f_num = 26
        # f_num = 6
        fracs = np.linspace(f_start, f_end, f_num)
        
        rand_frags = np.load(f'{outdir}/results/{graphtype}_rand_frag.npy')
        targ_frags = np.load(f'{outdir}/results/{graphtype}_targ_frag.npy')
        
        plot_real_frag_results(rand_frags, targ_frags, fracs, graphtype, outdir)


if __name__ == '__main__':
    main()
