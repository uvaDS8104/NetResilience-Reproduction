# Network Science Final Project -- Network Resilience
This repository contains a reproduction of certain results from [Error and attack tolerance of complex networks (2000)](https://www.nature.com/articles/35019019) by Albert, Jeong, and Barabasi.
Specifically, this work reproduces the results in Figure 2 (diameter as a function of fraction of nodes removed) and Figure 3 (size of the largest component and average size of all the other components as a function of fraction of nodes removed).
Note that the results for the exponential (Erdos-Renyi) networks were updated in a [future correction (2001)](https://www.nature.com/articles/35054111).

## The Code
There are three main files of interest. The first is `netres.py` which handles the computation of the statistics and saves the result in `outfiles/results/`.
The second is `plotters.py` which handles, as the name suggests, the plotting of the results to `outfiles/plots/`.
The third is `graph_saveload.py` which generated the graphs to compute the results on, saved in `data/<graphtype>/`.

## The Extension
This project extends the paper by exploring how the plots look for three variations of 2-layer random graphs where each layer has the same number of nodes.
The first has both layers be Erdos-Renyi graphs.
The second has both layers be scale-free graphs.
The third connects an Erdos-Renyi and scale-free graph.
The linking between the two layers is done by creating a one-to-one mapping of nodes between the layers, and selecting a random subset of those pairs.

