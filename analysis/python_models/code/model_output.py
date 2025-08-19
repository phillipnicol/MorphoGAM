import sys
import os
from collections import defaultdict
import pandas as pd
import scanpy as sc
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from glmpca import glmpca
from itertools import combinations
import torch

import sys
from importlib import reload

import gaston
from gaston import neural_net,cluster_plotting, dp_related, segmented_fit, restrict_spots, model_selection
from gaston import binning_and_plotting, isodepth_scaling, run_slurm_scripts, parse_adata
from gaston import spatial_gene_classification, plot_cell_types, filter_genes, process_NN_output

import seaborn as sns
import math


counts_mat = np.load('../data/counts_mat.npy', allow_pickle=True).T
coords_mat = np.load('../data/coords_mat.npy', allow_pickle=True)
gene_labels = np.load('../data/gene_labels.npy', allow_pickle=True)

# gaston_model, A, S= process_NN_output.process_files('colorectal_tumor_tutorial_outputs') # model trained above
gaston_model, A, S= process_NN_output.process_files('colorectal_tumor_tutorial_outputs') # MATCH PAPER FIGURES

num_layers=1 # CHANGE FOR YOUR APPLICATION: use number of layers from above!
gaston_isodepth, gaston_labels=dp_related.get_isodepth_labels(gaston_model,A,S,num_layers)

# DATASET-SPECIFIC: so domains are ordered with tumor being last
gaston_isodepth= np.max(gaston_isodepth) -1 * gaston_isodepth
gaston_labels=(num_layers-1)-gaston_labels

print(gaston_isodepth)
#Save the isodepth values as csv
np.savetxt('../data/gaston_isodepth.csv', gaston_isodepth, delimiter=',')
#np.save('../data/gaston_isodepth.npy', gaston_isodepth)

show_streamlines=True
rotate = np.radians(-90) # rotate coordinates by -90
arrowsize=2

cluster_plotting.plot_isodepth(gaston_isodepth, S, gaston_model, figsize=(7,6), streamlines=show_streamlines, 
                               rotate=rotate,arrowsize=arrowsize, 
                               neg_gradient=True)


plt.savefig('../plots/isodepth_plot.png', dpi=300, bbox_inches='tight')
plt.close()  # Close the figure to free memory


umi_thresh = 10
#reload(segmented_fit)
# compute piecewise linear fit for restricted spots
pw_fit_dict=segmented_fit.pw_linear_fit(counts_mat, gaston_labels, gaston_isodepth,
                                        None, [],  umi_threshold=umi_thresh, isodepth_mult_factor=0.01,)
# for plotting

print("Here")

binning_output=binning_and_plotting.bin_data(counts_mat, gaston_labels, gaston_isodepth, 
                         None, gene_labels, num_bins=50, umi_threshold=umi_thresh)

domain_colors=['dodgerblue', '#F44E3F']

discont_genes_layer=spatial_gene_classification.get_discont_genes(pw_fit_dict, binning_output,q=0.95)
cont_genes_layer=spatial_gene_classification.get_cont_genes(pw_fit_dict, binning_output,q=0.8)

print(discont_genes_layer)
print(cont_genes_layer)

gene_name='Ddx58'
print(f'gene {gene_name}: discontinuous after domain(s) {discont_genes_layer[gene_name]}') 
print(f'gene {gene_name}: continuous in domain(s) {cont_genes_layer[gene_name]}')

# display log CPM (if you want to do CP500, set offset=500)
offset=10**6

binning_and_plotting.plot_gene_pwlinear(gene_name, pw_fit_dict, gaston_labels, gaston_isodepth, 
                                        binning_output, cell_type_list=None, pt_size=50, 
                                        linear_fit=True, ticksize=15, figsize=(4,2.5), offset=offset, lw=3,
                                       domain_boundary_plotting=True)

plt.savefig('../plots/dhx58_gene.png', dpi=300, bbox_inches='tight')
plt.close()  # Close the figure to free memory


print("Through") 

gene_name='Apob'
print(f'gene {gene_name}: discontinuous after domain(s) {discont_genes_layer[gene_name]}') 
print(f'gene {gene_name}: continuous in domain(s) {cont_genes_layer[gene_name]}')

# display log CPM (if you want to do CP500, set offset=500)
offset=10**6

binning_and_plotting.plot_gene_pwlinear(gene_name, pw_fit_dict, gaston_labels, gaston_isodepth, 
                                        binning_output, cell_type_list=None, pt_size=50, 
                                        linear_fit=True, ticksize=15, figsize=(4,2.5), offset=offset, lw=3,
                                       domain_boundary_plotting=True)

plt.savefig('../plots/epob_gene.png', dpi=300, bbox_inches='tight')
plt.close()  # Close the figure to free memory
