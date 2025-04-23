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


data_folder='../data'
use_RGB=False # set to False if you do not want to use RGB as features
spot_umi_threshold=50 # filter out cells with low UMIs


#Load csv data
coords = pd.read_csv("../../data/mucosa_sub_coords.csv")
coords = coords.iloc[:, 1:]  # Remove the first column
np.save('../data/coords_mat.npy', coords)

print(coords.shape)

gene_expr = pd.read_csv("../../data/counts_sub_mucosa.csv")
#gene_expr = pd.read_csv("../data/mucosa_counts.csv")
gene_expr = gene_expr.iloc[:, 1:]  # Remove the first column
#Make sure gene_expr is int
#gene_expr = gene_expr.astype(int)
np.save('../data/counts_mat.npy', gene_expr)

print(gene_expr.shape)

# Load gene labels
gene_labels = pd.read_csv("../data/mucosa_gene_names.csv")
gene_labels = gene_labels['x'].values
np.save("../data/gene_labels.npy", gene_labels)

# save matrices
#np.save('../data/counts_mat.npy', counts_mat)
#np.save('../data/coords_mat.npy', coords_mat)
#np.save('colorectal_tumor_data/gene_labels.npy', gene_labels)


counts_mat = np.load('../data/counts_mat.npy', allow_pickle=True)
coords_mat = np.load('../data/coords_mat.npy', allow_pickle=True)
gene_labels = np.load('../data/gene_labels.npy', allow_pickle=True)

num_dims=5
penalty=20 # may need to increase if this is too small

# CHANGE THESE PARAMETERS TO REDUCE RUNTIME
num_iters=30
eps=1e-4
num_genes=counts_mat.shape[0]

counts_mat = counts_mat.T
counts_mat_glmpca=counts_mat[:,np.argsort(np.sum(counts_mat, axis=0))[-num_genes:]]
glmpca_res=glmpca.glmpca(counts_mat_glmpca.T, num_dims, fam="poi", penalty=penalty, verbose=True,
                        ctl = {"maxIter":num_iters, "eps":eps, "optimizeTheta":True})
A = glmpca_res['factors'] # should be of size N x num_dims, where each column is a PC

if use_RGB:
    A=np.hstack((A,rgb_mean)) # attach to RGB mean
np.save('../data/glmpca.npy', A)

print(A.shape)
rotated_coords=dp_related.rotate_by_theta(coords_mat, 0)
R=2
C=2
fig,axs=plt.subplots(R,C,figsize=(20,10))
for r in range(R):
    for c in range(C):
        i=r*C+c
        axs[r,c].scatter(rotated_coords[:,0], rotated_coords[:,1], c=A[:,i],cmap='Reds',s=3)
        if i < num_dims:
            axs[r,c].set_title(f'GLM-PC{i}')
        else:
            axs[r,c].set_title('RGB'[i-num_dims])

plt.tight_layout()  # Adjust the layout
plt.savefig('../plots/glmpca_components.png', dpi=300, bbox_inches='tight')
plt.close()  # Close the figure to free memory


S = coords_mat 

# z-score normalize S and A
S_torch, A_torch = neural_net.load_rescale_input_data(S,A)




# NEURAL NET PARAMETERS (USER CAN CHANGE)
# architectures are encoded as list, eg [20,20] means two hidden layers of size 20 hidden neurons
isodepth_arch=[20,20] # architecture for isodepth neural network d(x,y) : R^2 -> R 
expression_fn_arch=[20,20] # architecture for 1-D expression function h(w) : R -> R^G

num_epochs = 10000 # number of epochs to train NN (NOTE: it is sometimes beneficial to train longer)
checkpoint = 500 # save model after number of epochs = multiple of checkpoint
out_dir='colorectal_tumor_tutorial_outputs' # folder to save model runs
optimizer = "adam"
num_restarts=3

######################################

seed_list=range(num_restarts)
for seed in seed_list:
    print(f'training neural network for seed {seed}')
    out_dir_seed=f"{out_dir}/rep{seed}"
    os.makedirs(out_dir_seed, exist_ok=True)
    mod, loss_list = neural_net.train(S_torch, A_torch,
                          S_hidden_list=isodepth_arch, A_hidden_list=expression_fn_arch, 
                          epochs=num_epochs, checkpoint=checkpoint, 
                          save_dir=out_dir_seed, optim=optimizer, seed=seed, save_final=True)
    

# gaston_model, A, S= process_NN_output.process_files('colorectal_tumor_tutorial_outputs') # model trained above
gaston_model, A, S= process_NN_output.process_files('./colorectal_tumor_data/rep1') # MATCH PAPER FIGURES

num_layers=5 # CHANGE FOR YOUR APPLICATION: use number of layers from above!
gaston_isodepth, gaston_labels=dp_related.get_isodepth_labels(gaston_model,A,S,num_layers)

# DATASET-SPECIFIC: so domains are ordered with tumor being last
gaston_isodepth= np.max(gaston_isodepth) -1 * gaston_isodepth
gaston_labels=(num_layers-1)-gaston_labels

show_streamlines=True
rotate = 0 # rotate coordinates by -90
arrowsize=2

cluster_plotting.plot_isodepth(gaston_isodepth, S, gaston_model, figsize=(7,6), streamlines=show_streamlines, 
                               rotate=rotate,arrowsize=arrowsize, 
                               neg_gradient=True)


plt.savefig('../plots/isodepth_plot.png', dpi=300, bbox_inches='tight')
plt.close()  # Close the figure to free memory