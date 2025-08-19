import argparse
import os
import numpy as np
import pandas as pd
import anndata as ad

# Spateo (only) imports for graph + digitization + plotting
import spateo as st
from spateo.tools import find_neighbors as fn          # neighbor ops
from spateo.digitization import utils as digi_utils    # digitize_general
# plotting helpers (we'll try spateo, but also ship a matplotlib fallback)
from spateo.plotting.static import scatters as st_scatters

import matplotlib.pyplot as plt

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coords", default="../data/coords_mat.npy", 
                    help="Path to coordinates .npy file.")
    ap.add_argument("--axis", choices=["x", "y"], default="y",
                    help="Axis for auto-boundary selection (percentile).")
    ap.add_argument("--lower_pct", type=float, default=2.0,
                    help="Lower percentile for auto boundary (0-100).")
    ap.add_argument("--upper_pct", type=float, default=98.0,
                    help="Upper percentile for auto boundary (0-100).")
    ap.add_argument("--k", type=int, default=8,
                    help="Number of neighbors for KNN graph.")
    ap.add_argument("--outdir", default="../plots/spateo_output", 
                    help="Output directory.")
    ap.add_argument("--point_size", type=float, default=4.0, 
                    help="Marker size.")
    ap.add_argument("--lower_idx", type=int, nargs='+', default=None,
                    help="Explicit lower boundary indices (optional).")
    ap.add_argument("--upper_idx", type=int, nargs='+', default=None,
                    help="Explicit upper boundary indices (optional).")
    return ap.parse_args()

def build_knn_adjacency(coords: np.ndarray, k: int) -> np.ndarray:
    """
    Build a symmetric KNN adjacency using spateo.tools.find_neighbors.
    """
    # Use direct KNN calculation instead of affinity matrix
    try:
        # Try spateo's built-in KNN function first
        knn_idx, knn_w = fn.find_neighbors(coords, n_neighbors=k)
        adj_sparse = fn.knn_to_adj(knn_idx, knn_w)
        adj = adj_sparse.maximum(adj_sparse.T).toarray()
        return fn.normalize_adj(adj, exclude_self=True)
    except Exception as e:
        print(f"Warning: Spateo KNN failed ({str(e)}), falling back to sklearn")
        # Fallback to sklearn
        from sklearn.neighbors import NearestNeighbors
        nbrs = NearestNeighbors(n_neighbors=k, algorithm='ball_tree').fit(coords)
        distances, indices = nbrs.kneighbors(coords)
        
        # Convert to adjacency matrix
        n = coords.shape[0]
        adj = np.zeros((n, n))
        for i in range(n):
            adj[i, indices[i]] = 1
            adj[indices[i], i] = 1
        
        return fn.normalize_adj(adj, exclude_self=True)

def choose_boundaries(coords: np.ndarray,
                      axis: str,
                      lower_pct: float,
                      upper_pct: float,
                      lower_idx=None,
                      upper_idx=None):
    """
    Pick boundary indices either from explicit lists or by percentiles along chosen axis.
    """
    n = coords.shape[0]
    if lower_idx is not None and upper_idx is not None:
        lower = np.array(sorted(set([i for i in lower_idx if 0 <= i < n])), dtype=int)
        upper = np.array(sorted(set([i for i in upper_idx if 0 <= i < n])), dtype=int)
        if len(lower) == 0 or len(upper) == 0:
            raise ValueError("Provided boundary index lists are empty after sanity checks.")
        return lower, upper

    ax = 0 if axis == "x" else 1
    vals = coords[:, ax]
    lo_thr = np.percentile(vals, lower_pct)
    hi_thr = np.percentile(vals, upper_pct)
    lower = np.where(vals <= lo_thr)[0]
    upper = np.where(vals >= hi_thr)[0]
    # Safety: make sure we have at least one index on each side
    if len(lower) == 0:
        lower = np.array([int(np.argmin(vals))])
    if len(upper) == 0:
        upper = np.array([int(np.argmax(vals))])
    return lower, upper

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load coordinates from .npy file
    print("Loading coordinates...")
    coords = np.load(args.coords, allow_pickle=True)
    print(f"Processing {coords.shape[0]} points...")

    # Downsample if dataset is too large
    if coords.shape[0] > 10000:
        print(f"Large dataset detected ({coords.shape[0]} points). Downsampling...")
        idx = np.random.choice(coords.shape[0], size=10000, replace=False)
        coords = idx

    # Reduce memory usage in adjacency calculation
    print("Building KNN graph...")
    adj = build_knn_adjacency(coords, k=min(args.k, 20))  # Limit k to reduce memory

    # 2) Make an AnnData just to keep everything tidy
    adata = ad.AnnData(X=np.zeros((coords.shape[0], 1)))
    adata.obsm["spatial"] = coords
    adata.obs.index = [f"cell_{i}" for i in range(coords.shape[0])]

    # 3) Neighbor graph (Spateo-only)
    adj = build_knn_adjacency(coords, k=args.k)

    # 4) Get boundaries for both axes
    lower_idx, upper_idx = choose_boundaries(
        coords,
        axis='y',
        lower_pct=args.lower_pct,
        upper_pct=args.upper_pct,
        lower_idx=args.lower_idx,
        upper_idx=args.upper_idx
    )
        
    lower_idx_x, upper_idx_x = choose_boundaries(
        coords,
        axis='x',
        lower_pct=args.lower_pct,
        upper_pct=args.upper_pct,
        lower_idx=None,  # Add these as command line args if needed
        upper_idx=None
    )

    # Perform both layer (y-axis) and column (x-axis) digitization
    print("Performing layer digitization...")
    layer_heat = digi_utils.digitize_general(coords,
        adj, 
        boundary_lower=lower_idx,
        boundary_upper=upper_idx
    )
    
    print("Performing column digitization...")
    column_heat = digi_utils.digitize_general(coords,
        adj,
        boundary_lower=lower_idx_x,
        boundary_upper=upper_idx_x
    )

    # Plot and save layer digitization separately
    plt.figure(figsize=(7, 6))
    scatter1 = plt.scatter(coords[:, 0], coords[:, 1], c=layer_heat, 
                          s=args.point_size, cmap='viridis')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title('Layer Digitization')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.colorbar(scatter1)
    plt.tight_layout()
    
    # Save layer digitization plot
    layer_png = os.path.join(args.outdir, 'layer_digitization.png')
    plt.savefig(layer_png, dpi=250)
    plt.close()

    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    
    # Plot layer digitization
    scatter1 = ax1.scatter(coords[:, 0], coords[:, 1], c=layer_heat, 
                          s=args.point_size, cmap='viridis')
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_title('Layer Digitization')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    plt.colorbar(scatter1, ax=ax1)

    # Plot column digitization
    scatter2 = ax2.scatter(coords[:, 0], coords[:, 1], c=column_heat, 
                          s=args.point_size, cmap='viridis')
    ax2.set_aspect('equal', adjustable='box')
    ax2.set_title('Column Digitization')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    plt.colorbar(scatter2, ax=ax2)

    plt.tight_layout()
    
    # Save the plot
    png_out = os.path.join(args.outdir, 'digitization_comparison.png')
    plt.savefig(png_out, dpi=250)
    plt.close()

    # Save digitization values
    df_out = pd.DataFrame({
        'x': coords[:, 0],
        'y': coords[:, 1],
        'layer_digitization': layer_heat,
        'column_digitization': column_heat
    })
    csv_out = os.path.join(args.outdir, 'digitization_values.csv')
    df_out.to_csv(csv_out, index=False)

    print(f"[OK] Saved:\n - {csv_out}\n - {png_out}\n - {layer_png}")

if __name__ == "__main__":
    main()
