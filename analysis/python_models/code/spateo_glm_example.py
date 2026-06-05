import numpy as np
import pandas as pd
import anndata as ad
import spateo as st

# ============================================================
# Minimal fake dataset with t and r coordinates
# ============================================================

np.random.seed(1)

n_cells = 300
n_genes = 6

# ------------------------------------------------------------
# Simulated morphologically relevant coordinates
# ------------------------------------------------------------

t = np.linspace(0, 1, n_cells)

r = np.random.uniform(-1, 1, n_cells)

# ------------------------------------------------------------
# Simulate counts
# ------------------------------------------------------------

Y = np.random.poisson(1, size=(n_cells, n_genes))

# Gene0 varies only with t
Y[:, 0] += np.random.poisson(
    10 * np.sin(np.pi * t) ** 2
)

# Gene1 varies linearly with t
Y[:, 1] += np.random.poisson(
    8 * t
)

# Gene2 varies only with r
Y[:, 2] += np.random.poisson(
    10 * (r ** 2)
)

# Gene3 varies with both t and r
Y[:, 3] += np.random.poisson(
    6 * np.sin(np.pi * t) ** 2
    + 6 * (r > 0)
)

# Gene4 and Gene5 are mostly noise

# ------------------------------------------------------------
# Create AnnData
# ------------------------------------------------------------

adata = ad.AnnData(X=Y)

adata.var_names = [
    "Gene_t1",
    "Gene_t2",
    "Gene_r",
    "Gene_both",
    "Noise1",
    "Noise2",
]

# Spateo covariates must be in obs
adata.obs["t"] = t
adata.obs["r"] = r

# Required metadata for Spateo
adata.uns["__type"] = "UMI"

# ============================================================
# Run Spateo GLM
# ============================================================

st.tl.glm.glm_degs(
    adata,
    fullModelFormulaStr="~ cr(t, df=3) + cr(r, df=3)",
    reducedModelFormulaStr="~ 1",
    qval_threshold=None,
    llf_threshold=None,
)

# ============================================================
# Inspect output
# ============================================================

print("\nKeys in adata.uns:\n")
print(adata.uns.keys())

print("\nTop genes:\n")

res = adata.uns["glm_degs"]["glm_result"].sort_values("qval").reset_index(names="gene")

print(
    res[
        ["gene", "pval", "qval"]
    ].head(10)
)

print("\nAll columns:\n")
print(res.columns)

# ============================================================
# Optional: inspect all results
# ============================================================

print("\nFull result table:\n")
print(res)
