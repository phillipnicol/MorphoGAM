import numpy as np
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
import anndata as ad
import cv2
import spateo as st
from scipy.sparse import csr_matrix
from spateo.digitization.grid import digitize, gridit

# ============================================================
# Load CA3 data exported from R
# ============================================================

xy = pd.read_csv("../../moffitt_mucosa/data/ca3_xy.csv").to_numpy()

X = scipy.io.mmread(
    "../../moffitt_mucosa/data/ca3_counts.mtx"
).T.tocsr()

genes = pd.read_csv(
    "../../moffitt_mucosa/data/ca3_genes.txt",
    header=None
)[0].astype(str).values


adata = ad.AnnData(X=X)
adata.var_names = genes
adata.obsm["spatial"] = xy

# Required by Spateo's config/type checker
adata.uns["__type"] = "UMI"

# ============================================================
# User-specified Spateo boundary points
# ============================================================

pnt_xy  = np.array([3000, 3500])  # left lower
pnt_xY  = np.array([3000, 4000])  # left upper
pnt_Xy  = np.array([4800, 2000])  # right bottom
pnt_XY  = np.array([4800, 2500])  # right upper

print("\nChosen Spateo boundary points:\n")
print("pnt_xy  (left lower): ", pnt_xy)
print("pnt_xY  (left upper): ", pnt_xY)
print("pnt_Xy  (right bottom):", pnt_Xy)
print("pnt_XY  (right upper):", pnt_XY)

# ============================================================
# Build a binary image and extract contour
# ============================================================

pad = 100
max_x = int(np.ceil(max(xy[:, 0].max(), pnt_Xy[0], pnt_XY[0]))) + pad
max_y = int(np.ceil(max(xy[:, 1].max(), pnt_xY[1], pnt_XY[1]))) + pad

img = np.zeros((max_y + 1, max_x + 1), dtype=np.uint8)

# Draw each cell as a small disk to create a connected tissue mask.
# Increase radius if the contour is fragmented.
radius = 25

for x, y in xy:
    cv2.circle(
        img,
        (int(round(x)), int(round(y))),
        radius,
        255,
        -1
    )

# Morphological closing to fill small gaps.
kernel = np.ones((31, 31), np.uint8)
img_closed = cv2.morphologyEx(img, cv2.MORPH_CLOSE, kernel)

ctrs, hierarchy = cv2.findContours(
    img_closed,
    cv2.RETR_EXTERNAL,
    cv2.CHAIN_APPROX_NONE
)

if len(ctrs) == 0:
    raise RuntimeError("No contours found. Increase radius or check coordinates.")

areas = np.array([cv2.contourArea(c) for c in ctrs])
ctr_idx = int(np.argmax(areas))

print(f"\nFound {len(ctrs)} contours. Using largest contour: ctr_idx={ctr_idx}, area={areas[ctr_idx]:.1f}")

# ============================================================
# Plot contour and points
# ============================================================

fig, ax = plt.subplots(figsize=(7, 6))

ax.scatter(xy[:, 0], xy[:, 1], s=2, alpha=0.4)

contour = ctrs[ctr_idx][:, 0, :]
ax.plot(contour[:, 0], contour[:, 1], linewidth=1)

pts = np.vstack([pnt_xy, pnt_xY, pnt_Xy, pnt_XY])
labs = ["pnt_xy", "pnt_xY", "pnt_Xy", "pnt_XY"]

ax.scatter(pts[:, 0], pts[:, 1], s=120)

for lab, pt in zip(labs, pts):
    ax.text(pt[0], pt[1], lab, fontsize=11)

ax.set_aspect("equal")
ax.set_title("CA3 contour and Spateo boundary points")
plt.show()

# ============================================================
# Run Spateo digitization
# ============================================================

digitize(
    adata=adata,
    ctrs=ctrs,
    ctr_idx=ctr_idx,
    pnt_xy=pnt_xy,
    pnt_Xy=pnt_Xy,
    pnt_xY=pnt_xY,
    pnt_XY=pnt_XY,
    spatial_key="spatial",
    dgl_layer_key="digital_layer",
    dgl_column_key="digital_column",
    max_itr=100000,
)

gridit(
    adata=adata,
    layer_num=10,
    column_num=10,
    dgl_layer_key="digital_layer",
    dgl_column_key="digital_column",
    layer_label_key="layer_label",
    column_label_key="column_label",
    grid_label_key="grid_label",
)

# ============================================================
# Plot digitization values
# ============================================================

for col in ["digital_layer", "digital_column"]:
    fig, ax = plt.subplots(figsize=(7, 6))

    sc = ax.scatter(
        xy[:, 0],
        xy[:, 1],
        c=adata.obs[col].astype(float),
        s=3
    )

    ax.set_aspect("equal")
    ax.set_title(f"Spateo {col}")
    plt.colorbar(sc, ax=ax)
    plt.show()

# ============================================================
# Save output
# ============================================================

out = pd.DataFrame({
    "x": xy[:, 0],
    "y": xy[:, 1],
    "spateo_layer": np.asarray(adata.obs["digital_layer"]).astype(float),
    "spateo_column": np.asarray(adata.obs["digital_column"]).astype(float),
    "spateo_layer_label": np.asarray(adata.obs["layer_label"]),
    "spateo_column_label": np.asarray(adata.obs["column_label"]),
    "spateo_grid_label": np.asarray(adata.obs["grid_label"]),
})

out.to_csv("ca3_spateo_digitization_values.csv", index=False)

print("\nSaved: ca3_spateo_digitization_values.csv")