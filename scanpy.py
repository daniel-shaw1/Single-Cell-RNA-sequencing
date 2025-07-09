import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt


# === Load and label each dataset ===
adata1 = sc.read_10x_mtx('run1/filtered_feature_bc_matrix', var_names='gene_ids')
adata1.obs['batch'] = 'run1'

adata2 = sc.read_10x_mtx('run2/filtered_feature_bc_matrix', var_names='gene_ids')
adata2.obs['batch'] = 'run2'

adata3 = sc.read_10x_mtx('run3/filtered_feature_bc_matrix', var_names='gene_ids')
adata3.obs['batch'] = 'run3'

# === Combine datasets ===
adata = adata1.concatenate(adata2, adata3, batch_key='batch')

# === Save raw counts ===
adata.raw = adata

# === Preprocessing ===
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='batch')
sc.pl.highly_variable_genes(adata, save='_hvg.png', show=False)
plt.close()

adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

# === PCA and neighbors ===
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, save='_pca.png', show=False)
plt.close()

sc.pp.neighbors(adata)
sc.tl.umap(adata)

# === Clustering ===
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=['leiden', 'batch'], save='_clusters_batches.png', show=False)
plt.close()

# === Optional: Save marker genes heatmap or dotplot ===
# sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
# sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, save='_marker_genes.png', show=False)
