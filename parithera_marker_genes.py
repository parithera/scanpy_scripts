
import scanpy as sc
import anndata as ad
import sys
import hdf5plugin
import os
import json

if __name__=='__main__':
    output_path = sys.argv[1]
    sc.settings.figdir = output_path

    sc.settings.set_figure_params(dpi=50, facecolor="white")
    adata = ad.read_h5ad(output_path.replace("python", "out.h5ad"))
    
    # Obtain cluster-specific differentially expressed genes
    sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")
    sc.get.rank_genes_groups_df(adata, group="7").head(5)
    dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group="7").head(5)["names"]
    sc.pl.umap(
        adata,
        color=[*dc_cluster_genes, "leiden_res_0.50"],
        legend_loc="on data",
        frameon=False,
        ncols=3,
        save="graph.png"
    )

