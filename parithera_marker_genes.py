
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
    sc.pl.rank_genes_groups_dotplot(
        adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5, save="graph.png"
    )

