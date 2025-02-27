
import scanpy as sc
import anndata as ad
import sys
import os
import json
import hdf5plugin

if __name__=='__main__':
    output_path = sys.argv[1]
    sc.settings.figdir = output_path

    sc.settings.set_figure_params(dpi=50, facecolor="white")
    adata = ad.io.read_h5ad(output_path.replace("python", "out.h5ad"))
    
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

    # Export the data to generate this plot to a JSON file
    export_data = {
        "cluster_genes": dc_cluster_genes.tolist(),
        "adata_metadata": {
            "obs_names": adata.obs_names.tolist(),
            "var_names": adata.var_names.tolist(),
            "uns_keys": list(adata.uns.keys()),
            "obsm_keys": list(adata.obsm.keys())
        }
    }

    # Create the output directory if it does not exist
    os.makedirs(output_path, exist_ok=True)

    json_file_path = os.path.join(output_path, "exported_data.json")
    with open(json_file_path, 'w') as f:
        json.dump(export_data, f, indent=4)

