
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
    
    # Plot clusters
    sc.pl.umap(
        adata,
        color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
        wspace=0.5,
        ncols=2,
        legend_loc="on data",
        save='graph.png'
    )

    # Extract UMAP coordinates
    umap_data = adata.obsm['X_umap']

    # Create a list of dictionaries for each cell, containing its UMAP coordinates and any other metadata you want to include
    cells = []

    if 'leiden' in adata.obs:  # Example of including cluster labels (assuming leiden clustering was done)
        clusters = adata.obs['leiden'].astype(str).tolist()
    else:
        clusters = ['None'] * len(umap_data)

    samples = adata.obs['sample'].astype(str).tolist()

    for i, (x, y) in enumerate(umap_data):
        cell_info = {
            'id': i,
            'x': float(x),
            'y': float(y),
            'cluster': clusters[i],  # Add any other metadata you need here
            'sample': samples[i]     # Extract the sample from which the data is taken
        }
        cells.append(cell_info)

    output = {
        "cells": cells,
        "type": "cluster"
    }

    # Save the UMAP data as a JSON file
    json_output_path = os.path.join(output_path, "cluster_data.json")
    with open(json_output_path, 'w') as f:
        json.dump(output, f, indent=4)

