
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
    
    # Generating UMAP plot
    sc.pl.umap(adata, color='sample', save='graph.png', color_map="PuRd", size=50)

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
        "type": "umap"
    }

    # Save the UMAP data as a JSON file
    json_output_path = os.path.join(output_path, "umap_data.json")
    with open(json_output_path, 'w') as f:
        json.dump(output, f, indent=4)

