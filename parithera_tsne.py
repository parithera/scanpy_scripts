
import scanpy as sc
import anndata as ad
import sys
import os
import json
import hdf5plugin

class TSNE:
    def __init__(self, ouput_path):
        self.output_path = ouput_path

    def run(self):
        sc.settings.figdir = self.output_path

        sc.settings.set_figure_params(dpi=50, facecolor="white")
        adata = ad.io.read_h5ad(self.output_path.replace("python", "out.h5ad"))
        
        # Plot t-SNE
        sc.pl.tsne(adata, color='sample', color_map='PuRd', size=50, save='graph.png', show=False)

        # Extract UMAP coordinates
        umap_data = adata.obsm['X_tsne']

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
            "type": "tsne"
        }

        # Save the UMAP data as a JSON file
        json_output_path = os.path.join(self.output_path, "tsne_data.json")
        with open(json_output_path, 'w') as f:
            json.dump(output, f, indent=4)

        return output

if __name__=='__main__':
    output_path = sys.argv[1]
    tsne = TSNE(output_path)
    tsne.run()
