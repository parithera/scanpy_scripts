
import scanpy as sc
import anndata as ad
import sys
import os
import json
import hdf5plugin

class Leiden:
    def __init__(self, ouput_path, leiden_res):
        self.output_path = ouput_path
        self.leiden_res = leiden_res

    def run(self):
        sc.settings.figdir = self.output_path

        sc.settings.set_figure_params(dpi=50, facecolor="white")
        adata = ad.io.read_h5ad(self.output_path.replace("python", "out.h5ad"))
        
        # Plot clusters
        sc.pl.umap(
            adata,
            color=["leiden_res_"+self.leiden_res],
            legend_loc="on data",
            save='graph.png', show=False
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
            "data": cells,
            "type": "leiden"
        }

        # Save the UMAP data as a JSON file
        json_output_path = os.path.join(self.output_path, "leiden_data.json")
        with open(json_output_path, 'w') as f:
            json.dump(output, f)

        return output

if __name__=='__main__':
    output_path = sys.argv[1]
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Open and parse the arguments.json file
    arguments_file_path = os.path.join(script_dir, "arguments.json")
    with open(arguments_file_path, 'r') as f:
        arguments = json.load(f)
    leiden_res =  arguments["leiden_res"]
    leiden = Leiden(output_path, leiden_res)
    leiden.run()
