
import scanpy as sc
import anndata as ad
import sys
import os
import json
import hdf5plugin

class Markers:
    def __init__(self, ouput_path, cluster_name):
        self.output_path = ouput_path
        self.cluster_name = cluster_name

    def run(self):
        sc.settings.figdir = self.output_path

        sc.settings.set_figure_params(dpi=50, facecolor="white")
        adata = ad.io.read_h5ad(self.output_path.replace("python", "out.h5ad"))
        
        # Obtain cluster-specific differentially expressed genes
        sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")
        sc.get.rank_genes_groups_df(adata, group=self.cluster_name).head(5)
        dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group=self.cluster_name).head(5)["names"]
        sc.pl.umap(
            adata,
            color=[*dc_cluster_genes, "leiden_res_0.50"],
            legend_loc="on data",
            frameon=False,
            ncols=3,
            save="graph.png", show=False
        )

        # Extract UMAP coordinates
        umap_data = adata.obsm['X_umap']

        umaps = {}

        for marker in dc_cluster_genes.tolist():
            # Create a list of dictionaries for each cell, containing its UMAP coordinates and any other metadata you want to include
            cells = []
            if 'leiden' in adata.obs:  # Example of including cluster labels (assuming leiden clustering was done)
                clusters = adata.obs['leiden'].astype(str).tolist()
            else:
                clusters = ['None'] * len(umap_data)

            samples = adata.obs['sample'].astype(str).tolist()
            for i, (x, y) in enumerate(umap_data):
                cell_info = {
                    'x': float(x),
                    'y': float(y),
                    'marker_expression': float(adata.X[i, adata.var_names.tolist().index(marker)]),     # Extract the expression level for the marker
                    'cluster': clusters[i],  # Add any other metadata you need here
                    'sample': samples[i]     # Extract the sample from which the data is taken
                }
                cells.append(cell_info)
            umaps[marker] = cells

        # Export the data to generate this plot to a JSON file
        export_data = {
            "cluster_genes": dc_cluster_genes.tolist(),
            "umaps": umaps
        }

        output = {
            "data": export_data,
            "type": "marker"
        }

        # Create the output directory if it does not exist
        os.makedirs(self.output_path, exist_ok=True)

        json_file_path = os.path.join(self.output_path, "marker_genes_data.json")
        with open(json_file_path, 'w') as f:
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

    cluster_name = arguments["cluster_name"]
    markers = Markers(output_path, cluster_name)
    markers.run()
