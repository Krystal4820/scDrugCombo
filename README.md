# scDrugCombo: Predicting Effective Anti-Cancer Drug Combinations Using Patient-Derived Single-Cell RNA-Sequencing Data

_scDrugCombo is the updated version of [scDrug: From scRNA-seq to Drug Repositioning](https://github.com/ailabstw/scDrug) with a new function in the third part **Drug Combination Prediction** of the workflow, Utilizing the powerful potential of scRNA-seq in tumor heterogeneity analysis, combined with knowledge of drug MoA and machine learning methods, to efficiently screen anticancer drugs and predict responses, ultimately identifying effective anticancer drug combinations._

The scDrugCombo constructed a workflow for comprehensive analysis on single-cell RNA sequencing (scRNA-seq) data. It provided a powerful tool with various functions, from fundamental data analysis to drug response prediction, and treatment suggestions.

The scDrugCombo went through three parts on raw scRNA-seq data investigation: **Single-Cell Data Analysis**, **Drug Response Prediction**, and **Drug Combination Prediction**.

- **Single-Cell Data Analysis** performed data preprocessing, clustering, cell type annotation, Gene Set Enrichment Analysis (GSEA), and survival analysis.

- **Drug Response Prediction** estimated the half maximal inhibitory concentration (IC50) of cell clusters, and reported the cell death percentages as drug kill efficacy.

- **Drug Combination Prediction** listed drug combinations with different predict clusters in the NMF Leiden clustering results and their synergy scores obtained from testing on cell lines.

## Download and Installation

1.  Clone the repository to local directory, e.g., `./scDrug`.

    ```
    git clone https://github.com/Krystal4820/scDrugCombo.git ./scDrug
    ```

2.  Build the Docker image tagged `sc-drug`.

    ```
    docker build -t sc-drug ./scDrug
    ```

3.  Run the Docker container named `scDrug` with `/docker/path` mounted to `/server/path` to access files within the Docker container.

    ```
    docker run -it --name scDrug -v /server/path:/docker/path --privileged sc-drug
    ```

## Usage

Note: Refer to [example](example) for a detail illustration of the usage for the scDrug.

### Single-Cell Data Analysis

**Single-Cell Data Analysis** took the scRNA-seq data in a 10x-Genomics-formatted mtx directory or a CSV file as input, performed fundamental data analysis, and output a Scanpy Anndata object `scanpyobj.h5ad`, a UMAP `umap_cluster.png` and differentially expressed genes (DEGs) `cluster_DEGs.csv` of the clustering result, and a gene expression profile (GEP) file `GEP.txt`.

Optionally, **Single-Cell Data Analysis** carried out batch correction, cell type annotation and Gene Set Enrichment Analysis (GSEA), and provided additional UMAPs showing batch effects and cell types (`umap_batch.png` and `umap_cell_type.png`), and the GSEA result `GSEA_results.csv`. For cell type annotation, we used [scMatch: a single-cell gene expression profile annotation tool using reference datasets](https://github.com/asrhou/scMatch).

Furthermore, **Single-Cell Data Analysis** could take previously produced Anndata as input and applied sub-clustering on specified clusters.

- Run `python3 single_cell_analysis.py -h` to show the help messages as follow for **Single-Cell Data Analysis**.

```
usage: single_cell_analysis.py [-h] -i INPUT [-f FORMAT] [-o OUTPUT] [-r RESOLUTION] [--impute] [--auto-resolution] [-m METADATA]
                               [-b BATCH] [-c CLUSTERS] [--cname CNAME] [--GEP] [--annotation] [--gsea] [--cpus CPUS] [--survival]
                               [--tcga TCGA] [--id ID] [--prefix PREFIX] [--not_treated]

scRNA-seq data analysis

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input 10x directory or CSV file
  -f FORMAT, --format FORMAT
                        input format, 10x (default) | csv | h5ad (Anndata object for subclustering with --clusters CLUSTERS)
  -o OUTPUT, --output OUTPUT
                        path to output directory, default='./'
  -r RESOLUTION, --resolution RESOLUTION
                        resolution for clustering, default=0.6
  --impute              do imputation. default: no
  --auto-resolution     automatically determine resolution for clustering
  -m METADATA, --metadata METADATA
                        path to metadata CSV file for batch correction (index as input in first column)
  -b BATCH, --batch BATCH
                        column in metadata (or adata.obs) for batch correction, e.g. 'PatientID'
  -c CLUSTERS, --clusters CLUSTERS
                        perform single cell analysis only on specified clusters, e.g. '1,3,8,9'
  --cname CNAME         which variable should be used when selecting clusters; required when clusters are provided. Default:
                        'louvain'
  --GEP                 generate Gene Expression Profile file.
  --annotation          perform cell type annotation
  --gsea                perform gene set enrichment analysis (GSEA)
  --cpus CPUS           number of CPU used for auto-resolution and annotation, default=1
  --survival            perform survival analysis
  --tcga TCGA           path to TCGA data
  --id ID               Specify TCGA project id in the format "TCGA-xxxx", e.g., "TCGA-LIHC"
  --prefix PREFIX       Any prefix before matrix.mtx, genes.tsv and barcodes.tsv.
  --not_treated         only consider untreated samples from TCGA for survival analysis.
```

- Apply **Single-Cell Data Analysis** with batch correction, clustering resolution 1.0, cell type annotation and GSEA.

```
python3 single_cell_analysis.py --input INPUT --metadata METADATA --batch BATCH --resolution 1.0 --annotation --gsea
```

- **Single-Cell Data Analysis** for sub-clustering on specified clusters at automatically determined resolution run under 4 cpus.

```
python3 single_cell_analysis.py -f h5ad --input scanpyobj.h5ad --clusters CLUSTERS --auto-resolution --cpus 4
```

### Drug Response Prediction

**Drug Response Prediction** examined `scanpyobj.h5ad` generated in **Single-Cell Data Analysis**, reported clusterwise IC50 and cell death percentages to drugs in the GDSC database via [CaDRReS-Sc](https://github.com/CSB5/CaDRReS-SC) (a recommender system framework for _in silico_ drug response prediction), or drug sensitivity AUC in the PRISM database from [DepMap Portal PRISM-19Q4] (https://doi.org/10.1038/s43018-019-0018-6). The output the prediction results are `IC50_prediction.csv` and `drug_kill_prediction.csv` while using parameter `--model GDSC`, and `AUC_prediction.csv` whlie using parameter `--model PRISM`.

- Run `python3 drug_response_prediction.py -h` to show the help messages as follow for **Drug Response Prediction**.

```
usage: drug_response_prediction.py [-h] -i INPUT [-o OUTPUT] [-c CLUSTERS] [-m MODEL] [--n_drugs N_DRUGS]

Drug response prediction

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input Anndata object (h5ad file)
  -o OUTPUT, --output OUTPUT
                        path to output directory, default='./'
  -c CLUSTERS, --clusters CLUSTERS
                        perform sensitivity prediction on specified clusters, e.g. '1,3,8,9', default='All'
  -m MODEL, --model MODEL
                        the sensitivity screening is from GDSC ic50/PRISM auc, e.g. GDSC, PRISM
  --n_drugs N_DRUGS     the number of drugs to visualize for each cluster
```

- Predict drug response on specified clusters (here for default all clusters) with **Drug Response Prediction**.

```
python3 drug_response_prediction.py --input scanpyobj.sub.h5ad
```

### Drug Combination Prediction

- List drug combinations with different clusters from NMF Leiden results and their synergy scores tested on cell lines with **Drug Combination Prediction**.

```
python3 drug_combination_prediction(NMF_Leiden_Predict).py --input scanpyobj.sub.h5ad --c_pre NMF_Leiden_Clustering_Result.csv (--model = PRISM)

python3 drug_combination_prediction(NMF_Leiden_Predict).py --input scanpyobj.sub.h5ad --c_pre NMF_Leiden_Clustering_Result.csv --gdsc_smiles gdsc_drug_smiles_list.csv --model GDSC (--model = GDSC)

python3 drug_combination_prediction(MoA).py --input scanpyobj.sub.h5ad --depmap_prism depmap_prism_drug_info.csv (--model = PRISM)

python3 drug_combination_prediction(MoA).py --input scanpyobj.sub.h5ad --depmap_prism depmap_prism_drug_info.csv --depmap_prism_smiles depmap_prism_drug_smiles_list.csv --gdsc_smiles gdsc_drug_smiles_list.csv --model GDSC (--model = GDSC)
```
