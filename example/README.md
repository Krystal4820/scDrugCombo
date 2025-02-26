## Example Usage for The scDrugCombo

### Preprocessing

The example data is composed of a random 10% subdata from [GSE156625: Onco-fetal reprogramming of endothelial cells drives immunosuppressive macrophages in Hepatocellular Carcinoma (scRNA-seq)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156625).

Before going through the steps in the scDrug, download and uncompress the example data `data.zip ` from [figshare](https://figshare.com/articles/dataset/example_data_zip/20003180) and put the files under the `example` folder.

### Single-Cell Data Analysis

- First, we execute **Single-Cell Data Analysis** on the 10x-Genomics-formatted mtx directory `data/10x_mtx`, with batch correction of `PatientID` in the metadata `data/metadata.csv`, and clustering at resolution 0.6. Additionally, we assign arguments `--annotation` and `--gsea` to perform cell type annotation and Gene Set Enrichment Analysis (GSEA).
- Required memory: 2.5GB

Note: This step may take a few minutes.

```
mkdir write/clustering

python3 ../script/single_cell_analysis.py --input data/10x_mtx --output write/clustering \
--metadata data/metadata.csv --batch PatientID --resolution 0.6 \
--annotation --gsea
```

- Inspecting the preceding output stored in `write/clustering/scanpyobj.h5ad`, we regard the clusters with tumor cell percentages over twice the normal cell percentages, which consist of clusters 1, 5 and 9, as the tumor clusters. Then, we apply **Single-Cell Data Analysis** once again to carry out sub-clustering on the tumor clusters at automatically determined resolution.

Note: To accelerate the process of automatically determined resolution, increase the number of CPU with arguments `--cpus CPUS`.

```
mkdir write/subclustering

python3 script/single_cell_analysis.py --input write/clustering/scanpyobj.h5ad --output write/subclustering \
--format h5ad --clusters '1,5,9' --resolution 0.6
```

### Drug Response Prediction

- Based on the sub-clustering result `write/subclustering/scanpyobj.h5ad`, we run **Drug Response Prediction** to predict clusterwise IC50 and cell death percentages to drugs in the GDSC database.
- Required memory: 2GB

```
mkdir write/drug_response_prediction

python3 script/drug_response_prediction.py --input write/subclustering/scanpyobj.sub.h5ad --output write/drug_response_prediction
```

### Drug Combination Prediction

```
mkdir write/drug_combination_prediction

python3 script/drug_combination_prediction(NMF_Leiden_Predict).py --input write/subclustering/scanpyobj.sub.h5ad --output write/drug_combination_prediction --c_pre NMF_Leiden_Clustering_Result/NMF_Leiden_Clustering_Result.csv (--model = PRISM)

python3 script/drug_combination_prediction(NMF_Leiden_Predict).py --input write/subclustering/scanpyobj.sub.h5ad --output write/drug_combination_prediction --c_pre NMF_Leiden_Clustering_Result/NMF_Leiden_Clustering_Result.csv --gdsc_smiles GDSC_Drug_Smiles_List/gdsc_drug_smiles_list.csv --model GDSC (--model = GDSC)

python3 script/drug_combination_prediction(MoA).py --input write/subclustering/scanpyobj.sub.h5ad --output write/drug_combination_prediction --depmap_prism DepMap_PRISM_Drug/depmap_prism_drug_info.csv (--model = PRISM)

python3 script/drug_combination_prediction(MoA).py --input write/subclustering/scanpyobj.sub.h5ad --output write/drug_combination_prediction --depmap_prism DepMap_PRISM_Drug/depmap_prism_drug_info.csv --depmap_prism_smiles DepMap_PRISM_Drug/depmap_prism_drug_smiles_list.csv --gdsc_smiles GDSC_Drug_Smiles_List/gdsc_drug_smiles_list.csv --model GDSC (--model = GDSC)
```
