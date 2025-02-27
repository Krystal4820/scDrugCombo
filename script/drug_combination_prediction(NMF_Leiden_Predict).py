import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import itertools
from scipy.stats import gmean
import argparse
import sys
import os
from ipywidgets import widgets
import warnings
warnings.filterwarnings('ignore')

# Parse command-line arguments
# process arguments
parser = argparse.ArgumentParser(description='Drug combination prediction')

parser.add_argument('-i', '--input', required=True,
                    help='path to input Anndata object (h5ad file)')
parser.add_argument('-o', '--output', default='./',
                    help='path to output directory, default=\'./\'')
parser.add_argument(
    '--c_pre', help='result of NMF Leiden Clustering')
parser.add_argument(
    '--gdsc_smiles', help='GDSC Drug Smiles List')
parser.add_argument('-c', '--clusters', default='All', type=str,
                    help='perform IC50 prediction on specified clusters, e.g. \'1,3,8,9\', default=\'All\'')
parser.add_argument('-m', '--model', default='PRISM', type=str,
                    help='the sensitivity screening is from GDSC ic50/PRISM auc, e.g. GDSC, PRISM')
parser.add_argument('--n_drugs', default=10, type=int,
                    help='the number of drugs to visualize for each cluster')

args = parser.parse_args()

# check input
if not os.path.exists(args.input):
    sys.exit('The input path does not exist.')
if args.input[-5:] != '.h5ad':
    sys.exit('The input file is not a h5ad file.')

# check output
if not os.path.isdir(args.output):
    sys.exit('The output directory does not exist.')

if args.c_pre is None:
    sys.exit("Error: --c_pre argument is missing. Please provide a valid file path.")
if not os.path.exists(args.c_pre):
    sys.exit(f"Error: The file {args.c_pre} does not exist.")

# check model
if args.model == 'GDSC':
    if not args.input or not args.c_pre or not args.gdsc_smiles:
        sys.exit(
            "Error: GDSC model requires --input, --c_pre, and --gdsc_smiles parameters.")
else:
    if not args.input or not args.c_pre:
        sys.exit(
            "Error: PRISM model requires --input, and --c_pre parameters.")

if args.model.upper() == "PRISM" and args.gdsc_smiles:
    parser.error(
        "--gdsc_smiles should not be provided when --model is set to PRISM")

scriptpath = '../../opt/CaDRReS-Sc/'
sys.path.append(os.path.abspath(scriptpath))


class Drug_Combination:
    def __init__(self):
        self.load_model()
        self.drug_info()
        self.bulk_exp()
        self.sc_exp()
        self.kernel_feature_preparartion()
        self.sensitivity_prediction()
        if args.model == 'GDSC':
            self.masked_drugs = list(pd.read_csv(
                '../../data/masked_drugs.csv')['GDSC'].dropna().astype('int64').astype('str'))
            self.cell_death_proportion()
        else:
            self.masked_drugs = list(pd.read_csv(
                '../../data/masked_drugs.csv')['PRISM'])
        self.output_result()
        self.load_drug()                                   # --input
        if args.model == 'GDSC':
            self.select_366_drug_smiles_gdsc()             # --gdsc_smiles / --c_pre
            self.merge_pre_labels_gdsc()
        else:
            self.select_369_drug_name_prism()              # --c_pre
            self.merge_pre_labels_prism()
        self.compute_synergy_score()                       # --output

    def load_model(self):
        # IC50/AUC prediction
        # Read pre-trained model
        from cadrres_sc import pp, model, evaluation, utility
        model_dir = '../../CaDRReS-Sc-model/'
        obj_function = widgets.Dropdown(options=[
                                        'cadrres-wo-sample-bias', 'cadrres-wo-sample-bias-weight'], description='Objetice function')
        self.model_spec_name = obj_function.value
        if args.model == 'GDSC':
            model_file = model_dir + \
                '{}_param_dict_all_genes.pickle'.format(self.model_spec_name)
        elif args.model == 'PRISM':
            model_file = model_dir + \
                '{}_param_dict_prism.pickle'.format(self.model_spec_name)
        else:
            sys.exit('Wrong model name.')
        self.cadrres_model = model.load_model(model_file)

    def drug_info(self):
        # Read drug information
        if args.model == 'GDSC':
            self.drug_info_df = pd.read_csv(
                scriptpath + '/preprocessed_data/GDSC/drug_stat.csv', index_col=0)
            self.drug_info_df.index = self.drug_info_df.index.astype(str)
        else:
            self.drug_info_df = pd.read_csv(
                scriptpath + '/preprocessed_data/PRISM/PRISM_drug_info.csv', index_col='broad_id')

    def bulk_exp(self):
        # Read test data
        if args.model == 'GDSC':
            self.gene_exp_df = pd.read_csv(
                scriptpath + '/data/GDSC/GDSC_exp.tsv', sep='\t', index_col=0)
            self.gene_exp_df = self.gene_exp_df.groupby(
                self.gene_exp_df.index).mean()
        else:
            self.gene_exp_df = pd.read_csv(
                scriptpath + '/data/CCLE/CCLE_expression.csv', low_memory=False, index_col=0).T
            self.gene_exp_df.index = [gene.split(
                sep=' (')[0] for gene in self.gene_exp_df.index]

    def sc_exp(self):
        # Load cluster-specific gene expression profile
        self.adata = sc.read(args.input)
        if args.clusters == 'All':
            clusters = sorted(self.adata.obs['louvain'].unique(), key=int)
        else:
            clusters = [x.strip() for x in args.clusters.split(',')]

        self.cluster_norm_exp_df = pd.DataFrame(
            columns=clusters, index=self.adata.raw.var.index)
        for cluster in clusters:
            self.cluster_norm_exp_df[cluster] = self.adata.raw.X[self.adata.obs['louvain'] == cluster].mean(axis=0).T \
                if np.sum(self.adata.raw.X[self.adata.obs['louvain'] == cluster]) else 0.0

    def kernel_feature_preparartion(self):
        from cadrres_sc import pp, model, evaluation, utility
        # Read essential genes list
        if args.model == 'GDSC':
            ess_gene_list = self.gene_exp_df.index.dropna().tolist()
        else:
            ess_gene_list = utility.get_gene_list(
                scriptpath + '/preprocessed_data/PRISM/feature_genes.txt')

        # Calculate fold-change
        cell_line_log2_mean_fc_exp_df, cell_line_mean_exp_df = pp.gexp.normalize_log2_mean_fc(
            self.gene_exp_df)

        self.adata_exp_mean = pd.Series(self.adata.raw.X.mean(
            axis=0).tolist()[0], index=self.adata.raw.var.index)
        cluster_norm_exp_df = self.cluster_norm_exp_df.sub(
            self.adata_exp_mean, axis=0)

        # Calculate kernel feature
        self.test_kernel_df = pp.gexp.calculate_kernel_feature(
            cluster_norm_exp_df, cell_line_log2_mean_fc_exp_df, ess_gene_list)

    def sensitivity_prediction(self):
        from cadrres_sc import pp, model, evaluation, utility
        # Drug response prediction
        if args.model == 'GDSC':
            print('Predicting drug response for using CaDRReS(GDSC): {}'.format(
                self.model_spec_name))
            self.pred_ic50_df, P_test_df = model.predict_from_model(
                self.cadrres_model, self.test_kernel_df, self.model_spec_name)
            print('done!')
        else:
            print('Predicting drug response for using CaDRReS(PRISM): {}'.format(
                self.model_spec_name))
            self.pred_auc_df, P_test_df = model.predict_from_model(
                self.cadrres_model, self.test_kernel_df, self.model_spec_name)
            print('done!')

    def cell_death_proportion(self):
        # Drug kill prediction
        ref_type = 'log2_median_ic50'
        self.drug_list = [
            x for x in self.pred_ic50_df.columns if not x in self.masked_drugs]
        self.drug_info_df = self.drug_info_df.loc[self.drug_list]
        self.pred_ic50_df = self.pred_ic50_df.loc[:, self.drug_list]

        # Predict cell death percentage at the ref_type dosage
        pred_delta_df = pd.DataFrame(
            self.pred_ic50_df.values - self.drug_info_df[ref_type].values, columns=self.pred_ic50_df.columns)
        pred_cv_df = 100 / (1 + (np.power(2, -pred_delta_df)))
        self.pred_kill_df = 100 - pred_cv_df

    def output_result(self):
        if args.model == 'GDSC':
            drug_df = pd.DataFrame({'Drug ID': self.drug_list,
                                    'Drug Name': [self.drug_info_df.loc[drug_id]['Drug Name'] for drug_id in self.drug_list]})
            self.pred_ic50_df = (self.pred_ic50_df.T-self.pred_ic50_df.min(axis=1))/(
                self.pred_ic50_df.max(axis=1)-self.pred_ic50_df.min(axis=1))
            self.pred_ic50_df = self.pred_ic50_df.T
            self.pred_ic50_df.columns = pd.MultiIndex.from_frame(drug_df)
            self.pred_ic50_df.round(3).to_csv(
                os.path.join(args.output, 'IC50_prediction.csv'))
            self.pred_kill_df.columns = pd.MultiIndex.from_frame(drug_df)
            self.pred_kill_df.round(3).to_csv(os.path.join(
                args.output, 'drug_kill_prediction.csv'))
        else:
            drug_list = list(self.pred_auc_df.columns)
            drug_list = [d for d in drug_list if d not in self.masked_drugs]
            drug_df = pd.DataFrame({'Drug ID': drug_list,
                                    'Drug Name': [self.drug_info_df.loc[d, 'name'] for d in drug_list]})
            self.pred_auc_df = self.pred_auc_df.loc[:, drug_list].T
            self.pred_auc_df = (self.pred_auc_df-self.pred_auc_df.min()) / \
                (self.pred_auc_df.max()-self.pred_auc_df.min())
            self.pred_auc_df = self.pred_auc_df.T
            self.pred_auc_df.columns = pd.MultiIndex.from_frame(drug_df)
            self.pred_auc_df.round(3).to_csv(
                os.path.join(args.output, 'PRISM_prediction.csv'))

    def load_drug(self):
        if args.model == 'GDSC':
            self.df = pd.read_csv(os.path.join(
                args.output, 'IC50_prediction.csv'))
        else:
            self.df = pd.read_csv(os.path.join(
                args.output, 'PRISM_prediction.csv'))
        self.df_T = self.df.T
        self.df_T.columns = self.df_T.iloc[0]
        self.df_T = self.df_T[1:].reset_index(drop=False)
        self.df_T.rename(columns={'index': 'Drug ID'}, inplace=True)

        print(self.df_T)

    def select_366_drug_smiles_gdsc(self):
        cleaned_gdsc_alldrug_smiles_df = pd.read_csv(args.gdsc_smiles)
        drug_369_df = pd.read_csv(args.c_pre)
        self.cleaned_drug_369_df = drug_369_df[drug_369_df['SMILES']
                                               != 'Error: 404']
        drug_smiles_set = set(self.cleaned_drug_369_df['SMILES'].to_list())
        self.cleaned_gdsc_in369_all_drug_df = cleaned_gdsc_alldrug_smiles_df[cleaned_gdsc_alldrug_smiles_df['SMILES'].isin(
            drug_smiles_set)]

        print(self.cleaned_gdsc_in369_all_drug_df.shape)  # 34,2 (name,smiles)
        print(self.cleaned_gdsc_in369_all_drug_df)

    def select_369_drug_name_prism(self):
        self.one_name_one_broad_id_df = pd.read_csv(args.c_pre)
        drug_name_set = set(self.one_name_one_broad_id_df['name'].to_list())
        self.all_drug_df = self.df_T[self.df_T['Drug Name'].isin(
            drug_name_set)]

        print(self.all_drug_df)

    def merge_pre_labels_gdsc(self):
        predict_label_df = self.cleaned_gdsc_in369_all_drug_df.merge(
            self.cleaned_drug_369_df, left_on='SMILES', right_on='SMILES', how='left')

        predict_label_merge_df = predict_label_df.merge(
            self.df_T, left_on='Compound Name', right_on='Drug Name', how='left')

        df_1 = self.df[1:]
        self.cell_num_list = df_1['Drug ID'].to_list()

        self.select_df_1 = predict_label_merge_df[[
            'Drug ID', 'Drug Name', 'predicted_labels']]
        self.select_df_2 = predict_label_merge_df[self.cell_num_list]

        self.select_allcell_df = pd.concat(
            [self.select_df_1, self.select_df_2], axis=1)

        print(self.select_allcell_df.shape)
        print(self.select_allcell_df)

    def merge_pre_labels_prism(self):
        self.predict_label_df = self.all_drug_df.merge(
            self.one_name_one_broad_id_df, left_on='Drug Name', right_on='name', how='left')

        df_1 = self.df[1:]
        self.cell_num_list = df_1['Drug ID'].to_list()

        self.select_df_1 = self.predict_label_df[[
            'Drug ID', 'Drug Name', 'predicted_labels']]
        self.select_df_2 = self.predict_label_df[self.cell_num_list]

        self.select_allcell_df = pd.concat(
            [self.select_df_1, self.select_df_2], axis=1)

        print(self.select_allcell_df.shape)
        print(self.select_allcell_df)

    def compute_synergy_score(self):

        drug_combinations = list(itertools.combinations(
            self.select_allcell_df["Drug Name"], 2))
        synergy_scores = []

        for drug1, drug2 in drug_combinations:
            avg_scores = [drug1, drug2]
            for i in range(len(self.cell_num_list)):
                score1 = self.select_allcell_df.loc[self.select_allcell_df["Drug Name"] == drug1, str(
                    i)].values[0]
                score2 = self.select_allcell_df.loc[self.select_allcell_df["Drug Name"] == drug2, str(
                    i)].values[0]
                avg_score = gmean([float(score1), float(score2)])
                avg_scores.append(avg_score)

            synergy_scores.append(avg_scores)

        columns = ["Drug Name 1", "Drug Name 2"] + \
            [f"Synergy_Score_cell_{i}" for i in range(len(self.cell_num_list))]
        synergy_df = pd.DataFrame(synergy_scores, columns=columns)

        predict_label_merge_df_drug1 = synergy_df.merge(
            self.select_allcell_df, left_on='Drug Name 1', right_on='Drug Name', how='left')
        predict_label_merge_df_drug2 = predict_label_merge_df_drug1.merge(
            self.select_allcell_df, left_on='Drug Name 2', right_on='Drug Name', how='left')
        predict_label_merge_df_drug2.rename(columns={
                                            'predicted_labels_x': 'predicted_labels_1', 'predicted_labels_y': 'predicted_labels_2'}, inplace=True)

        diff_predicted_labels_synergy_score_allcell_df_noselect = predict_label_merge_df_drug2[
            predict_label_merge_df_drug2['predicted_labels_1'] != predict_label_merge_df_drug2['predicted_labels_2']
        ]

        columns_to_keep = ["Drug Name 1", "predicted_labels_1", "Drug Name 2",
                           "predicted_labels_2"] + [f"Synergy_Score_cell_{i}" for i in range(len(self.cell_num_list))]
        diff_predicted_labels_synergy_score_allcell_df = diff_predicted_labels_synergy_score_allcell_df_noselect[
            columns_to_keep]

        print(diff_predicted_labels_synergy_score_allcell_df)

        diff_predicted_labels_synergy_score_allcell_df.to_csv(os.path.join(
            args.output, 'diff_predicted_labels_synergy_score_allcell_df.csv'), index=False)

        print("output file done！")


job = Drug_Combination()
