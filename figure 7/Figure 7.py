import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
config = {'method': 'Bonferroni', 'alpha': 0.05, 'logFC_threshold': 1, 'up_or_down_or_both': 'both_up_and_down_regulated', 'highlight_genes': [], 'top_genes': 20, 'upreg_criteria': None, 'downreg_criteria': None, 'upregulated_color': '#FF0000', 'downregulated_color': '#0000FF', 'not_significant_color': '#808080', 'not_significant_between_logFC_threshold_color': '#A0A0A0', 'pie_chart_size': 0.2, 'pie_chart_position': (0.85, 0.85), 'ranking_method': 'MAS', 'Magnitude (M)': 0, 'Altitude (A)': 1, 'g100': (-1, 0.2), 'g010': (1, 0.24), 'g001': (-0.4, -0.5), 'g110': (-0.6, 0.24), 'g101': (-0.8, -0.3), 'g011': (0.7, -0.1), 'g111': (0.4, -0.6)}
import numpy as np

def MAS(df, M, A, method, alpha, logFC_threshold, up_or_down_or_both):
    if method == 'PValue':
        y_column = 'PValue'
    elif method == 'BH':
        df['BH_adjusted_p-value'] = BH_adjusted_pvalue(df['PValue'].values, alpha)
        y_column = 'BH_adjusted_p-value'
    elif method == 'Bonferroni':
        n_tests = len(df)
        df['Bonferroni_adjusted_p-value'] = np.minimum(df['PValue'] * n_tests, 1.0)
        y_column = 'Bonferroni_adjusted_p-value'
    elif method == 'BY':
        df['BY_adjusted_p-value'] = BY_adjusted_pvalue(df['PValue'].values, alpha)
        y_column = 'BY_adjusted_p-value'
    elif method == 'Holm':
        df['Holm_adjusted_p-value'] = Holm_adjusted_pvalue(df['PValue'].values, alpha)
        y_column = 'Holm_adjusted_p-value'
    elif method == 'Hochberg':
        df['Hochberg_adjusted_p-value'] = Hochberg_adjusted_pvalue(df['PValue'].values, alpha)
        y_column = 'Hochberg_adjusted_p-value'
    else:
        y_column = 'PValue'
    df[y_column] = df[y_column].clip(upper=1)
    adjusted_alpha = alpha / len(df) if method == 'Bonferroni' else alpha
    if up_or_down_or_both == 'Only_upregulated':
        df['MAS_Score'] = np.where((df['logFC'] > logFC_threshold) & (df[y_column] <= adjusted_alpha), np.abs(df['logFC']) ** M * np.abs(-np.log10(df[y_column])) ** A, 0)
    elif up_or_down_or_both == 'Only_downregulated':
        df['MAS_Score'] = np.where((df['logFC'] < -logFC_threshold) & (df[y_column] <= adjusted_alpha), np.abs(df['logFC']) ** M * np.abs(-np.log10(df[y_column])) ** A, 0)
    elif up_or_down_or_both == 'both_up_and_down_regulated':
        df['MAS_Score'] = np.where((df['logFC'].abs() > logFC_threshold) & (df[y_column] <= adjusted_alpha), np.abs(df['logFC']) ** M * np.abs(-np.log10(df[y_column])) ** A, 0)
    df = df.sort_values(by='MAS_Score', ascending=False).reset_index(drop=True)
    df['MAS_rank'] = df.index + 1
    return df

def BH_adjusted_pvalue(p, alpha):
    sort = np.argsort(p)
    rank = np.zeros(len(p))
    j = 1
    for i in range(len(sort)):
        rank[sort[i]] = j
        j = j + 1
    new_alpha = rank / len(p) * alpha
    L = []
    for j in range(len(rank)):
        if p[np.where(rank == len(rank) - j)[0][0]] < new_alpha[np.where(rank == len(rank) - j)[0][0]]:
            L.append(np.where(rank == len(rank) - j)[0][0])
    new_p = p * (len(p) / rank)
    H = np.zeros(len(p))
    H[np.where(rank == len(p))[0][0]] = new_p[np.where(rank == len(p))[0][0]]
    for k in range(1, len(p)):
        if new_p[np.where(rank == len(p) - k)[0][0]] > new_p[np.where(rank == len(p) - k + 1)[0][0]]:
            H[np.where(rank == len(p) - k)[0][0]] = new_p[np.where(rank == len(p) - k + 1)[0][0]]
        else:
            H[np.where(rank == len(p) - k)[0][0]] = new_p[np.where(rank == len(p) - k)[0][0]]
    BH_adjusted_p_values = np.where(H > 1, 1, H)
    return BH_adjusted_p_values
from matplotlib.ticker import MaxNLocator
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from itertools import combinations

def rank_genes(df, gene_set, ranking_method='MAS'):
    if ranking_method not in ['MAS', 'corner_distance']:
        raise ValueError("Ranking method should be 'MAS' or 'corner_distance'")
    score_field = 'MAS_Score' if ranking_method == 'MAS' else 'corner_distance'
    ranked_genes = df[df['Gene Symbol'].isin(gene_set)].copy()
    if ranking_method == 'MAS':
        ranked_genes = ranked_genes.nlargest(10, score_field)
    else:
        ranked_genes = ranked_genes.nsmallest(10, score_field)
    return ranked_genes['Gene Symbol'].tolist()
import pandas as pd
from itertools import combinations
import pandas as pd
from itertools import combinations

data1 = pd.read_csv('GSE152418-original.csv')
data1.set_index('Gene Symbol', inplace=True)

new_columns = []
control_count = 1
sarscov2_count = 1
for col in data1.columns:
    if col.startswith('Healthy'):
        new_columns.append(f'Control ({control_count})')
        control_count += 1
    else:
        new_columns.append(f'SARS-CoV-2 ({sarscov2_count})')
        sarscov2_count += 1
data1.columns = new_columns

control_cols = [col for col in data1.columns if col.startswith('Control')]
sars_cols = [col for col in data1.columns if not col.startswith('Control')]
data1 = data1[control_cols + sars_cols]

data1.to_csv('GSE152418-Clean.csv', index=True)

data2 = pd.read_csv('GSE161731-original.csv')
data2.set_index('Gene Symbol', inplace=True)

group_map = data2.columns.str.extract('^([\\w\\s]+?)\\s*\\(')[0]
group_map.index = data2.columns
grouped_dataframes = {group: data2.loc[:, group_map[group_map == group].index].copy() for group in group_map.unique()}

combined_df = pd.concat([grouped_dataframes['HLTY'], grouped_dataframes['SARS_CoV_2_Yes']], axis=1)
combined_df.index.name = 'Gene Symbol'

new_columns = []
control_count = 1
sarscov2_count = 1
for col in combined_df.columns:
    if col.startswith('HLTY'):
        new_columns.append(f'Control ({control_count})')
        control_count += 1
    else:
        new_columns.append(f'SARS-CoV-2 ({sarscov2_count})')
        sarscov2_count += 1
combined_df.columns = new_columns

combined_df.to_csv('GSE161731-Clean.csv', index=True)

data3 = pd.read_csv('GSE171110-original.csv')
data3.set_index('Gene Symbol', inplace=True)

new_columns = []
control_count = 1
sarscov2_count = 1
for col in data3.columns:
    if col.startswith('HLTY'):
        new_columns.append(f'Control ({control_count})')
        control_count += 1
    else:
        new_columns.append(f'SARS-CoV-2 ({sarscov2_count})')
        sarscov2_count += 1
data3.columns = new_columns

data3.to_csv('GSE171110-Clean.csv', index=True)

data4 = pd.read_csv('PMC8202013-original.csv')
data4.set_index('Gene Symbol', inplace=True)

group_map = data4.columns.str.extract('^([\\w\\s]+?)\\s*\\(')[0]
group_map.index = data4.columns
grouped_dataframes = {group: data4.loc[:, group_map[group_map == group].index].copy() for group in group_map.unique()}

combined_df = pd.concat([grouped_dataframes['HLTY'], grouped_dataframes['SARS_CoV_2']], axis=1)
combined_df.index.name = 'Gene Symbol'

new_columns = []
control_count = 1
sarscov2_count = 1
for col in combined_df.columns:
    if col.startswith('HLTY'):
        new_columns.append(f'Control ({control_count})')
        control_count += 1
    else:
        new_columns.append(f'SARS-CoV-2 ({sarscov2_count})')
        sarscov2_count += 1
combined_df.columns = new_columns

combined_df.to_csv('PMC8202013-Clean.csv', index=True)

import pandas as pd
filenames = ['GSE152418-Clean.csv', 'GSE161731-Clean.csv', 'GSE171110-Clean.csv', 'PMC8202013-Clean.csv']
dataframes = {}
for fname in filenames:
    df = pd.read_csv(fname)
    df = df.drop_duplicates(subset='Gene Symbol')
    df.set_index('Gene Symbol', inplace=True)
    dataframes[fname] = df
common_genes = set.intersection(*(set(df.index) for df in dataframes.values()))
for fname, df in dataframes.items():
    filtered_df = df.loc[sorted(common_genes)]
    new_name = fname.replace('-Clean', '-Ready')
    filtered_df.to_csv(new_name)

import subprocess

def run_r_script():
    command = ['Rscript', 'TMM_GLMQL.R']
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if result.stdout:
        print('R script output:', result.stdout)
    if result.stderr:
        print('R script errors:', result.stderr)
if __name__ == '__main__':
    run_r_script()

import subprocess

def run_r_script():
    command = ['Rscript', 'DESeq2.R']
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if result.stdout:
        print('R script output:', result.stdout)
    if result.stderr:
        print('R script errors:', result.stderr)
if __name__ == '__main__':
    run_r_script()

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, precision_score, recall_score
datasets = ['GSE152418', 'GSE161731', 'GSE171110', 'PMC8202013']

def load_sig_genes(tool, dataset):
    df = pd.read_csv(f'{tool}-{dataset}.csv')
    df = df[(df['Bonferroni_pvalue'] < 0.05) & (abs(df['logFC']) > 1)]
    return set(df['GeneSymbol'])

def evaluate_pca(common_genes, raw_file):
    raw_df = pd.read_csv(raw_file)
    raw_df.set_index('Gene Symbol', inplace=True)
    raw_df = raw_df.loc[raw_df.index.intersection(common_genes)]
    data = raw_df.T
    labels = data.index.str.contains('SARS-CoV-2').astype(int)
    scaler = StandardScaler()
    scaled = scaler.fit_transform(data)
    pcs = PCA(n_components=2).fit_transform(scaled)
    pc1 = pcs[:, 0]
    pred = (pc1 > np.median(pc1)).astype(int)
    auc = roc_auc_score(labels, pc1)
    precision = precision_score(labels, pred)
    recall = recall_score(labels, pred)
    return (auc, precision, recall)
results = []
for tool in ['DESeq2', 'EdgeR']:
    for test in datasets:
        train = [d for d in datasets if d != test]
        common_genes = set.intersection(*[load_sig_genes(tool, d) for d in train])
        auc, prec, rec = evaluate_pca(common_genes, f'{datasets.index(test) + 1}-{test}-Ready.csv')
        results.append({'Tool': tool, 'TestSet': test, 'NumGenes': len(common_genes), 'AUC': auc, 'Precision': prec, 'Recall': rec})
results_df = pd.DataFrame(results)
results_df.to_csv('CrossValidation_EdgeR_DESeq2_PCA.csv', index=False)
print(results_df)

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, precision_score, recall_score
datasets = ['GSE152418', 'GSE161731', 'GSE171110', 'PMC8202013']

def load_sig_genes(tool, dataset):
    df = pd.read_csv(f'{tool}-{dataset}.csv')
    df = df[(df['FDR'] < 0.05) & (abs(df['logFC']) > 1)]
    return set(df['GeneSymbol'])

def evaluate_pca(genes, raw_file):
    raw_df = pd.read_csv(raw_file)
    raw_df.set_index('Gene Symbol', inplace=True)
    raw_df = raw_df.loc[raw_df.index.intersection(genes)]
    if raw_df.shape[0] < 2:
        return (np.nan, np.nan, np.nan, raw_df.shape[0])
    data = raw_df.T
    labels = data.index.str.contains('SARS-CoV-2').astype(int)
    scaler = StandardScaler()
    scaled = scaler.fit_transform(data)
    pcs = PCA(n_components=2).fit_transform(scaled)
    pc1 = pcs[:, 0]
    pred = (pc1 > np.median(pc1)).astype(int)
    auc = roc_auc_score(labels, pc1)
    precision = precision_score(labels, pred)
    recall = recall_score(labels, pred)
    return (auc, precision, recall, raw_df.shape[0])
results = []
for tool in ['DESeq2', 'EdgeR']:
    other = 'EdgeR' if tool == 'DESeq2' else 'DESeq2'
    for test_set in datasets:
        train_sets = [d for d in datasets if d != test_set]
        tool_genes = set.intersection(*[load_sig_genes(tool, d) for d in train_sets])
        other_genes = set.intersection(*[load_sig_genes(other, d) for d in train_sets])
        unique_genes = tool_genes - other_genes
        test_index = datasets.index(test_set) + 1
        raw_file = f'{test_index}-{test_set}-Ready.csv'
        auc, prec, rec, n_genes = evaluate_pca(unique_genes, raw_file)
        results.append({'Tool': tool, 'TestSet': test_set, 'NumGenes': n_genes, 'AUC': auc, 'Precision': prec, 'Recall': rec})
df_results = pd.DataFrame(results)
df_results.to_csv('UniqueCommonGenes_EdgeR_DESeq2.csv', index=False)
print(df_results)

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from scipy import interp
tools = ['DESeq2', 'EdgeR']
fpr_grid = np.linspace(0, 1, 100)
roc_data = {tool: [] for tool in tools}
for tool in tools:
    other = 'EdgeR' if tool == 'DESeq2' else 'DESeq2'
    for test_set in datasets:
        train_sets = [d for d in datasets if d != test_set]
        tool_genes = set.intersection(*[load_sig_genes(tool, d) for d in train_sets])
        other_genes = set.intersection(*[load_sig_genes(other, d) for d in train_sets])
        unique_genes = tool_genes - other_genes
        test_index = datasets.index(test_set) + 1
        raw_file = f'{test_index}-{test_set}-Ready.csv'
        auc, prec, rec, n_genes, y_true, y_score = evaluate_pca(unique_genes, raw_file)
        if y_true is not None:
            fpr, tpr, _ = roc_curve(y_true, y_score)
            tpr_interp = np.interp(fpr_grid, fpr, tpr)
            roc_data[tool].append(tpr_interp)
plt.figure(figsize=(8, 6))
colors = {'DESeq2': 'darkred', 'EdgeR': 'blue'}
for tool in tools:
    tprs = np.array(roc_data[tool])
    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)
    plt.plot(fpr_grid, mean_tpr, label=tool, color=colors[tool], lw=2)
    plt.fill_between(fpr_grid, mean_tpr - std_tpr, mean_tpr + std_tpr, color=colors[tool], alpha=0.2)
plt.plot([0, 1], [0, 1], 'k--', lw=1)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC AUC Curves Across Folds (Tool-Unique Genes)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, precision_score, recall_score, roc_curve
datasets = ['GSE152418', 'GSE161731', 'GSE171110', 'PMC8202013']

def load_sig_genes(tool, dataset):
    df = pd.read_csv(f'{tool}-{dataset}.csv')
    df = df[(df['FDR'] < 0.05) & (abs(df['logFC']) > 1)]
    return set(df['GeneSymbol'])

def evaluate_pca(genes, raw_file):
    raw_df = pd.read_csv(raw_file)
    raw_df.set_index('Gene Symbol', inplace=True)
    raw_df = raw_df.loc[raw_df.index.intersection(genes)]
    if raw_df.shape[0] < 2:
        return (np.nan, np.nan, np.nan, raw_df.shape[0], None, None)
    data = raw_df.T
    labels = data.index.str.contains('SARS-CoV-2').astype(int)
    scaler = StandardScaler()
    scaled = scaler.fit_transform(data)
    pcs = PCA(n_components=2).fit_transform(scaled)
    pc1 = pcs[:, 0]
    pred = (pc1 > np.median(pc1)).astype(int)
    auc = roc_auc_score(labels, pc1)
    precision = precision_score(labels, pred)
    recall = recall_score(labels, pred)
    return (auc, precision, recall, raw_df.shape[0], labels, pc1)
tools = ['DESeq2', 'EdgeR']
fpr_grid = np.linspace(0, 1, 100)
roc_data = {tool: [] for tool in tools}
for tool in tools:
    other = 'EdgeR' if tool == 'DESeq2' else 'DESeq2'
    for test_set in datasets:
        train_sets = [d for d in datasets if d != test_set]
        tool_genes = set.intersection(*[load_sig_genes(tool, d) for d in train_sets])
        other_genes = set.intersection(*[load_sig_genes(other, d) for d in train_sets])
        unique_genes = tool_genes - other_genes
        test_index = datasets.index(test_set) + 1
        raw_file = f'{test_index}-{test_set}-Ready.csv'
        auc, prec, rec, n_genes, y_true, y_score = evaluate_pca(unique_genes, raw_file)
        if y_true is not None:
            fpr, tpr, _ = roc_curve(y_true, y_score)
            tpr_interp = np.interp(fpr_grid, fpr, tpr)
            tpr_interp[0] = 0.0
            tpr_interp[-1] = 1.0
            roc_data[tool].append(tpr_interp)
plt.figure(figsize=(8, 6))
colors = {'DESeq2': 'darkred', 'EdgeR': 'blue'}
for tool in tools:
    tprs = np.array(roc_data[tool])
    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)
    plt.plot(fpr_grid, mean_tpr, label=tool, color=colors[tool], lw=2)
    plt.fill_between(fpr_grid, mean_tpr - std_tpr, mean_tpr + std_tpr, color=colors[tool], alpha=0.2)
plt.plot([0, 1], [0, 1], 'k--', lw=1)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC AUC Curves (Unique Common Genes across Folds)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, precision_score, recall_score, accuracy_score, roc_curve
datasets = ['GSE152418', 'GSE161731', 'GSE171110', 'PMC8202013']

def load_sig_genes(tool, dataset):
    df = pd.read_csv(f'{tool}-{dataset}.csv')
    df = df[(df['FDR'] < 0.05) & (abs(df['logFC']) > 1)]
    return set(df['GeneSymbol'])

def evaluate_pca(genes, raw_file):
    raw_df = pd.read_csv(raw_file)
    raw_df.set_index('Gene Symbol', inplace=True)
    raw_df = raw_df.loc[raw_df.index.intersection(genes)]
    if raw_df.shape[0] < 2:
        return (np.nan, np.nan, np.nan, np.nan, raw_df.shape[0], None, None)
    data = raw_df.T
    labels = data.index.str.contains('SARS-CoV-2').astype(int)
    scaler = StandardScaler()
    scaled = scaler.fit_transform(data)
    pcs = PCA(n_components=2).fit_transform(scaled)
    pc1 = pcs[:, 0]
    pred = (pc1 > np.median(pc1)).astype(int)
    auc = roc_auc_score(labels, pc1)
    precision = precision_score(labels, pred)
    recall = recall_score(labels, pred)
    accuracy = accuracy_score(labels, pred)
    return (auc, precision, recall, accuracy, raw_df.shape[0], labels, pc1)
tools = ['DESeq2', 'EdgeR']
fpr_grid = np.linspace(0, 1, 100)
roc_data = {tool: [] for tool in tools}
metrics_data = {tool: {'auc': [], 'precision': [], 'recall': [], 'accuracy': []} for tool in tools}
for tool in tools:
    other = 'EdgeR' if tool == 'DESeq2' else 'DESeq2'
    for test_set in datasets:
        train_sets = [d for d in datasets if d != test_set]
        tool_genes = set.intersection(*[load_sig_genes(tool, d) for d in train_sets])
        other_genes = set.intersection(*[load_sig_genes(other, d) for d in train_sets])
        unique_genes = tool_genes - other_genes
        test_index = datasets.index(test_set) + 1
        raw_file = f'{test_index}-{test_set}-Ready.csv'
        auc, prec, rec, acc, n_genes, y_true, y_score = evaluate_pca(unique_genes, raw_file)
        if y_true is not None:
            fpr, tpr, _ = roc_curve(y_true, y_score)
            tpr_interp = np.interp(fpr_grid, fpr, tpr)
            tpr_interp[0] = 0.0
            tpr_interp[-1] = 1.0
            roc_data[tool].append(tpr_interp)
            metrics_data[tool]['auc'].append(auc)
            metrics_data[tool]['precision'].append(prec)
            metrics_data[tool]['recall'].append(rec)
            metrics_data[tool]['accuracy'].append(acc)
plt.figure(figsize=(8, 6))
colors = {'DESeq2': 'darkred', 'EdgeR': 'blue'}
for tool in tools:
    tprs = np.array(roc_data[tool])
    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)
    aucs = np.array(metrics_data[tool]['auc'])
    precs = np.array(metrics_data[tool]['precision'])
    recalls = np.array(metrics_data[tool]['recall'])
    accs = np.array(metrics_data[tool]['accuracy'])
    legend_text = f'{tool}  AUC={aucs.mean():.2f}±{aucs.std():.2f},  Acc={accs.mean():.2f}±{accs.std():.2f},  \n             Prec={precs.mean():.2f}±{precs.std():.2f},  Rec={recalls.mean():.2f}±{recalls.std():.2f}'
    plt.plot(fpr_grid, mean_tpr, label=legend_text, color=colors[tool], lw=2)
    upper = np.minimum(mean_tpr + std_tpr, 1.0)
    lower = mean_tpr - std_tpr
    plt.fill_between(fpr_grid, lower, upper, color=colors[tool], alpha=0.2)
plt.plot([0, 1], [0, 1], 'k--', lw=1.2)
plt.xlabel('False Positive Rate', fontsize=16)
plt.ylabel('True Positive Rate', fontsize=16)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.title('ROC AUC with Mean ± Std (Unique Genes per Tool)', fontsize=15)
plt.legend(loc='lower right', fontsize=14)
plt.grid(True)
plt.gca().spines['bottom'].set_linewidth(1.8)
plt.gca().spines['left'].set_linewidth(1.8)
plt.tight_layout()
plt.show()

for tool in tools:
    aucs = np.array(metrics_data[tool]['auc'])
    precs = np.array(metrics_data[tool]['precision'])
    recalls = np.array(metrics_data[tool]['recall'])
    accs = np.array(metrics_data[tool]['accuracy'])
    legend_text = f'{tool}  AUC={aucs.mean():.2f}±{aucs.std():.2f},  Acc={accs.mean():.2f}±{accs.std():.2f},  \n             Prec={precs.mean():.2f}±{precs.std():.2f},  Rec={recalls.mean():.2f}±{recalls.std():.2f}'
    print(legend_text)

def evaluate_pca(genes, raw_file):
    raw_df = pd.read_csv(raw_file)
    raw_df.set_index('Gene Symbol', inplace=True)
    raw_df = raw_df.loc[raw_df.index.intersection(genes)]
    if raw_df.shape[0] < 2:
        return (np.nan, np.nan, np.nan, np.nan, raw_df.shape[0], None, None, None, None)
    data = raw_df.T
    labels = data.index.str.contains('SARS-CoV-2').astype(int)
    scaler = StandardScaler()
    scaled = scaler.fit_transform(data)
    pcs = PCA(n_components=2).fit_transform(scaled)
    pc1 = pcs[:, 0]
    pc2 = pcs[:, 1]
    pred = (pc1 > np.median(pc1)).astype(int)
    auc = roc_auc_score(labels, pc1)
    precision = precision_score(labels, pred)
    recall = recall_score(labels, pred)
    accuracy = accuracy_score(labels, pred)
    return (auc, precision, recall, accuracy, raw_df.shape[0], labels, pc1, pc2, data.index)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, precision_score, recall_score, accuracy_score
from scipy.spatial import ConvexHull
import matplotlib
matplotlib.rcParams['font.family'] = 'Arial'

def draw_convex_hull(points, ax, color, alpha=0.2):
    if len(points) < 3:
        return
    hull = ConvexHull(points)
    vertices = np.append(hull.vertices, hull.vertices[0])
    ax.fill(points[vertices, 0], points[vertices, 1], color=color, alpha=alpha, lw=0)
test_set = 'GSE152418'
test_index = datasets.index(test_set) + 1
raw_file = f'{test_index}-{test_set}-Ready.csv'
fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
for ax, tool in zip(axes, ['DESeq2', 'EdgeR']):
    other = 'EdgeR' if tool == 'DESeq2' else 'DESeq2'
    train_sets = [d for d in datasets if d != test_set]
    tool_genes = set.intersection(*[load_sig_genes(tool, d) for d in train_sets])
    other_genes = set.intersection(*[load_sig_genes(other, d) for d in train_sets])
    unique_genes = tool_genes - other_genes
    auc, prec, rec, acc, n_genes, y_true, pc1, pc2, sample_names = evaluate_pca(unique_genes, raw_file)
    if y_true is None:
        ax.set_title(f'{tool}: Not enough genes', fontsize=14)
        continue
    raw_df = pd.read_csv(raw_file)
    raw_df.set_index('Gene Symbol', inplace=True)
    raw_df = raw_df.loc[raw_df.index.intersection(unique_genes)]
    data = raw_df.T
    labels = data.index.str.contains('SARS-CoV-2').astype(int)
    scaler = StandardScaler()
    scaled = scaler.fit_transform(data)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(scaled)
    pc1, pc2 = (pcs[:, 0], pcs[:, 1])
    evr = pca.explained_variance_ratio_
    for label_value, label_name, color in [(0, 'Control', 'green'), (1, 'SARS-CoV-2', 'red')]:
        idx = labels == label_value
        ax.scatter(pc1[idx], pc2[idx], label=label_name, color=color, s=300, edgecolor='black', alpha=0.85)
        draw_convex_hull(np.vstack([pc1[idx], pc2[idx]]).T, ax, color=color, alpha=0.15)
    box_text = f'AUC: {auc:.3f}, Accuracy: {acc:.3f}, \nPrecision: {prec:.3f}, Recall: {rec:.3f}.'
    ax.text(0.02, 0.02, box_text, transform=ax.transAxes, fontsize=15, verticalalignment='bottom', horizontalalignment='left', bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='black', alpha=0.9))
    ax.set_title(f'{tool} (GSE152418)', fontsize=16)
    ax.set_xlabel(f'PC1 ({evr[0] * 100:.1f}%)', fontsize=14)
    if ax is axes[0]:
        ax.set_ylabel(f'PC2 ({evr[1] * 100:.1f}%)', fontsize=14)
    ax.tick_params(labelsize=12)
    ax.grid(True)
    ax.set_facecolor('#f8f8f8')
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=2, fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.92])
plt.show()

gene_list = sorted(unique_genes)
gene_df = pd.DataFrame(gene_list, columns=['UniqueGenes'])
gene_df.to_csv(f'{tool}_UniqueGenes_GSE152418.csv', index=False)
print(f'\n{tool} unique genes for GSE152418 ({len(gene_list)} genes):')
print(gene_list)

import pandas as pd
datasets = ['GSE152418', 'GSE161731', 'GSE171110', 'PMC8202013']

def load_sig_genes(tool, dataset):
    df = pd.read_csv(f'{tool}-{dataset}.csv')
    df = df[(df['FDR'] < 0.05) & (abs(df['logFC']) > 1)]
    return set(df['GeneSymbol'])

def compute_overlap_metrics(A, B):
    intersection = A & B
    union = A | B
    do_A_B = len(intersection) / len(B) if len(B) > 0 else float('nan')
    do_B_A = len(intersection) / len(A) if len(A) > 0 else float('nan')
    jaccard = len(intersection) / len(union) if len(union) > 0 else float('nan')
    return {'DO(DESeq2, edgeR)': round(do_A_B, 4), 'DO(edgeR, DESeq2)': round(do_B_A, 4), 'Jaccard': round(jaccard, 4), 'DESeq2_n': len(A), 'edgeR_n': len(B), 'Common': len(intersection), 'DESeq2_Unique': len(A - B), 'edgeR_Unique': len(B - A)}
results = []
for test_set in datasets:
    train_sets = [d for d in datasets if d != test_set]
    deseq2_genes = set.intersection(*[load_sig_genes('DESeq2', d) for d in train_sets])
    edger_genes = set.intersection(*[load_sig_genes('EdgeR', d) for d in train_sets])
    metrics = compute_overlap_metrics(deseq2_genes, edger_genes)
    metrics['TestSet'] = test_set
    results.append(metrics)
df_overlap = pd.DataFrame(results)
df_overlap = df_overlap[['TestSet', 'DO(DESeq2, edgeR)', 'DO(edgeR, DESeq2)', 'Jaccard', 'DESeq2_n', 'edgeR_n', 'Common', 'DESeq2_Unique', 'edgeR_Unique']]
df_overlap.to_csv('DirectionalOverlap_Jaccard_UniqueGenes.csv', index=False)
print(df_overlap)

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
df = pd.read_csv('DirectionalOverlap_Jaccard_UniqueGenes.csv')
melted = df.melt(id_vars='TestSet', value_vars=['DO(DESeq2, edgeR)', 'DO(edgeR, DESeq2)', 'Jaccard'], var_name='Metric', value_name='Score')
plt.figure(figsize=(10, 6))
sns.barplot(data=melted, x='TestSet', y='Score', hue='Metric', palette=['darkred', 'blue', 'gray'])
plt.title('Directional Overlap and Jaccard Index Across Folds', fontsize=16)
plt.ylabel('Score', fontsize=14)
plt.xlabel('Test Set', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(title='', fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.show()

import matplotlib.pyplot as plt
import seaborn as sns
df_count = df[['TestSet', 'DESeq2_Unique', 'edgeR_Unique', 'Common']].copy()
df_count['TestSet'] = df_count['TestSet'].astype(str)
df_count.set_index('TestSet', inplace=True)
plt.figure(figsize=(10, 6))
df_count.plot(kind='bar', stacked=True, color=['green', 'steelblue', 'darkgray'], edgecolor='black', width=0.7, linewidth=0.6)
plt.title('Common and Unique Genes Across Folds', fontsize=18, weight='bold')
plt.ylabel('Number of Genes', fontsize=15)
plt.xlabel('Test Set', fontsize=15)
plt.xticks(rotation=0, fontsize=13)
plt.yticks(fontsize=13)
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.legend(['DESeq2 Unique', 'edgeR Unique', 'Common'], fontsize=13, title='Gene Type', title_fontsize=14)
plt.tight_layout()
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
datasets = ['GSE152418', 'GSE161731', 'GSE171110', 'PMC8202013']

def load_sig_genes(tool, dataset):
    df = pd.read_csv(f'{tool}-{dataset}.csv')
    df = df[(df['FDR'] < 0.05) & (abs(df['logFC']) > 1)]
    return set(df['GeneSymbol'])
training_labels = []
common_counts = []
deseq2_unique_counts = []
edger_unique_counts = []
for i in range(len(datasets)):
    test_set = datasets[i]
    train_set = [d for j, d in enumerate(datasets) if j != i]
    label = '\n'.join(train_set)
    deseq2_sets = [load_sig_genes('DESeq2', d) for d in train_set]
    edger_sets = [load_sig_genes('EdgeR', d) for d in train_set]
    deseq2_common = set.intersection(*deseq2_sets)
    edger_common = set.intersection(*edger_sets)
    common = deseq2_common & edger_common
    deseq2_only = deseq2_common - edger_common
    edger_only = edger_common - deseq2_common
    training_labels.append(label)
    common_counts.append(len(common))
    deseq2_unique_counts.append(len(deseq2_only))
    edger_unique_counts.append(len(edger_only))
df = pd.DataFrame({'Training Sets': training_labels, 'Common': common_counts, 'DESeq2 Unique': deseq2_unique_counts, 'edgeR Unique': edger_unique_counts})
df_plot = df.set_index('Training Sets')[['DESeq2 Unique', 'edgeR Unique', 'Common']]
df_plot.plot(kind='bar', stacked=True, figsize=(8, 6), color=['darkred', 'steelblue', 'darkgray'], edgecolor='black')
plt.title('Common and Unique Genes in Training Combinations', fontsize=16)
plt.ylabel('Number of significant Genes', fontsize=16)
plt.xlabel('Training Sets', fontsize=18)
plt.xticks(rotation=0, fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16, title_fontsize=18, ncol=3)
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
datasets = ['GSE152418', 'GSE161731', 'GSE171110', 'PMC8202013']

def load_sig_genes(tool, dataset):
    df = pd.read_csv(f'{tool}-{dataset}.csv')
    df = df[(df['FDR'] < 0.05) & (abs(df['logFC']) > 1)]
    return set(df['GeneSymbol'])
training_labels = []
common_counts = []
deseq2_unique_counts = []
edger_unique_counts = []
for i in range(len(datasets)):
    test_set = datasets[i]
    train_set = [d for j, d in enumerate(datasets) if j != i]
    label = '\n'.join(train_set)
    deseq2_sets = [load_sig_genes('DESeq2', d) for d in train_set]
    edger_sets = [load_sig_genes('EdgeR', d) for d in train_set]
    deseq2_common = set.intersection(*deseq2_sets)
    edger_common = set.intersection(*edger_sets)
    common = deseq2_common & edger_common
    deseq2_only = deseq2_common - edger_common
    edger_only = edger_common - deseq2_common
    training_labels.append(label)
    common_counts.append(len(common))
    deseq2_unique_counts.append(len(deseq2_only))
    edger_unique_counts.append(len(edger_only))
df = pd.DataFrame({'Training Sets': training_labels, 'Common': common_counts, 'DESeq2 Unique': deseq2_unique_counts, 'edgeR Unique': edger_unique_counts})
df_plot = df.set_index('Training Sets')[['DESeq2 Unique', 'edgeR Unique', 'Common']]
df_plot.plot(kind='bar', stacked=True, figsize=(8, 6), color=['darkred', 'steelblue', 'darkgray'], edgecolor='black')
plt.title('Common and Unique Genes in Training Combinations', fontsize=16)
plt.ylabel('Number of significant Genes', fontsize=16)
plt.xlabel('Training Sets', fontsize=18)
plt.xticks(rotation=0, fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16, title_fontsize=18, ncol=3)
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

deseq2_sets = [load_sig_genes('DESeq2', d) for d in datasets]
edger_sets = [load_sig_genes('EdgeR', d) for d in datasets]
common_deseq2 = set.intersection(*deseq2_sets)
common_edger = set.intersection(*edger_sets)
print(f'DESeq2 Common Genes Across All 4: {len(common_deseq2)}')
print(f'edgeR Common Genes Across All 4: {len(common_edger)}')

import matplotlib.pyplot as plt
tools = ['DESeq2', 'edgeR']
common_counts = [len(common_deseq2), len(common_edger)]
plt.figure(figsize=(6, 5))
bars = plt.bar(tools, common_counts, color=['darkred', 'steelblue'], edgecolor='black')
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2, height + 2, f'{height}', ha='center', fontsize=13)
plt.ylabel('Number of Common DE Genes', fontsize=14)
plt.title('Common Significant Genes Across All 4 Datasets', fontsize=15, weight='bold')
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.show()
