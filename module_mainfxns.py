#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib import cm, colors
import scipy.cluster.hierarchy as sch
import seaborn as sns
import datetime
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import KMeans

# Here are the lists of genes for our main pathways of interest.

wnts = ['Wnt1', 'Wnt2', 'Wnt2b', 'Wnt3', 'Wnt3a', 'Wnt4', 'Wnt5a',
       'Wnt5b', 'Wnt6', 'Wnt7a', 'Wnt7b', 'Wnt8a', 'Wnt8b', 'Wnt9a', 'Wnt9b',
       'Wnt10a', 'Wnt10b', 'Wnt11', 'Wnt16']

wntr = ['Lrp5', 'Lrp6', 'Dvl1', 'Dvl2', 'Dvl3', 'Fzd1', 'Fzd2', 'Fzd3',
       'Fzd4', 'Fzd5', 'Fzd6', 'Fzd7', 'Fzd8', 'Fzd9', 'Fzd10']

bmps = ['Bmp3', 'Bmp4', 'Bmp5', 'Bmp6', 'Bmp7', 'Bmp8a', 'Bmp2',
            'Bmp10', 'Bmp11', 'Bmp15', 'Gdf6', 'Gdf7', 'Gdf5', 'Gdf10',
            'Gdf11']

bmpr = ["Bmpr1a" ,"Bmpr1b" ,"Acvr1"  ,"Acvrl1" ,"Acvr1b" ,"Tgfbr1" , "Acvr1c",
            "Acvr2a", "Acvr2b", "Bmpr2" ,"Tgfbr2"]

notch = ["Dll1", "Dll3","Dll4", "Jag1", "Jag2", "Notch1", "Notch2", 
             "Notch3", "Notch4", "Mfng", "Rfng", "Lfng"]


def get_genes(adata, genes):
    '''This function gets genes of interest that have not been filtered out.
    Input: 
    
    adata: AnnData object whose genes have been filtered
    genes: list of gene names
    
    Output:
    
    list of genes of interest that are actually in the filtered dataset
    '''
    list_genes = []
    for i in genes:
        if True in adata.var_names.str.endswith(i):
            list_genes.append(i)
    return list_genes

def vis_pre_processing(adata, genes_range = (0,10000), counts_range = (0, 400000), title=""):
    '''A histogram of genes/cell and counts/cell, a boxplot of 15 highest
    expressed genes, a scatterplot of genes against counts, finally violin
    plots of genes and total counts. This is to visualize the data before
    we filter genes and cells.'''
    fig, ax = plt.subplots(2, 2, figsize = (9, 9))
    ax[0,0].hist(adata.obs['n_genes_per_cell'][:], bins = 100, range = genes_range)
    ax[0,0].axvline(x=2000, color='r', linestyle='dashed', linewidth=2)
    ax[0,0].grid(False)
    ax[0,0].set_title('Histogram of Number of Genes per Cell')
    ax[0,0].set_xlabel('Number of Genes')
    ax[0,0].set_ylabel('Frequency (# of Cells)')
    
    ax[0,1].hist(adata.obs['n_total_counts_per_cell'][:], bins = 100, range=counts_range)
    ax[0,1].axvline(x=20000, color='r', linestyle='dashed', linewidth=2)
    ax[0,1].grid(False)
    ax[0,1].set_title('Histogram of Counts per Cell')
    ax[0,1].set_xlabel('Log Counts per Cell')
    ax[0,1].set_ylabel('Frequency')
    ax[0,1].set_xscale('log')
    
    sc.pl.highest_expr_genes(adata, n_top=15, ax=ax[1,0], show=False)
    ax[1,0].set_title('15 Highest Expressed Genes')
    sc.pl.scatter(adata, x='n_total_counts_per_cell', y='n_genes', size = 20, 
                  title = 'Counts vs. Genes', ax=ax[1,1], show=False)
    ax[1,1].set_xlabel('Counts')
    ax[1,1].set_ylabel('Genes')
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    return fig

def filter_data(adata, min_counts=2000, min_genes=2000, min_cells=3):
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata

def vis_post_processing(adata, genes_range = (0,10000), counts_range = (0, 400000),title=""):
    '''Histograms of genes and total counts, and finally a scatter plot
    of genes against counts.
    Input: 
    Can specify the range of the histograms for genes and counts per cell.
    '''
    fig,ax = plt.subplots(2, 2, figsize = (8,8))
    ax[0,0].hist(adata.obs['n_genes_per_cell'][:], bins = 100, range = genes_range)
    ax[0,0].axvline(x=2000, color='r', linestyle='dashed', linewidth=2)
    ax[0,0].grid(False)
    ax[0,0].set_title('Histogram of Number of Genes per Cell')
    ax[0,0].set_xlabel('Number of Genes')
    ax[0,0].set_ylabel('Frequency (# of Cells)')
    ax[0,1].hist(adata.obs['n_total_counts_per_cell'][:], bins = 100, range = counts_range)
    ax[0,1].axvline(x=20000, color='r', linestyle='dashed', linewidth=2)
    ax[0,1].grid(False)
    ax[0,1].set_title('Histogram of Counts per Cell')
    ax[0,1].set_xlabel('Log Counts per Cell')
    ax[0,1].set_ylabel('Frequency (# of Cells)')
    ax[0,1].set_xscale('log')
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    sc.pl.scatter(adata, x='n_counts', y='n_genes', size=20,ax=ax[1,0], title='Counts vs. Genes',show=False)
    ax[1,0].set_xlabel('Counts')
    ax[1,0].set_ylabel('Genes')
    ax[1,1].axis('off')
    plt.show()
    return fig

def normalize_data(adata,count,**kwargs):
    '''This function normalizes the data, does a log(x+1) transformation, and sets a raw attribute
    of the anndata object. It also sets the highly-variable genes attribute of the anndata observation 
    parameters.
    
    Inputs:
    adata: AnnData object
    count: value to normalize with
    **kwargs: any other arguments to normalize the total with (applied to sc.pp.normalize_total fxn)
    '''
    sc.pp.normalize_total(adata, target_sum=count, **kwargs)
    sc.pp.log1p(adata)
    adata.raw=adata
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    return adata

def merge_genes(adata, all_genes):
    '''This function "densifies" the anndata object. It preserves highly variable genes and all genes of
    interest.
    Input: 
    
    adata: the AnnData object
    all_genes: list of lists of genes of interest
    
    Output:
    
    AnnData object
    '''
    joint_genes=adata.var.highly_variable
    for i in all_genes:
        for j in i:
            joint_genes = joint_genes|adata.var_names.str.startswith(j)
    adata=adata[:,joint_genes]
    return adata

def scale_data(adata):
    '''This function regresses out the AnnData object againist total counts per cell, and scales the 
    gene expression matrix so that each gene has zero mean and unit variance.
    Input:
    adata: AnnData object with ['n_total_counts_per_cell'] parameter in observations
    
    Output:
    AnnData object
    '''
    sc.pp.regress_out(adata, ['n_total_counts_per_cell'])
    sc.pp.scale(adata)
    return adata

def marker_gene_expression(anndata, marker_dict, gene_symbol_key=None, partition_key='leiden'):
    """
    A function to get mean z-score expressions of marker genes
     
    Inputs:
        anndata         - An AnnData object containing the data set and a partition
        marker_dict     - A dictionary with cell-type markers. The markers should be stores as anndata.var_names or 
                          an anndata.var field with the key given by the gene_symbol_key input
        gene_symbol_key - The key for the anndata.var field with gene IDs or names that correspond to the marker 
                          genes
        partition_key   - The key for the anndata.obs field where the cluster IDs are stored. The default is
                          'louvain_r1' 
    """

    #Test inputs
    if partition_key not in anndata.obs.columns.values:
        print('KeyError: The partition key was not found in the passed AnnData object.')
        print('   Have you done the clustering? If so, please tell pass the cluster IDs with the AnnData object!')
        raise

    if (gene_symbol_key != None) and (gene_symbol_key not in anndata.var.columns.values):
        print('KeyError: The provided gene symbol key was not found in the passed AnnData object.')
        print('   Check that your cell type markers are given in a format that your anndata object knows!')
        raise
        
    if gene_symbol_key:
        gene_ids = anndata.var[gene_symbol_key]
    else:
        gene_ids = anndata.var_names

    clusters = anndata.obs[partition_key].cat.categories
    n_clust = len(clusters)
    marker_exp = pd.DataFrame(columns=clusters)
    marker_exp['cell_type'] = pd.Series({}, dtype='str')
    marker_names = []
    
    z_scores = sc.pp.scale(anndata, copy=True)

    i = 0
    for group in marker_dict:
        # Find the corresponding columns and get their mean expression in the cluster
        for gene in marker_dict[group]:
            ens_idx = np.in1d(gene_ids, gene) #Note there may be multiple mappings
            if np.sum(ens_idx) == 0:
                continue
            else:
                z_scores.obs[ens_idx[0]] = z_scores.X[:,ens_idx].mean(1) #works for both single and multiple mapping
                ens_idx = ens_idx[0]

            clust_marker_exp = z_scores.obs.groupby(partition_key)[ens_idx].apply(np.mean).tolist()
            clust_marker_exp.append(group)
            marker_exp.loc[i] = clust_marker_exp
            marker_names.append(gene)
            i+=1

    #Replace the rownames with informative gene symbols
    marker_exp.index = marker_names

    return(marker_exp)

def gene_expression(anndata, marker_list, gene_symbol_key=None, partition_key='leiden'):
    """A function to get mean z-score expressions of marker genes
     
     Inputs:
        anndata         - An AnnData object containing the data set and a partition
        marker_dict     - A dictionary with cell-type markers. The markers should be stores as anndata.var_names or 
                          an anndata.var field with the key given by the gene_symbol_key input
        gene_symbol_key - The key for the anndata.var field with gene IDs or names that correspond to the marker 
                          genes
        partition_key   - The key for the anndata.obs field where the cluster IDs are stored. The default is
                          'louvain_r1' """

    #Test inputs
    if partition_key not in anndata.obs.columns.values:
        print('KeyError: The partition key was not found in the passed AnnData object.')
        print('   Have you done the clustering? If so, please tell pass the cluster IDs with the AnnData object!')
        raise

    if (gene_symbol_key != None) and (gene_symbol_key not in anndata.var.columns.values):
        print('KeyError: The provided gene symbol key was not found in the passed AnnData object.')
        print('   Check that your cell type markers are given in a format that your anndata object knows!')
        raise
        
    if gene_symbol_key:
        gene_ids = anndata.var[gene_symbol_key]
    else:
        gene_ids = anndata.var_names

    clusters = anndata.obs[partition_key].cat.categories
    n_clust = len(clusters)
    marker_exp = pd.DataFrame(columns=clusters)
    marker_names = []

    i = 0
    
    for gene in marker_list:
        ens_idx = np.in1d(gene_ids, gene) #Note there may be multiple mappings
        if np.sum(ens_idx) == 0:
            continue
        else:
            anndata.obs[ens_idx[0]] = anndata.X[:,ens_idx].mean(1) #works for both single and multiple mapping
            ens_idx = ens_idx[0]

        clust_marker_exp = anndata.obs.groupby(partition_key)[ens_idx].apply(np.mean).tolist()
        marker_exp.loc[i] = clust_marker_exp
        marker_names.append(gene)
        i+=1

    #Replace the rownames with informative gene symbols
    marker_exp.index = marker_names

    return(marker_exp)

def gene_expression_norm(anndata, marker_list, gene_symbol_key=None, partition_key='leiden'):
    """A function to get normalized expressions of marker genes
     
     Inputs:
        anndata         - An AnnData object containing the data set and a partition
        marker_dict     - A dictionary with cell-type markers. The markers should be stores as anndata.var_names or 
                          an anndata.var field with the key given by the gene_symbol_key input
        gene_symbol_key - The key for the anndata.var field with gene IDs or names that correspond to the marker 
                          genes
        partition_key   - The key for the anndata.obs field where the cluster IDs are stored. The default is
                          'louvain_r1' """

    #Test inputs
    if partition_key not in anndata.obs.columns.values:
        print('KeyError: The partition key was not found in the passed AnnData object.')
        print('   Have you done the clustering? If so, please tell pass the cluster IDs with the AnnData object!')
        raise

    if (gene_symbol_key != None) and (gene_symbol_key not in anndata.var.columns.values):
        print('KeyError: The provided gene symbol key was not found in the passed AnnData object.')
        print('   Check that your cell type markers are given in a format that your anndata object knows!')
        raise
        
    if gene_symbol_key:
        gene_ids = anndata.var[gene_symbol_key]
    else:
        gene_ids = anndata.var_names

    clusters = anndata.obs[partition_key].cat.categories
    n_clust = len(clusters)
    marker_exp = pd.DataFrame(columns=clusters)
    marker_names = []
    
    i = 0
    
    ann_data = anndata
    
    for gene in marker_list:
        ens_idx = np.in1d(gene_ids, gene) #Note there may be multiple mappings
        if np.sum(ens_idx) == 0:
            continue
        else:
            ann_data.obs[ens_idx[0]] = ann_data.raw[:, ann_data.var.index].X[:,ens_idx].mean(1) #works for both single and multiple mapping
            ens_idx = ens_idx[0]

        clust_marker_exp = ann_data.obs.groupby(partition_key)[ens_idx].apply(np.mean).tolist()
        marker_exp.loc[i] = clust_marker_exp
        marker_names.append(gene)
        i+=1        
        
    #Replace the rownames with informative gene symbols
    marker_exp.index = marker_names

    return(marker_exp)

#Define cluster score for all markers
def evaluate_partition(anndata, marker_dict, gene_symbol_key=None, partition_key='leiden'):
    ''' This function gives a cell-type score for each partition key (i.e. Leiden clusters)
    Inputs:
    
    anndata         - An AnnData object containing the data set and a partition
    marker_dict     - A dictionary with cell-type markers. The markers should be stores as anndata.var_names or 
                    an anndata.var field with the key given by the gene_symbol_key input
    gene_symbol_key - The key for the anndata.var field with gene IDs or names that correspond to the marker 
                      genes
    partition_key   - The key for the anndata.obs field where the cluster IDs are stored. The default is
                      'louvain_r1'
    Returns:
    
    A dataframe with a score for each cell type.
    '''

    if partition_key not in anndata.obs.columns.values:
        print('KeyError: The partition key was not found in the passed AnnData object.')
        print('   Have you done the clustering? If so, please tell pass the cluster IDs with the AnnData object!')
        raise

    if (gene_symbol_key != None) and (gene_symbol_key not in anndata.var.columns.values):
        print('KeyError: The provided gene symbol key was not found in the passed AnnData object.')
        print('   Check that your cell type markers are given in a format that your anndata object knows!')
        raise
        

    if gene_symbol_key:
        gene_ids = anndata.var[gene_symbol_key]
    else:
        gene_ids = anndata.var_names

    clusters = np.unique(anndata.obs[partition_key])
    n_clust = len(clusters)
    n_groups = len(marker_dict)
    
    marker_res = np.zeros((n_groups, n_clust))
    z_scores = sc.pp.scale(anndata, copy=True)

    i = 0
    for group in marker_dict:
        # Find the corresponding columns and get their mean expression in the cluster
        j = 0
        for clust in clusters:
            cluster_cells = np.in1d(z_scores.obs[partition_key], clust)
            marker_genes = np.in1d(gene_ids, marker_dict[group])
            marker_res[i,j] = z_scores.X[np.ix_(cluster_cells,marker_genes)].mean()
            j += 1
        i+=1

    variances = np.nanvar(marker_res, axis=0)
    if np.all(np.isnan(variances)):
        print("No variances could be computed, check if your cell markers are in the data set.")
        print("Maybe the cell marker IDs do not correspond to your gene_symbol_key input or the var_names")
        raise

    marker_res_df = pd.DataFrame(marker_res, columns=clusters, index=marker_dict.keys())

    #Return the median of the variances over the clusters
    return(marker_res_df)

def silhouette_analysis(range_n_clusters, X):
    '''This function takes as input a matrix X and a list of a range of
    clusters range_n_clusters (that should be from 2 - (n-1) where n is 
    the total number of clusters in the dataset) and yields as output
    a list of the average silhouette score from 100 trials.'''
    scores = []
    for n_clusters in range_n_clusters:
        i=0
        cluster_avg = []
        while i < 100:
            clusterer = KMeans(n_clusters=n_clusters, random_state=10)
            cluster_labels = clusterer.fit_predict(X)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
            silhouette_avg = silhouette_score(X, cluster_labels, metric='cosine')
            cluster_avg.append(silhouette_avg)
            i+=1
        scores.append((n_clusters, cluster_avg))
    return scores

def silhouette_plots(adata, pathway_names, pathway_genes, norm = False, ax=None):
    '''This function gives silhouette scores and boxplots of the scores
    sampled from 100 trials for different numbers of clusters on our 
    pathways.
    
    Input: 
    adata : AnnData object
    pathway_names : string name of the pathways you're evaluating
    pathway_genes : a list of pathway genes.
    norm : whether or not to used z-score or normalized data. z-score is default.
    '''
    range_n_clusters = list(range(2,len(adata.obs['leiden'].unique())))
    if ax == None:
        fig, ax = plt.subplots(figsize=(4,3.85))
    else:
        fig = Figure()
    scores = pd.DataFrame()
    if norm == False:
        df = gene_expression(adata, pathway_genes)
    else:
        df = gene_expression_norm(adata, pathway_genes)
    df = df.T
    df.reset_index(inplace=True)
    X = df[df.columns[1:]].values
    score = silhouette_analysis(range_n_clusters, X)
    scores[pathway_names] = [np.mean(s[1]) for s in score]
    plot = [s[1] for s in score]
    ax.boxplot(plot)
    ax.set_title(pathway_names)
    ax.set_xticklabels(range(2,len(scores)+2))
    ax.tick_params(axis="x", labelsize=10)
    ax.set_xlabel('Number of Clusters')
    ax.set_ylabel('Silhouette Score')
    new_dict = dict(zip(scores.index.values, range_n_clusters))
    scores.rename(new_dict)
    fig.tight_layout()
    return fig, scores

def heatmap(adata, pathway_genes, num_clust, name, norm = False,
                     leg_axes = (1.3, 1.3), leg_cols = 1):
    '''We group the leiden clusters based on similarity of expression of 
    specific genes in a pathway. 
    
    Arguments :
    
    adata : the AnnData gene expression matrix
    pathway_genes : a list of the genes in the pathway
    num_clust : the optimal number of clusters based on silhouette score on
    cosine distance
    name : the name with which we want to label the clusters of this pathway
    leg_axes : we can change the coordinates of the legend
    
    Returns: 
    
    AnnData object labeled with the pathway clusters.
    Return at Index 0: Clustermap Figure you can later save
    Return at Index 1: A dataframe of all the gene expression values
    '''
    if norm:
        df = gene_expression_norm(adata, pathway_genes)
    else:
        df = gene_expression(adata, pathway_genes)
    d = sch.distance.pdist(df.transpose(), metric='cosine')
    L=sch.linkage(d)
    linkage = sch.fcluster(L, num_clust,'maxclust')
    str_linkage = []
    for i in linkage:
        str_linkage.append(str(i))
    new_dict = dict(zip([str(i) for i in range(0,len(str_linkage))], str_linkage))
    adata.obs[name] = adata.obs['leiden'].replace(new_dict)
    if norm:
        df = gene_expression_norm(adata, pathway_genes)
    else:
        df = gene_expression(adata, pathway_genes)
    cols = {}
    for j in list(df.columns):
        cols[j] = str(adata[adata.obs['leiden'] == j].obs[name][0])
    cols=pd.Series(data=cols, name='Clusters')
    labels = adata.obs[name].unique()
    labels = list(map(str, labels))
    cmap = plt.get_cmap('Paired')
    colors = cmap(np.linspace(0, 1, len(labels)))
    lut1 = dict(zip(labels, colors))
    cols_to_return = []
    keys_for_colors = list(lut1.keys())
    keys_for_colors.sort()
    for k in keys_for_colors:
        cols_to_return.append(lut1[k])
    adata.uns[name+'_colors'] = cols_to_return
    row_colors1 = cols.map(lut1)
    g = sns.clustermap(df, metric='cosine', row_cluster=False, cmap='viridis',
                      col_linkage=L, col_colors=row_colors1,figsize=(6,6));
    ax = g.ax_heatmap
    legend_elements=[]
    keys= list(lut1.keys())
    keys.sort()
    for j in keys:
        legend_elements.append(Line2D([0], [0], marker = 's', 
                                      label = j, color = lut1[j]))
    #ax.legend(handles = legend_elements, title = 'Clusters', fontsize='small',
              #loc='upper right', bbox_to_anchor=(1.4,1.3), ncol=leg_cols)

    g.fig.suptitle((name + ' with ' + str(num_clust) + ' clusters'), y=1.0,x=0.5,fontsize='large') 
    ax.set_xlabel('Leiden clustering', x=0.5)
    return g.fig, df

def exp_across_clusters(df):
    '''Plots a bar chart of expression summed across different Leiden clusters.
    This is useful to visualize which clusters we can remove from our heatmaps.
    
    Input: A dataframe whose rows are different genes, and whose columns are the Leiden
    cluster labels.
    
    Output: a bar chart with the total expression for each cluster.'''
    
    cols = df.sum(axis=0)
    return cols.plot.bar(title='Z-Score Expression Sum Across Leiden Clusters');

def exp_across_genes(df):
    '''Plots a bar chart of expression summed across different genes.
    This is useful to visualize which genes we can remove from our heatmaps.
    
    Input: 
    df: A dataframe whose rows are different genes, and whose columns are the Leiden
    cluster labels.
    
    Output: a bar chart with the total expression for each cluster.'''
    
    rows = df.sum(axis=1)
    return rows.plot.bar(title='Z-Score Expression Sum Across Genes');

def exp_above_threshold(df, axis, threshold):
    '''This function visualizes counts of gene expression above a certain threshold.
    Inputs: 
    
    df: A pandas dataframe whose rows are genes, and whose columns are Leiden 
    cluster labels.
    axis: 0=rows, 1=columns
    threshold: define a threshold value for gene expression.
    
    Output: a bar chart
    '''
    ax = df.iloc[:,:].ge(threshold).sum(axis).plot.bar(
        title=('Histogram of Counts Above Threshold Value'))
    ax.set_xlabel('Sum along axis '+ str(axis))
    ax.set_ylabel('Counts Above Threshold Value')
    return ax

