import numpy as np
import pandas as pd
import requests
import matplotlib.pyplot as plt
import seaborn as sns
import re
import requests
import scipy.stats as stats
import networkx as nx


def parse_features(features, row_idx=0):
    '''A function for parsing "feature" (from genome annotation file) column in files generated by `bedtools intersect` command
    
    Parameters:
        features (iterable): list or pd.Series (column) with features
        row_idx (int): index of reference row with items from "feature" column
        
    Returns:
        features_df (pd.DataFrame): pd.DataFrame with parsed features
    '''
    
    
    features_dict = {}
    a = []
    for item in features.iloc[row_idx].split(";"):
        key, value = item.split("=")
        features_dict[key] = a

    for i in range(features.shape[0]): 
        for row in features.iloc[i].split(";"):
            row = row.split(";")
            for el in row:
                k, v = el.split("=")
                if k in features_dict.keys():
                    if features_dict[k]:
                        features_dict[k].append(v)
                    else:
                        features_dict.update({k: [v]})

    features_df = pd.DataFrame(features_dict)
    return features_df


def parse_string(s: str) -> str:
    '''A function for compacting a string on a graph by inserting "\n"
    
    Parameters:
        s (str): string to compact
    
    Returns:
        string with inserted "\n"
    '''
    
    
    split = s.split()
    if len(split) <= 2:
        return s
    return f"{split[0]} {split[1]}" + "\n" + " ".join(split[2:])


def generate_1kb_regions(dataframe, coord_file):
    '''A function for generating 1kb region around SNPs for genome assembly remapping in chrN:start-stop format and 
    writing these coordinates in file for genome assembly remapping
    
    Parameters:
        dataframe (pd.DataFrame): dataframe with SNP coordinates
        coord_file (str): prefix for output file with generated coordinates
        
    Returns:
        snps_unique_df (pd.DataFrame): dataframe containing +-500bp coordinates around SNPs
    '''
    
    
    # list with chromosomes in format from MIDESP output
    chrom_ = [f"chr_0{i}" if i < 10 else f"chr_{i}" for i in range(1, 21)]
    
    # list with chromosomes in format for remapping
    chrom_num = [f"chr{i}" for i in range(1, 21)]

    # dictionary with chromosomes for renaming
    chrom = dict(zip(chrom_, chrom_num))
    
    # combining all SNPs from pairs
    unique_snps = pd.concat([dataframe["SNP1"], dataframe["SNP2"]], axis=0).unique()
    
    # generating dataframe for coordinates
    snps_unique_df = pd.DataFrame(unique_snps, columns=["SNP_ID"])
    
    # parsing SNP IDs to obtain chr and position
    snps_unique_df[["chr", "position"]] = snps_unique_df["SNP_ID"].str.split(":", expand=True)
    
    # rename chromosomes in format for remapping
    snps_unique_df.loc[:, "chr"] = snps_unique_df["chr"].apply(lambda x: chrom[x])
    
    # generating coordinates +-500bp from SNP position for remapping on genome assembly v2.1
    snps_unique_df.loc[:, "position-500"] = snps_unique_df["position"].astype(int) - 500
    snps_unique_df.loc[:, "position+500"] = snps_unique_df["position"].astype(int) + 500
    pos_columns = ["position-500", "position+500"]
    snps_unique_df.loc[:, "coordinate"] = snps_unique_df[pos_columns].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
    snps_unique_df.loc[:, "coordinate"] = snps_unique_df[["chr", "coordinate"]].apply(
        lambda row: ':'.join(row.values.astype(str)), axis=1
    )    
    
    # writing coordinates in file
    snps_unique_df.to_csv(
        f"{coord_file}.txt", sep="\n", 
        columns=["coordinate"], header=False, index=False
    )
    return snps_unique_df


def generate_bed_file_aa(dataframe, bed_path):
    '''A function for replacing chromosome names with RefSeq IDs used in genome assembly v2.1 annotation for amino acide content dataset
    and generating .bed file for intersection with annotation
    
        This file containes four columns:
        1) chromosome RefSeq ID;
        2) start coordinate;
        3) stop coordinate;
        4) SNP ID
    
    Parameters:
        dataframe (pd.DataFrame): dataframe with SNP region coordinates
        chr_names_dict (dict): dictionary with chromosome names and respective IDs for renaming
        bed_path (str): prefix for output .bed file
        
    Returns:
        dataframe (pd.DataFrame): initial dataframe with replaced chromosome names and .bed coordinates
    '''
    
    # getting chromosome names in genome assembly v2.1
    request = requests.get("https://www.ncbi.nlm.nih.gov/assembly/GCF_000004515.5")

    # url contains a few tables, pd.read_html is really useful in this case
    tables = pd.read_html(request.content)

    # chromosome IDs in forth table
    chrom_names_v21 = tables[3]["RefSeq sequence"]

    # list with v2.1 chromosome IDs, the last element is nan, skipping it
    chrom_names_v21 = list(chrom_names_v21[:-1])
    
    # list with chromosomes in format for remapping
    chrom_num = [f"chr{i}" for i in range(1, 21)]

    # dictionary with "chrN" and respective chromosome RefSeq IDs in genome assembly v2.1
    chr_names = dict(zip(chrom_num, chrom_names_v21))
    
    # replacing chromosome names with RefSeq IDs
    dataframe["chr_ID"] = dataframe.loc[:, "mapped_id"].replace(chr_names)
    
    # generating .bed coordinates
    dataframe[["mapped_start", "mapped_stop"]] = dataframe[["mapped_start", "mapped_stop"]].astype(int) - 1
    
    # generating .bed file
    dataframe.to_csv(f"{bed_path}.bed", sep='\t',
                      columns=["chr_ID", "mapped_start", "mapped_stop", "SNP_IDs"], 
                      header=False, index=False)
    return dataframe


def generate_bed_file_unique(dataframe, bed_path):
    '''A function for generating .bed file for intersection with genome assembly v2.1 annotation from 
    EnsemblePlants for commercial dataset
    
        This file containes four columns:
        1) chromosome name as in annotation;
        2) start coordinate;
        3) stop coordinate;
        4) SNP ID
    
    Parameters:
        dataframe (pd.DataFrame): dataframe with SNP region coordinates
        bed_path (str): prefix for output .bed file
        
    Returns:
        unique_snps_dataframe (pd.DataFrame): dataframe with replaced chromosome names and .bed coordinates for unique SNPs from all pairs
    '''
    
    
    # combining all SNPs from pairs
    unique_snps = pd.concat([dataframe["SNP1"], dataframe["SNP2"]], axis=0).unique()
    
    # generating dataframe with unique SNP IDs
    unique_snps_df = pd.DataFrame(unique_snps, columns=["SNP_ID"])
    
    # generating columns with coordinates parsing SNP IDs
    unique_snps_df["chr"] = unique_snps_df["SNP_ID"].apply(lambda x: re.split(r':', x)[0])
    unique_snps_df["pos_stop"] = unique_snps_df["SNP_ID"].apply(lambda x: re.split(r':', x)[1])
    
    # renaming chromosome as in annotation file (just numbers or scaffold names)
    unique_snps_df["chr"] = unique_snps_df["chr"].apply(
        lambda x: x[3:] if x.startswith("chr") else x
    )
    unique_snps_df["pos_start"] = unique_snps_df["pos_stop"].astype(int) - 1
    
    # generating .bed file with unique SNPs
    unique_snps_df.to_csv(
        f"{bed_path}.bed", sep='\t',
        columns=["chr", "pos_start", "pos_stop", "SNP_ID"],
        header=False, index=False)
    return unique_snps_df
    

def assign_gene_for_snp(dataframe, n_snps, snp_gene_dict):
    '''A function for assigning gene ids for SNPs
    
    Parameters:
        dataframe (pd.DataFrame): dataframe with SNP pairs
        n_snps (int): number of interacting SNPs
        snp_gene_dict (dict): dictionary with SNP-gene pairs
        
    Returns:
        pairs_df (pd.DataFrame): dataframe with SNP and respective gene pairs
    '''
    
    
    for idx in range(n_snps):
        # from SNPs to genes
        dataframe[f"gene{idx+1}"] = dataframe[f"SNP{idx+1}"].apply(lambda x: snp_gene_dict[x] if x in snp_gene_dict.keys() else np.nan)    
        
    # removing pairs/triplets for which respective genes were not found
    pairs_df = dataframe.dropna(axis=0)
    
    # transforming each gene pair/triplet into frozenset to filter mirrored pairs and 
    # those pairs/triplets with genes interacting with themselves
    pairs_df["setted"] = [
        frozenset(i) for i in zip(pairs_df["gene1"], pairs_df["gene2"])
    ]
    
    # sorting dataframe by "MI_APC value
    pairs_df = pairs_df.sort_values(by="MI_APC", ascending=False)
    
    # dropping duplicated gene-gene pairs/triplets with less "MI_APC" value
    pairs_df = pairs_df.drop_duplicates(subset="setted", keep="first").reset_index(drop=True)
    
    # dropping pairs/triplets with genes interacting with themselves
    # length of a set < 2
    pairs_df = pairs_df[pairs_df["setted"].apply(lambda x: len((x)) >= 2)]
    pairs_df = pairs_df.drop(columns=["setted"])
    return pairs_df


def parse_gene_info(gene_info_df, column):
    '''A function for parsing file with Gene Info results generated for convertion of gene names starts with "GLYMA" to respective Gene IDs using 
    bioloigical DataBase network and its utility for database to database conversions ([db2db](https://biodbnet.abcc.ncifcrf.gov/db/db2db.php))
    
    Parameters:
        gene_info_df (pd.DataFrame): dataframe with Gene Info conversion results
        column (str): column name
        
    Returns: 
        gene_info_df (pd.DataFrame): initial dataframe with new colomn containing parsed Gene IDs
    '''
    
    
    gene_info_df["gene_ids"] = gene_info_df[column].apply(lambda x: x.split()[0])
    gene_info_df = gene_info_df[gene_info_df["gene_ids"] != "-"]
    return gene_info_df


def parse_kegg_results(kegg_dataframe, dataset_name):
    '''A function for parsing KEGG enrichment results obtained via ShinyGo v0.77
    
    Parameters:
        kegg_dataframe (pd.DataFrame): dataframe with KEGG enrichment results
        dataset_name (str): name of analysed dataset to display on plot
        
    Returns:
        kegg_dataframe (pd.DataFrame): initial dataframe with modified and additional columns for visualization
    '''
    
    
    # parsing KEGG Pathway url to obtain Pathway ID
    kegg_dataframe.loc[:, "KEGG_ID"] = kegg_dataframe["URL"].str.split("?", expand=True)[1]

    # calculating -log10(FDR)
    kegg_dataframe.loc[:, "Enrichment FDR"] = kegg_dataframe["Enrichment FDR"].apply(lambda x: round(-np.log10(x), 1))

    # generating strings with Pathway name and its ID
    kegg_dataframe.loc[:, "KEGG_ID"] = kegg_dataframe[["Pathway", "KEGG_ID"]].agg(": ".join, axis=1)

    # inserting "\n" into long strings for prettier visualization
    kegg_dataframe.loc[:, "KEGG_ID"] = kegg_dataframe["KEGG_ID"].apply(parse_string)
    
    # barplot with KEGG Enrichment analysis results for Alanine content dataset
    fig, ax = plt.subplots(figsize=(8,8), dpi=300)
    sns.barplot(data=kegg_dataframe, x="Fold Enrichment", y="KEGG_ID", hue="Enrichment FDR", palette="RdBu", ax=ax, dodge=False)
    ax.set_ylabel("")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, ncol=4, loc="lower right", 
              title="$\mathbf{-log_{10}(FDR)}$")
    ax.set_title(f"Results of KEGG enrichment analysis for all gene\nfound for {dataset_name}",
                fontsize=12, fontweight="bold")
    return kegg_dataframe


def plot_graph(graph, edge_attr, title, pos, labels=None):
    '''A function for visualize a network with all found interacting genes during MIDESP analysis
    
    Parameters:
        graph (nx.classes.graph.Graph): graph to visualize
        edge_attr (list): list with attributes for edges
        title (str): plot title
        pos (func): networkx layout
        labels (dict): dictionary with nodes and their names to display on network
        
    Returns: None
    '''
    
    
    fig, ax = plt.subplots(figsize=(20,20), dpi=300)

    # to visualize graph with specified layout
    pos = pos(graph)    
    
    if labels:
        # drawing network nodes
        nx.draw_networkx_nodes(graph, pos, node_size=200, alpha=0.6, ax=ax, label=False)
        # drawing nodes labels
        nx.draw_networkx_labels(graph, pos, labels, font_size=14, font_color='r')
    else:
        # drawing network nodes
        nx.draw_networkx_nodes(graph, pos, node_size=25, ax=ax, label=False)
    
    # drowing network edges
    # MI_APC value displayed by edge width
    nx.draw_networkx_edges(graph, pos, width=edge_attr, ax=ax)

    ax.set_title(f"{title}", fontsize=20, fontweight="bold")

    
def filter_graph_by_degree(graph, return_filt_graph=False):
    '''A function to filter graph by z-score of node degrees
    
    Parameters:
        graph (nx.classes.graph.Graph): graph for filtering
        return_filt_graph (bool): True - to return filtered both resulted dataframe and filtered graph
        
    Returns:
        degrees_df_filt (pd.DataFrame): dataframe with filtered input graph by z-score of nodes degrees
        filtered_G (nx.classes.graph.Graph): optional, filtered graph
    '''
    
    
    # dataframe with graph nodes and their degrees  
    degrees_df= pd.DataFrame(graph.degree, columns=["nodes","degrees"])
    
    # calculating z-scores of node degrees
    degrees_df["z_scores"] = stats.zscore(degrees_df["degrees"], ddof=1)
    
    # keeping nodes with z-score of degrees >= 3
    degrees_df = degrees_df[np.abs(degrees_df["z_scores"]) >= 3].sort_values(by="z_scores", ascending=False)
    
    # subgraph with filtered nodes by z-score of degrees
    filtered_G = graph.subgraph(degrees_df["nodes"]).copy()
    
    # dataframe with filtered nodes and their degrees
    degrees_df_filtered = pd.DataFrame(dict(filtered_G.degree()).items(), columns=["nodes", "degrees"])
    
    # filtering dataframe by degrees for further extracting top nodes
    degrees_df_filtered.sort_values(by="degrees", ascending=False, inplace=True)
    
    # removing nodes without edges from dataframe and filtered graph
    degrees_df_filt = degrees_df_filtered[degrees_df_filtered["degrees"] >= 1].reset_index(drop=True)
    
    if return_filt_graph:
        # removing nodes without edges from filtered graph
        nodes_to_drop = degrees_df_filtered[degrees_df_filtered["degrees"] == 0]["nodes"]
        filtered_G.remove_nodes_from(nodes_to_drop)
        return degrees_df_filt, filtered_G
    else:
        return degrees_df_filt
    
def filter_graph_by_MIAPC(genes_df):
    '''A function for calculating z-score of MI_APC values, filtering gene pairs with z-score of MI_APC >= 3, and
    returning top10 gene pairs by these score
    
    Parameters:
        genes_df (pd.DataFrame): dataframe with gene pairs and their MI_APC values
        
    Returns:
        genes_df (pd.DataFrame): filtered input dataframe with new column "z-scores"
        top_pairs (list): list with top 10 gene pairs by MI_APC z-score
    '''
    
    
    # calculating z-score of MI_APC
    genes_df["z_scores"] = stats.zscore(genes_df["MI_APC"], ddof=1)
    
    # keeping nodes with z-score of degrees >= 3
    genes_df = genes_df.sort_values(by="z_scores", ascending=False)
    genes_df = genes_df[genes_df["z_scores"] >= 3]
    
    # finding top 10 gene pairs by z-score MI_APC
    top_pairs = list(pd.concat([genes_df["gene1"][:10], genes_df["gene2"][:10]], axis=0))
    return genes_df, top_pairs