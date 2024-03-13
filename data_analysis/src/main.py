import matplotlib.pyplot as plt
import networkx as nx
import os
import pandas as pd
from utils.string import fetch_STRING_db, merge_STRING_into_BioID


SAINT_FILE = "../data/Supplementary Table 1 - 291 (B) SAINT and MiST analysis for mTb.txt"
NETWORK_INPUT_FILE = "../results/network_input.txt"
BFDR_THRESHOLD = 0.01
STRING_THRESHOLD = 0.7


def load_mTb_dataset(file_path):
    df_saint = pd.read_csv(file_path, sep="\t", index_col=0)
    df_saint = df_saint.reset_index(drop=True)
    columns = list(df_saint)
    df_saint[["BaitGene", "CellCycle"]
             ] = df_saint.Bait.str.split("_", expand=True)
    columns.insert(1, "BaitGene")
    columns.insert(2, "CellCycle")
    df_saint = df_saint[columns]

    return df_saint


def pool_saint_results(df):
    df = df.sort_values(by=["AvgSpec"])
    df_pool = df.drop_duplicates(subset=["BaitGene", "PreyGene"], keep='last')

    return df_pool


def generate_network(df_pool, output_file, bfdr=0.01, string_score=0.4):
    df = df_pool[["BaitGene", "PreyGene", "AvgSpec", "BFDR"]]
    df = df[df.BFDR <= bfdr]

    # To fetch STRING database based on the pooled results.
    bait_list = df.BaitGene.drop_duplicates().tolist()
    prey_list = df.PreyGene.drop_duplicates().tolist()
    gene_list = set(bait_list + prey_list)

    df_string = fetch_STRING_db(gene_list, string_score)
    df_network = merge_STRING_into_BioID(df, df_string, output_file)
    return df_network


def min_max_normalization(value, min_value, max_value):
    return (value - min_value) / (max_value - min_value)


def create_weights():
    return None


if __name__ == "__main__":
    # To load mTb dataset: Columns for BaitGene and Cell cycle are created during loading.
    df_saint = load_mTb_dataset(SAINT_FILE)

    # To pool the dataset based on cell cycle
    df_pool = pool_saint_results(df_saint)

    # To generate input for network input
    if not os.path.exists(NETWORK_INPUT_FILE):
        df_network = generate_network(df_pool, output_file=NETWORK_INPUT_FILE,
                                      bfdr=BFDR_THRESHOLD, string_score=STRING_THRESHOLD)
    else:
        df_network = pd.read_csv(NETWORK_INPUT_FILE, sep="\t")

    import numpy as np
    
    # -1을 
    df_network["logAvgSpec"] = np.log10(df_network['AvgSpec'])
    min_avg_spec = df_network.AvgSpec.min()
    max_avg_spec = df_network.AvgSpec.max()
    min_string_score = df_network["STRING score"].min()
    max_string_score = df_network["STRING score"].max()
    df_network['NormalizedAvgSpec'] = df_network['AvgSpec'].apply(
        min_max_normalization, args=(min_avg_spec, max_avg_spec))
    
    df_network['NormalizedStringScore'] = df_network['STRING score'].apply(
        min_max_normalization, args=(min_string_score, max_string_score))
    
    df_network["Weights"] = (df_network['NormalizedAvgSpec'] + df_network['NormalizedStringScore']) / 2

    G = nx.Graph()

    for i in range(len(df_network)):
        row = df_network.iloc[i]
        bait_gene, prey_gene = row[0], row[1]
        G.add_edge(bait_gene, prey_gene)

    plt.figure(figsize=(10, 8))
    nx.draw(G, node_color='lightblue', edge_color='gray')
    plt.title("Gene Interaction Graph")
    plt.show()
