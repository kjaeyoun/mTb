import pandas as pd
from utils.string_db import fetch_STRING_db, merge_STRING_into_BioID


mTb_gene_correct = {"FGFR1OP": "CEP43", "MPP5": "PALS1", 
                    "FAM21A": "WASHC2A", "FAM21B": "WASHC2A"}


def load_mTb_dataset(file_path):
    '''Load and preprocess the SAINT dataset'''
    # Load dataset and directly split 'Bait' column into 'BaitGene' and "CellCycle'
    df_saint = pd.read_csv(file_path, sep="\t", index_col=0)
    df_saint[["BaitGene", "CellCycle"]
             ] = df_saint['Bait'].str.split("_", expand=True)

    # Reorder columns to place "BaitGene" and "CellCycle" after "Bait"
    columns = ["Bait", "BaitGene", "CellCycle"]
    columns += [col for col in df_saint.columns if col not in columns]
    df_saint = df_saint[columns]
    
    # To correct the gene name
    corrected_gene_names = [mTb_gene_correct[x] if x in mTb_gene_correct else x for x in df_saint.PreyGene]
    df_saint["PreyGene"] = corrected_gene_names
    return df_saint


def pool_saint_results(df, bfdr=0.01):
    '''Pool SAINT results by keeping the last duplicate based on AvgSpec'''
    baits = df.BaitGene.drop_duplicates().tolist()
    df_pool = None
    for bait in baits:
        temp_df = df[(df.BaitGene == bait) & (df.BFDR <= bfdr)]
        # Sorting and deduplilcation
        temp_df = temp_df.sort_values(by=["AvgSpec"]).drop_duplicates(
            subset=["BaitGene", "PreyGene"], keep="last")
        
        if df_pool is None:
            df_pool = temp_df.copy()
        else:
            df_pool = pd.concat([df_pool, temp_df])
    

    return df_pool


def generate_network_input(df_pool, output_file, bfdr=0.01, string_score=0.4):
    '''Generate network input by filtering based on BFDR and merging STRING database'''
    
    # Filter based on BFDR then select required columns
    selected_columns = ["BaitGene", "PreyGene", "AvgSpec", "BFDR"]
    df_filtered = df_pool[df_pool['BFDR'] <= bfdr][selected_columns]
    
    # Fetch STRING database information
    gene_list = set(df_filtered["BaitGene"].tolist() + df_filtered["PreyGene"].tolist())
    df_string = fetch_STRING_db(gene_list, string_score)
    
    # Merge STRING data into BioID data
    df_network = merge_STRING_into_BioID(df_filtered, df_string)
    
    # To correct the gene name
    corrected_gene_a = [mTb_gene_correct[x] if x in mTb_gene_correct else x for x in df_network.Protein_A]
    corrected_gene_b = [mTb_gene_correct[x] if x in mTb_gene_correct else x for x in df_network.Protein_B]

    df_network["Protein_A"] = corrected_gene_a
    df_network["Protein_B"] = corrected_gene_b
    
    # To save the network
    df_network.to_csv(output_file, sep="\t", index=None)
    return df_network



if __name__ == "__main__":
    SAINT_FILE = "../../data/Supplementary Table 1 - 291 (B) SAINT and MiST analysis for mTb.txt"
    df_saint = load_mTb_dataset(SAINT_FILE)
    
    df_pool = pool_saint_results(df_saint)
    
    sig_preys_pool = df_pool[df_pool.BFDR <= 0.01]["PreyGene"].drop_duplicates()
    sig_preys_saint = df_saint[df_saint.BFDR <= 0.01]["PreyGene"].drop_duplicates()
    
    if len(sig_preys_pool) != len(sig_preys_saint):
        print("There is a problem in pooling")
    else:
        print("BioID results are pooled correctly.")