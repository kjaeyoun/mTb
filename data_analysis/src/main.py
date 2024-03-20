import csv
import os
import math
import numpy as np
import pandas as pd
from utils import bioid, safe_utils

SAINT_FILE = "../data/Supplementary Table 1 - 291 (B) SAINT and MiST analysis for mTb.txt"
NETWORK_INPUT_FILE = "../results/networks/network_input.txt"
NETWORK_FILE_CYJS = "../results/networks/mTb_network.cyjs"
NETWORK_FILE_EDGE = "../results/networks/mTb_network default edge.csv"
GMT_FILE_CC = "../data/gprofiler_hsapiens.name/hsapiens.GO_CC.name.gmt"
SAFE_OUTPUT = "../results/safe/"

BFDR_THRESHOLD = 0.01
STRING_THRESHOLD = 0.7

NODE_COLORS = {0: "#808080", 1: "#F2D56F", 2: "#E71733", 3: "#3690C0", 4: "#41AB5D",
               5: "#705592", 6: "#E797A5", 7: "#66B2A2", 8: "#f9a100", 9: "#003679"}


class SafeAnalysisParams:
    def __init__(self, saint_file, network_input_file, cyjs_network_file, network_edge, gmt_file, safe_output, bfdr_threshold=0.01, string_threshold=0.7):
        self.saint_file = saint_file
        self.network_input_file = network_input_file
        self.cyjs_network_file = cyjs_network_file
        self.network_edge = network_edge
        self.gmt_file = gmt_file
        self.safe_output = safe_output
        self.safe_node_attribute = os.path.join(
            safe_output, "node_properties_annotation.txt")
        self.bfdr_threshold = bfdr_threshold
        self.string_threshold = string_threshold
        self.output_file_pdf = os.path.join(safe_output, "SAFE_results.pdf")
        return None


def min_max_normalization(value, min_value, max_value):
    """Perform min-max normalization"""
    return (value - min_value) / (max_value - min_value)


def create_weights(df_network):
    """Create weights based on normalization log-scaled AvgSpec and STRING scores."""
    # AvgSpec을 log10 혹은 log2 스케일로 변환해야할 것 같음.
    # 각각의 Bait에서의 결과에서 scale을 정해야할 것 같음.
    # Spectral counts와 String score를 서로 결합할 방법을 고민해야 할 것 같음.
    # 아래 코드는 다 지우고 다시 작성해야 할 것으로 보임.
    df_network["logAvgSpec"] = np.log10(
        df_network['AvgSpec'] + 1)  # +1 to avoid log(0)
    df_network['NormalizedLogAvgSpec'] = df_network['logAvgSpec'].apply(
        min_max_normalization, args=(df_network['logAvgSpec'].min(), df_network['logAvgSpec'].max()))
    df_network['NormalizedStringScore'] = df_network['STRING score'].apply(
        min_max_normalization, args=(df_network["STRING score"].min(), df_network["STRING score"].max()))
    df_network["Weights"] = (
        df_network['NormalizedLogAvgSpec'] + df_network['NormalizedStringScore']) / 2
    return None




def gene_ontology_from_safe():
    return None


def main():
    params = SafeAnalysisParams(
        saint_file=SAINT_FILE,
        network_input_file=NETWORK_INPUT_FILE,
        network_edge=NETWORK_FILE_EDGE,
        cyjs_network_file=NETWORK_FILE_CYJS,
        gmt_file=GMT_FILE_CC,
        safe_output=SAFE_OUTPUT,
        bfdr_threshold=0.01,
        string_threshold=0.7)

    df_saint = bioid.load_mTb_dataset(SAINT_FILE)

    safe_utils.run_safe_analysis(params)
    safe_utils.define_node_attributes(params, df_saint)
    safe_utils.define_edge_attributes(params)
    return None


if __name__ == "__main__":
    # To load mTb dataset: Columns for BaitGene and Cell cycle are created during loading.
    df_saint = bioid.load_mTb_dataset(SAINT_FILE)

    # To pool the dataset based on cell cycle
    df_pool = bioid.pool_saint_results(df_saint)

    # To generate input for network input
    main()


params = SafeAnalysisParams(
    saint_file=SAINT_FILE,
    network_input_file=NETWORK_INPUT_FILE,
    cyjs_network_file=NETWORK_FILE_CYJS,
    network_edge=NETWORK_FILE_EDGE,
    gmt_file=GMT_FILE_CC,
    safe_output=SAFE_OUTPUT,
    bfdr_threshold=0.01,
    string_threshold=0.7)

df_edge = safe_utils.define_edge_attributes(params)
df_node = safe_utils.define_node_attributes(params, df_saint)


def retrive_gmt_file(gmt_file, terms):
    output_dict = {}
    with open(gmt_file, "r", newline="") as tsv_input:
        file_read = csv.reader(tsv_input, delimiter="\t")
        
        for term in terms:
            for row in file_read:
                go_name = row[1]
                go_genes = row[2:]
    
                if term in go_name:
                    output_dict[go_name] = go_genes

    return output_dict


def merge_avg_spec():
    return None


significant_preys = df_pool.PreyGene.drop_duplicates().tolist()

df_saint = df_saint[df_saint.PreyGene.isin(significant_preys)]
df_saint = df_saint[["BaitGene", "CellCycle", "PreyGene", "AvgSpec"]]
df_g1 = df_saint[(df_saint.CellCycle == "G1")].groupby("PreyGene").sum()
df_s = df_saint[(df_saint.CellCycle == "S")].groupby("PreyGene").sum()
df_g2 = df_saint[(df_saint.CellCycle == "G2")].groupby("PreyGene").sum()
df_m = df_saint[(df_saint.CellCycle == "M")].groupby("PreyGene").sum()


df_g1 = df_g1.rename(columns={"AvgSpec": "G1"})
df_s = df_s.rename(columns={"AvgSpec": "S"})
df_g2 = df_g2.rename(columns={"AvgSpec": "G2"})
df_m = df_m.rename(columns={"AvgSpec": "M"})


df = pd.merge(df_g1, df_s, left_index=True, right_index=True)
df = pd.merge(df, df_g2, left_index=True, right_index=True)
df = pd.merge(df, df_m, left_index=True, right_index=True)


log2_fc = []
for i in range(len(df)):
    row = df.iloc[i]
    log2_fc.append(math.log2(max(row) / min(row)))
df["log2FC"] = log2_fc


go_cc = retrive_gmt_file(params.gmt_file, ["centrosom", "centriol", "pericentriolar"])
go_bp = retrive_gmt_file(
    "../data/gprofiler_hsapiens.name/hsapiens.GO_BP.name.gmt",
    ["centrosom", "centriol", "pericentriol"])

centrosomal_proteins = [y for x in go_cc.values() for y in x ]
bp = [y for x in go_bp.values() for y in x ]
centrosomal_proteins = set(centrosomal_proteins + bp)

df = df[df.index.isin(centrosomal_proteins) == False]

df["Max"] = df.sum(axis=1)

df.to_csv("../results/novel/test.txt", sep="\t")
