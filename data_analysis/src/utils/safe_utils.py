import csv
import json
import math
import numpy as np
import networkx as nx
import os
import pandas as pd 
from . import bioid

from safepy import safe

NODE_COLORS = {0: "#808080", 1: "#F2D56F", 2: "#E71733", 3: "#3690C0", 4: "#41AB5D",
               5: "#705592", 6: "#E797A5", 7: "#66B2A2", 8: "#f9a100", 9: "#003679"}



# This function makes a gpickle file from a cyjs file.
def cyjs_to_gpickle(cyjs_file):
    # Generate the output file name by replacing .cyjs with .gpickle
    output_file = os.path.join(os.path.dirname(cyjs_file), 
                               os.path.basename(cyjs_file).replace(".cyjs", ".gpickle"))
    
    if os.path.exists(output_file):
        G = nx.read_gpickle(output_file)
        return G, output_file
    
    # Load the JSON data from the file
    with open(cyjs_file, "r") as json_file:
        json_data = json.load(json_file)
    
    # Extract nodes and edges from the JSON data
    nodes = json_data["elements"]["nodes"]
    edges = json_data["elements"]["edges"]
    
    # Create a list of nodes with their attributes
    node_list = [(int(node["data"]["id"]), {"label": node["data"]["shared_name"],
                                            "label_orf": node["data"]["name"],
                                            "x": node["position"]["x"],
                                            "y": node["position"]["y"]}) for node in nodes]
    
    # Create a dictionary mapping original node IDs to new node IDs
    node_index = {node_id: index for index, (node_id, _) in enumerate(node_list)}
    
    # Update the node IDs in the node list
    node_list = [(node_index[node_id], attrs) for node_id, attrs in node_list]
    
    # Create a graph and add the nodes
    G = nx.Graph()
    G.add_nodes_from(node_list)

    # Add the edges to the graph
    for edge in edges:
        source, target = int(edge["data"]["source"]), int(edge["data"]["target"])
        source_cord, target_cord = G.nodes[node_index[source]], G.nodes[node_index[target]]
        length = math.hypot(source_cord["x"] - target_cord["x"], source_cord["y"] - target_cord["y"])
        G.add_edge(node_index[source], node_index[target], weight=1, length=length)
    
    # Write the graph to a file and return it
    nx.write_gpickle(G, output_file)
    return G, output_file
    
    

# Function to create an annotation file for SAFE analysis from GMT files
def gmt_to_annotation(gmt_file):
    # Generate the output file name by replacing .gmt with .txt
    output_file = os.path.basename(gmt_file).replace(".gmt", ".txt")
    output_file = os.path.join(os.path.dirname(gmt_file), output_file)
    
    if os.path.exists(output_file):
        print("You already have an annotation file.")
        return output_file

    # Parsing the GMT file to collect gene lists and terms
    gene_list = []
    term_dict = {}
    with open(gmt_file, "r", newline="", encoding="utf-8") as tsv_input:
        file_reader = csv.reader(tsv_input, delimiter="\t")
        for row in file_reader:
            go_term = row[1]
            genes = row[2:]
            gene_list.extend(genes)
            term_dict[go_term] = genes
        gene_list = list(set(gene_list))

    # Writing the ananotation file
    with open(output_file, "w", newline="", encoding="utf-8") as tsv_output:
        file_writer = csv.writer(tsv_output, delimiter="\t")
        term_list = list(term_dict.keys())
        header = ["ORF"] + term_list  # Prepare header with ORF and GO terms
        file_writer.writerow(header)
        
        # Initialize progress counter for user feedback
        total_row = len(gene_list)
        counter = 0
    
        for gene in gene_list:
            counter += 1
            row = [gene] + [1 if gene in term_dict[term] else 0 for term in term_list]
            file_writer.writerow(row)
            
            # Print progress update
            progress = round((counter / total_row) * 100, 2)
            print(f"\rMaking an annotation file for SAFE analysis: {progress}% completed.", end="")
    
    print(f"\nattribute file of {os.path.basename(output_file)} has been made")
    return output_file



def run_safe_analysis(params):
    """Run SAFE analysis pipeline."""
    if not os.path.exists(params.safe_output):
        os.makedirs(params.safe_output)

    OUTPUT_FILE_PDF = os.path.join(params.output_file_pdf)
    if os.path.exists(OUTPUT_FILE_PDF):
        print("SAFE analysis has already been performed.")
        return

    df_saint = bioid.load_mTb_dataset(params.saint_file)
    df_pool = bioid.pool_saint_results(df_saint)

    if not os.path.exists(params.network_input_file):
        bioid.generate_network_input(
            df_pool, params.network_input_file, params.bfdr_threshold, params.string_threshold)

    # If the CYJS file does not exist, tell the user how to create the CYJS file
    if not os.path.exists(params.cyjs_network_file):
        print("CYJS network file not found. Please use Cytoscape and your network_input_file to create a CYJS network file named 'mTb_network.cyjs' and place it in:", params.cyjs_network_file)
        print("Refer to the provided guide for creating a CYJS file using Cytoscape.")
        return

    G, gpickle_network_file = cyjs_to_gpickle(
        params.cyjs_network_file)
    annotation_file = gmt_to_annotation(params.gmt_file)

    sf = safe.SAFE()
    sf.multiple_testing = True
    sf.output_dir = params.safe_output
    sf.load_network(network_file=gpickle_network_file)
    sf.plot_network()
    sf.load_attributes(attribute_file=annotation_file)
    sf.define_neighborhoods(
        node_distance_metric="shortpath_weighted_layout", neighborhood_radius=0.15)
    sf.compute_pvalues()
    sf.define_top_attributes()
    sf.define_domains(attribute_distance_threshold=0.65)
    sf.trim_domains()
    sf.plot_composite_network(show_each_domain=True,
                              save_fig=params.output_file_pdf)
    sf.print_output_files()
    sf.save()

    return None


def define_node_attributes(params, df_saint):
    network_node_annotation = "../results/network_annotations/node.txt"

    safe_node = params.safe_node_attribute
    df_safe_node = pd.read_csv(safe_node, sep="\t")
    del df_safe_node["Unnamed: 0"]
    del df_safe_node["id"]
    del df_safe_node["label"]

    nodes = df_safe_node.key.drop_duplicates().tolist()
    baits = df_saint.BaitGene.drop_duplicates().tolist()

    df_safe_node["node_bait"] = [
        True if key in baits else False for key in df_safe_node.key.tolist()]
    df_safe_node["node_color"] = [NODE_COLORS[x]
                                  for x in df_safe_node.domain.tolist()]

    # To calculate AvgSpecSum from all baits.
    df = df_saint[df_saint.PreyGene.isin(nodes)]
    df = df[["BaitGene", "PreyGene", "AvgSpec"]]
    df = df.groupby(by="PreyGene").sum()
    df['logAvgSpec'] = np.log10(df.AvgSpec)
    df['node_size'] = df['logAvgSpec'] * 20 + 10
    df = df.reset_index(drop=False).rename(columns={"PreyGene": "key"})
    df_node = pd.merge(df_safe_node, df, how='left')

    # Handling minimum node size for CETN2
    min_bait_node = (df_node[df_node.node_bait == True].node_size.min() / 2) + 10
    df_node.iloc[df_node[df_node.key == "CETN2"].index,
                 list(df_node).index("node_size")] = min_bait_node

    df_node.to_csv(network_node_annotation, sep="\t", index=None)
    return df_node


def define_edge_attributes(params):
    edge_file = params.network_edge
    if not os.path.exists(edge_file):
        print("Please export the edge information from the Cytoscape")
        print("File name should be [mTb_network default edge.csv]")
        return None

    df_edge = pd.read_csv(edge_file)

    del df_edge["interaction"]
    del df_edge["shared name"]
    del df_edge["shared interaction"]
    del df_edge["selected"]

    df_edge["edge_size"] = [
        round(math.log2(x) * 3) if x > 0 else 2 for x in df_edge.AvgSpec.tolist()]
    df_edge["edge_color"] = "#C0C0C0"
    df_edge["edge_line_type"] = ["Solid" if x >
                                 0 else "Equal Dash" for x in df_edge.AvgSpec.tolist()]

    df_edge.to_csv("../results/network_annotations/edge.txt",
                   sep="\t", index=None)
    return df_edge
