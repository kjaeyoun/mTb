import pandas as pd
import sys
import requests


def fetch_STRING_db(gene_list=["PCNT", "CDK5RAP2", "PLK1"], score_offset=0.7):
    if len(gene_list) > 2000:
        print("The number of gene is over 2000. Please make sure that the total number of gene is less then 2000")
        sys.exit()
    
    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"
    request_url = "/".join([string_api_url, output_format, method])

    my_genes = gene_list
    params = {
        "identifiers": "%0d".join(my_genes),  # your protein
        "species": 9606,  # species NCBI identifier - Currently Human
        "network_flavot": "confidence",
        "caller_identity": "TEDC2 analysis",  # your app name
    }
    response = requests.post(request_url, data=params)

    output_dict = {}
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")
        # print(l)
        p1, p2 = l[2], l[3]
        interaction = [p1, p2]
        interaction.sort()
        interaction = "{}::{}".format(interaction[0], interaction[1])
        # filter the interaction according to experimental score
        string_score = float(l[5])
        if string_score > score_offset:
            if interaction not in output_dict:
                output_dict[interaction] = string_score
            elif output_dict[interaction] < string_score:
                output_dict[interaction] = string_score
            # print("\t".join([p1, p2, "STRING score (prob. %.3f)" % string_score]))

    # To make a dataframe from the output results
    protein_a = [x.split("::")[0] for x in output_dict.keys()]
    protein_b = [x.split("::")[1] for x in output_dict.keys()]
    string_score = [x for x in output_dict.values()]

    string_df = pd.DataFrame(
        {"Protein_A": protein_a, "Protein_B": protein_b, "STRING score": string_score}
    )

    return string_df  # l



# To make an cytoscsape input file to draw a network using the BioID and STRING reuslts
def merge_STRING_into_BioID(df_pool, df_string, output_file):
    # To merge the BioID and STRING results
    df_pool = df_pool[["BaitGene", "PreyGene", "AvgSpec", "BFDR"]]
    df_pool = df_pool.rename(
        columns={"BaitGene": "Protein_A", "PreyGene": "Protein_B"})
    output_df = pd.concat([df_pool, df_string])

    # To change NaN (empty) values to -1
    output_df = output_df.fillna(-1)
    output_df.to_csv(output_file, sep="\t", index=None)
    
    return output_df




if __name__ == "__main__":
    print('string.py')
