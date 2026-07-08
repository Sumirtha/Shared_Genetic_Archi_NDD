#!/usr/bin/env python3

import os
import pandas as pd
# ----------------------------
# Folder locations
# ----------------------------

folders = {
    "ASD": "/home/ibab/WES_Project/ChrX_Pilot/Annotation/ASD",
    "ADHD": "/home/ibab/WES_Project/ChrX_Pilot/Annotation/ADHD",
    "IDDD": "/home/ibab/WES_Project/ChrX_Pilot/Annotation/IDDD"
}

# Store data
sample_genes = {}
sample_labels = {}
all_genes = set()

# ----------------------------
# Read one VCF file
# ----------------------------

def parse_vcf(vcf_file):

    genes = set()
    symbol_index = None

    with open(vcf_file, "r") as file:

        for line in file:

            # Find the SYMBOL column in the CSQ annotation
            if line.startswith("##INFO=<ID=CSQ"):

                format_text = line.split("Format:")[1]
                format_text = format_text.replace('">', "").strip()

                fields = format_text.split("|")

                if "SYMBOL" in fields:
                    symbol_index = fields.index("SYMBOL")

            # Skip header lines
            if line.startswith("#"):
                continue

            columns = line.strip().split("\t")

            if len(columns) < 8:
                continue

            info = columns[7]

            if "CSQ=" not in info:
                continue

            csq = info.split("CSQ=")[1].split(";")[0]

            annotations = csq.split(",")

            for ann in annotations:

                values = ann.split("|")

                if symbol_index is None:
                    continue

                if len(values) <= symbol_index:
                    continue

                gene = values[symbol_index].strip()

                if gene != "":
                    genes.add(gene)

    return genes

# ----------------------------
# Read all VCF files
# ----------------------------

label_map = {
    "ASD": 0,
    "ADHD": 1,
    "IDDD": 2
}

for disease in folders:

    folder = folders[disease]

    print(f"\nReading {disease} files...")

    for file in os.listdir(folder):

        if file.endswith(".vcf"):

            sample = file.replace(".vcf", "")

            filepath = os.path.join(folder, file)

            genes = parse_vcf(filepath)

            sample_genes[sample] = genes
            sample_labels[sample] = label_map[disease]

            all_genes.update(genes)

print("\nTotal unique genes:", len(all_genes))

# ----------------------------
# Build Feature Matrix
# ----------------------------

gene_list = sorted(list(all_genes))

data = []

samples = []

for sample in sample_genes:

    row = []

    for gene in gene_list:

        if gene in sample_genes[sample]:
            row.append(1)
        else:
            row.append(0)

    data.append(row)
    samples.append(sample)

feature_matrix = pd.DataFrame(
    data,
    columns=gene_list,
    index=samples
)

feature_matrix.index.name = "Sample"

# ----------------------------
# Save Feature Matrix
# ----------------------------

feature_matrix.to_csv("feature_matrix.csv")

print("\nFeature matrix saved as feature_matrix.csv")

# ----------------------------
# Save Labels
# ----------------------------

labels = pd.DataFrame({
    "Sample": samples,
    "Label": [sample_labels[s] for s in samples]
})

labels.to_csv("labels.csv", index=False)

print("Labels saved as labels.csv")

print("\nFeature Matrix Shape:", feature_matrix.shape)

print(feature_matrix.head())