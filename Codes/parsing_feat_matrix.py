#!/usr/bin/env python3

import os
import pandas as pd

# --------------------------------------------------
# Annotation folders
# --------------------------------------------------

BASE_DIR = "/home/ibab/WES_Project/ChrX_Pilot/Annotation"

folders = {
    "ASD": os.path.join(BASE_DIR, "ASD"),
    "ADHD": os.path.join(BASE_DIR, "ADHD"),
    "IDDD": os.path.join(BASE_DIR, "IDDD")
}

# --------------------------------------------------
# Storage
# --------------------------------------------------

sample_genes = {}
sample_labels = {}
all_genes = set()

# --------------------------------------------------
# Parse one annotated VCF
# --------------------------------------------------

def parse_vcf(vcf_file):

    genes = set()
    symbol_index = None

    with open(vcf_file) as f:

        for line in f:

            if line.startswith("##INFO=<ID=CSQ"):

                fields = line.split("Format:")[1]
                fields = fields.replace('">', "").strip().split("|")

                if "SYMBOL" in fields:
                    symbol_index = fields.index("SYMBOL")

            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")

            if len(cols) < 8:
                continue

            info = cols[7]

            if "CSQ=" not in info:
                continue

            csq = info.split("CSQ=")[1].split(";")[0]

            for ann in csq.split(","):

                values = ann.split("|")

                if symbol_index is None:
                    continue

                if len(values) <= symbol_index:
                    continue

                gene = values[symbol_index].strip()

                if gene != "":
                    genes.add(gene)

    return genes

# --------------------------------------------------
# Labels
# --------------------------------------------------

label_map = {
    "ASD": 0,
    "ADHD": 1,
    "IDDD": 2
}

# --------------------------------------------------
# Read all annotated VCFs
# --------------------------------------------------

for disease in folders:

    cohort_folder = folders[disease]

    print(f"\nReading {disease} samples...")

    for sample in sorted(os.listdir(cohort_folder)):

        sample_path = os.path.join(cohort_folder, sample)

        if not os.path.isdir(sample_path):
            continue

        vcf_file = os.path.join(
            sample_path,
            f"{sample}_Annotated.vcf"
        )

        if not os.path.exists(vcf_file):
            print(f"Missing: {vcf_file}")
            continue

        genes = parse_vcf(vcf_file)

        sample_genes[sample] = genes
        sample_labels[sample] = label_map[disease]

        all_genes.update(genes)

# --------------------------------------------------
# Build Feature Matrix
# --------------------------------------------------

gene_list = sorted(all_genes)

rows = []
samples = []

for sample in sorted(sample_genes.keys()):

    row = [1 if gene in sample_genes[sample] else 0
           for gene in gene_list]

    rows.append(row)
    samples.append(sample)

feature_matrix = pd.DataFrame(
    rows,
    columns=gene_list,
    index=samples
)

feature_matrix.index.name = "Sample"

# --------------------------------------------------
# Save Feature Matrix
# --------------------------------------------------

feature_matrix.to_csv("feature_matrix.csv")

# --------------------------------------------------
# Save Labels
# --------------------------------------------------

labels = pd.DataFrame({
    "Sample": samples,
    "Label": [sample_labels[s] for s in samples]
})

labels.to_csv("labels.csv", index=False)

# --------------------------------------------------
# Summary
# --------------------------------------------------

print("\n====================================")
print("Feature matrix saved : feature_matrix.csv")
print("Labels saved         : labels.csv")
print("====================================")

print("Samples processed :", len(samples))
print("Unique genes      :", len(gene_list))
print("Matrix shape      :", feature_matrix.shape)

print("\nFirst five samples:")
print(feature_matrix.head())