#!/usr/bin/env python3

"""
===========================================================
Feature Extraction Pipeline
Author : Rahul N
Purpose:
    Extract summary statistics from annotated VEP VCF files
    for ASD vs ADHD vs IDDD classification.

Output:
    feature_matrix.csv
===========================================================
"""

import os
import gzip
import logging
from pathlib import Path
from statistics import median
from collections import Counter, defaultdict

import pandas as pd

##############################################################
# PROJECT PATHS
##############################################################

BASE_DIR = Path("/home/ibab/WES_Project/ChrX_Pilot/Annotation")

OUTPUT_DIR = BASE_DIR / "/home/ibab/WES_Project/ChrX_Pilot/output"

OUTPUT_DIR.mkdir(exist_ok=True)

##############################################################
# LOGGING
##############################################################

logging.basicConfig(

    filename=OUTPUT_DIR / "feature_extraction.log",

    level=logging.INFO,

    format="%(asctime)s - %(levelname)s - %(message)s"

)

logging.info("Feature extraction started")

##############################################################
# DATASET CLASSES
##############################################################

CLASSES = {

    "ASD": "ASD",

    "ADHD": "ADHD",

    "IDDD": "IDDD"

}

##############################################################
# OPEN NORMAL OR GZIPPED VCF
##############################################################

def open_vcf(vcf_path):

    if str(vcf_path).endswith(".gz"):

        return gzip.open(vcf_path, "rt")

    return open(vcf_path, "r")

##############################################################
# FIND ALL VCF FILES
##############################################################

def discover_vcfs(base_dir):

    samples = []

    for disease in CLASSES:

        disease_dir = base_dir / disease

        if not disease_dir.exists():

            logging.warning(f"{disease_dir} not found")

            continue

        for root, dirs, files in os.walk(disease_dir):

            for file in files:

                if file.endswith("_Annotated.vcf"):

                    sample_name = file.replace("_Annotated.vcf","")

                    full_path = Path(root) / file

                    samples.append({

                        "sample": sample_name,

                        "label": disease,

                        "vcf": full_path

                    })

    return samples

##############################################################
# CONSEQUENCE SEVERITY RANKING
##############################################################

SEVERITY_RANK = {

    "transcript_ablation": 1,
    "splice_acceptor_variant": 2,
    "splice_donor_variant": 3,
    "stop_gained": 4,
    "frameshift_variant": 5,
    "stop_lost": 6,
    "start_lost": 7,
    "missense_variant": 8,
    "protein_altering_variant": 9,
    "inframe_insertion": 10,
    "inframe_deletion": 11,
    "synonymous_variant": 12,
    "splice_region_variant": 13,
    "5_prime_UTR_variant": 14,
    "3_prime_UTR_variant": 15,
    "intron_variant": 16,
    "upstream_gene_variant": 17,
    "downstream_gene_variant": 18,
    "intergenic_variant": 19

}

##############################################################
# RETURN MOST SEVERE CONSEQUENCE
##############################################################

def get_most_severe_annotation(csq_string, csq_fields):

    best = None
    best_rank = 999

    for record in csq_string.split(","):

        values = record.split("|")

        if len(values) != len(csq_fields):
            continue

        ann = dict(zip(csq_fields, values))

        consequences = ann["Consequence"].split("&")

        for consequence in consequences:

            rank = SEVERITY_RANK.get(consequence, 999)

            if rank < best_rank:

                best_rank = rank
                best = ann.copy()
                best["Consequence"] = consequence

    return best

##############################################################
# EXTRACT SUMMARY FEATURES FROM ONE ANNOTATED VCF
##############################################################

def extract_summary_features(vcf_path):

    features = {}

    ##########################################################
    # Variant Statistics
    ##########################################################

    total_variants = 0
    snv_count = 0
    insertion_count = 0
    deletion_count = 0
    indel_count = 0

    ##########################################################
    # Quality Statistics
    ##########################################################

    qual_scores = []

    pass_count = 0
    filtered_count = 0

    ##########################################################
    # Consequence Statistics
    ##########################################################

    consequence_counter = Counter()

    ##########################################################
    # Impact Statistics
    ##########################################################

    impact_counter = Counter()

    ##########################################################
    # Existing Variant Statistics
    ##########################################################

    known_variants = 0
    novel_variants = 0
    rsid_variants = 0

    ##########################################################
    # CSQ HEADER
    ##########################################################

    csq_fields = None

    ##########################################################
    # READ FILE
    ##########################################################

    with open_vcf(vcf_path) as file:

        for line in file:

            ##################################################
            # Read CSQ Header
            ##################################################

            if line.startswith("##INFO=<ID=CSQ"):

                desc = line.split("Format: ")[1]

                desc = desc.split('">')[0]

                csq_fields = desc.split("|")

                continue

            ##################################################
            # Skip Header
            ##################################################

            if line.startswith("#"):

                continue

            cols = line.strip().split("\t")

            ref = cols[3]
            alt = cols[4]
            qual = cols[5]
            filt = cols[6]
            info = cols[7]

            total_variants += 1

            ##################################################
            # Variant Type
            ##################################################

            if len(ref) == 1 and len(alt) == 1:

                snv_count += 1

            elif len(ref) < len(alt):

                insertion_count += 1

            elif len(ref) > len(alt):

                deletion_count += 1

            else:

                indel_count += 1

            ##################################################
            # QUAL
            ##################################################

            if qual != ".":

                qual_scores.append(float(qual))

            ##################################################
            # FILTER
            ##################################################

            if filt == "." or filt == "PASS":

                pass_count += 1

            else:

                filtered_count += 1

            ##################################################
            # Skip if CSQ missing
            ##################################################

            if "CSQ=" not in info:

                continue

            ##################################################
            # Parse only FIRST annotation
            ##################################################

            csq = info.split("CSQ=")[1].split(";")[0]

            ann = get_most_severe_annotation(csq, csq_fields)

            if ann is None:
                continue

            ##################################################
            # Consequence
            ##################################################

            for c in ann["Consequence"].split("&"):

                consequence_counter[c] += 1

            ##################################################
            # IMPACT
            ##################################################

            impact_counter[ann["IMPACT"]] += 1

            ##################################################
            # Existing Variant
            ##################################################

            existing = ann["Existing_variation"]

            if existing == "":

                novel_variants += 1

            else:

                known_variants += 1

                if "rs" in existing:

                    rsid_variants += 1

    ##########################################################
    # STORE FEATURES
    ##########################################################

    features["Total_variants"] = total_variants

    features["SNV_count"] = snv_count
    features["Insertion_count"] = insertion_count
    features["Deletion_count"] = deletion_count
    features["Indel_count"] = indel_count

    ##########################################################
    # Consequences
    ##########################################################

    consequence_list = [

        "missense_variant",
        "synonymous_variant",
        "frameshift_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "splice_region_variant",
        "intron_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "intergenic_variant"

    ]

    for c in consequence_list:

        features[c] = consequence_counter[c]

    ##########################################################
    # IMPACT
    ##########################################################

    for impact in ["HIGH","MODERATE","LOW","MODIFIER"]:

        features[impact] = impact_counter[impact]

    ##########################################################
    # QUALITY
    ##########################################################

    if len(qual_scores):

        features["Average_QUAL"] = round(sum(qual_scores)/len(qual_scores),2)

        features["Median_QUAL"] = sorted(qual_scores)[len(qual_scores)//2]

    else:

        features["Average_QUAL"] = 0

        features["Median_QUAL"] = 0

    features["PASS_count"] = pass_count

    features["Filtered_count"] = filtered_count

    if total_variants:

        features["PASS_percentage"] = round(

            (pass_count/total_variants)*100,

            2

        )

    else:

        features["PASS_percentage"] = 0

    ##########################################################
    # Existing Variant Statistics
    ##########################################################

    features["Known_variants"] = known_variants

    features["Novel_variants"] = novel_variants

    features["rsID_variants"] = rsid_variants

    return features

##############################################################
# TEST DISCOVERY
##############################################################

##############################################################
# TEST
##############################################################

if __name__ == "__main__":

    samples = discover_vcfs(BASE_DIR)

    sample = samples[0]

    print()

    print("Reading")

    print(sample["sample"])

    print()

    features = extract_summary_features(sample["vcf"])

    print("="*60)

    for k,v in features.items():

        print(f"{k:30} : {v}")