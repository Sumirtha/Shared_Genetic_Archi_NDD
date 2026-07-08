#!/usr/bin/env python3

import os
import subprocess

# --------------------------------------------------
# Paths
# --------------------------------------------------

BASE_DIR = "/home/ibab/WES_Project/ChrX_Pilot"

VCF_DIR = os.path.join(BASE_DIR, "vcf")
ANNOTATION_DIR = os.path.join(BASE_DIR, "Annotation")

VEP = "/home/ibab/ensembl-vep/vep"

# --------------------------------------------------
# Samples
# --------------------------------------------------

SAMPLES = {
    "ASD": [f"ASD_X_{i:02d}" for i in range(1, 11)],
    "ADHD": [f"ADHD_X_{i:02d}" for i in range(1, 11)],
    "IDDD": [f"IDDD_X_{i:02d}" for i in range(1, 11)]
}

# --------------------------------------------------
# Annotate
# --------------------------------------------------

for cohort in SAMPLES:

    print(f"\n========== {cohort} ==========")

    for sample in SAMPLES[cohort]:

        input_vcf = os.path.join(
            VCF_DIR,
            sample,
            "filtered",
            "qual_dp",
            f"{sample}_filtered.vcf.gz"
        )

        if not os.path.exists(input_vcf):
            print(f"Skipping {sample} (Filtered VCF not found)")
            continue

        # Create sample folder
        sample_folder = os.path.join(
            ANNOTATION_DIR,
            cohort,
            sample
        )

        os.makedirs(sample_folder, exist_ok=True)

        output_vcf = os.path.join(
            sample_folder,
            f"{sample}_Annotated.vcf"
        )

        stats_file = os.path.join(
            sample_folder,
            f"{sample}_Annotated.vcf_summary.html"
        )

        print(f"Annotating {sample}...")

        cmd = [
            "perl",
            VEP,

            "-i", input_vcf,
            "-o", output_vcf,

            "--vcf",
            "--cache",
            "--offline",
            "--assembly", "GRCh38",

            "--stats_file", stats_file,

            "--fork", "4",

            "--force_overwrite"
        ]

        try:
            subprocess.run(cmd, check=True)
            print(f"{sample} Done")

        except subprocess.CalledProcessError:
            print(f"{sample} Failed")

print("\n================================")
print(" All samples annotated successfully ")
print("================================")