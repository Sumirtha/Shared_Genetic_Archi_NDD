#!/usr/bin/env python3

import os
import subprocess

BASE_DIR = "/home/ibab/WES_Project/ChrX_Pilot/vcf"

samples = {
    "ASD": [f"ASD_X_{i:02d}" for i in range(1,11)],
    "ADHD": [f"ADHD_X_{i:02d}" for i in range(1,11)],
    "IDDD": [f"IDDD_X_{i:02d}" for i in range(1,11)]
}

for cohort in samples:

    print(f"\n========== {cohort} ==========")

    for sample in samples[cohort]:

        sample_dir = os.path.join(BASE_DIR, sample)

        input_vcf = os.path.join(
            sample_dir,
            f"{sample}_final_variants.vcf"
        )

        output_dir = os.path.join(
            sample_dir,
            "filtered",
            "qual_dp"
        )

        os.makedirs(output_dir, exist_ok=True)

        output_vcf = os.path.join(
            output_dir,
            f"{sample}_filtered.vcf.gz"
        )

        cmd = [
            "bcftools",
            "view",
            "-i",
            "QUAL>=20 && INFO/DP>=10",
            "-Oz",
            "-o",
            output_vcf,
            input_vcf
        ]

        subprocess.run(cmd, check=True)

        subprocess.run([
            "tabix",
            "-p",
            "vcf",
            output_vcf
        ], check=True)

        print(f"{sample} Filtered")

print("\nAll 30 samples filtered successfully.")