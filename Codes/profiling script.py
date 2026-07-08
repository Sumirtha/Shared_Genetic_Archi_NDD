import os
import subprocess
import time
import pandas as pd
from datetime import datetime

# =============================================================================
# CONFIGURATION & PARAMETERS
# =============================================================================
# Threads to speed up computationally heavy processes
THREADS = "16"

# Reference Genome Paths
REF_FASTA = "/path/to/reference/genome.fa"
STAR_INDEX_DIR = "/path/to/STAR_index"

# Input Raw Data
FASTQ_1 = "/path/to/data/sample_R1.fastq.gz"
FASTQ_2 = "/path/to/data/sample_R2.fastq.gz"

# Output Destination Directories
OUTPUT_DIR = "./pipeline_outputs"
LOG_EXCEL_FILE = "rna_seq_pipeline_profiling.xlsx"

os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# PROFILING HELPER FUNCTION
# =============================================================================
def run_and_profile_step(step_name, command, input_files, output_files):
    """
    Executes a shell command, measures runtime, and logs performance to Excel.
    """
    print(f"\n[STARTING] Step: {step_name}")
    print(f"Running command: {' '.join(command)}")

    start_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    start_time = time.perf_counter()

    try:
        # Run the heavy command and wait for completion
        # shell=False is preferred for security and process tracking
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        status = "Success"
        error_msg = "None"
    except subprocess.CalledProcessError as e:
        status = "Failed"
        error_msg = e.stderr if e.stderr else str(e)
        print(f"[ERROR] Step '{step_name}' failed!\n{error_msg}")

    end_time = time.perf_counter()
    duration_seconds = end_time - start_time
    duration_minutes = duration_seconds / 60

    # Structure data row for profiling
    log_entry = {
        "Step Name": step_name,
        "Start Timestamp": start_timestamp,
        "Input File(s)": ", ".join(input_files) if isinstance(input_files, list) else input_files,
        "Output File(s)": ", ".join(output_files) if isinstance(output_files, list) else output_files,
        "Duration (Seconds)": round(duration_seconds, 2),
        "Duration (Minutes)": round(duration_minutes, 2),
        "Status": status,
        "Command Executed": " ".join(command),
        "Error/Notes": error_msg[:200]  # Limit message length in Excel
    }

    # Save directly to Excel sheet (appends or creates new)
    df_new = pd.DataFrame([log_entry])
    if os.path.exists(LOG_EXCEL_FILE):
        try:
            df_existing = pd.read_excel(LOG_EXCEL_FILE)
            df_combined = pd.concat([df_existing, df_new], ignore_index=True)
            df_combined.to_excel(LOG_EXCEL_FILE, index=False)
        except Exception as e:
            print(f"Could not update Excel sheet: {e}. Writing to emergency backup CSV.")
            df_new.to_csv("backup_log.csv", mode='a', header=not os.path.exists("backup_log.csv"), index=False)
    else:
        df_new.to_excel(LOG_EXCEL_FILE, index=False)

    print(f"[COMPLETED] {step_name} in {duration_minutes:.2f} minutes. Status: {status}")

    if status == "Failed":
        raise RuntimeError(f"Pipeline stopped due to critical failure at step: {step_name}")


# =============================================================================
# PIPELINE EXECUTION
# =============================================================================
def main():
    print("Initializing RNA-Seq Variant Calling Pipeline...")
    pipeline_start = time.perf_counter()

    # Define intermediate filenames
    bam_prefix = os.path.join(OUTPUT_DIR, "Aligned.out.")
    aligned_bam = f"{bam_prefix}bam"  # STAR output if using direct BAM generation
    sorted_bam = os.path.join(OUTPUT_DIR, "aligned_sorted.bam")
    mpileup_bcf = os.path.join(OUTPUT_DIR, "raw_variants.bcf")
    final_vcf = os.path.join(OUTPUT_DIR, "final_variants.vcf")

    # --- STEP 1: Alignment using STAR ---
    # Converts FASTQ files directly into an unsorted BAM file
    star_cmd = [
        "STAR",
        "--runThreadN", THREADS,
        "--genomeDir", STAR_INDEX_DIR,
        "--readFilesIn", FASTQ_1, FASTQ_2,
        "--readFilesCommand", "zcat",  # Assumes .gz compressed fastq inputs
        "--outFileNamePrefix", bam_prefix,
        "--outSAMtype", "BAM", "Unsorted"
    ]
    run_and_profile_step(
        step_name="Alignment (STAR Fastq to SAM/BAM)",
        command=star_cmd,
        input_files=[FASTQ_1, FASTQ_2],
        output_files=[aligned_bam]
    )

    # --- STEP 2: Coordinate Sorting using Samtools ---
    # Sorts the aligned reads by genomic location for variant calling compatibility
    sort_cmd = [
        "samtools", "sort",
        "-@", THREADS,
        "-o", sorted_bam,
        aligned_bam
    ]
    run_and_profile_step(
        step_name="BAM Coordinate Sorting (Samtools)",
        command=sort_cmd,
        input_files=[aligned_bam],
        output_files=[sorted_bam]
    )

    # --- STEP 3: BAM Indexing using Samtools ---
    # Generates a .bai index file required for quick random access by downstream tools
    index_cmd = [
        "samtools", "index",
        "-@", THREADS,
        sorted_bam
    ]
    run_and_profile_step(
        step_name="BAM Indexing (Samtools)",
        command=index_cmd,
        input_files=[sorted_bam],
        output_files=[f"{sorted_bam}.bai"]
    )

    # --- STEP 4: Mpileup Generation using BCFtools ---
    # Collects summary information across the genome from the alignments
    mpileup_cmd = [
        "bcftools", "mpileup",
        "--threads", THREADS,
        "-f", REF_FASTA,
        "-O", "b",  # Output compressed BCF format for processing efficiency
        "-o", mpileup_bcf,
        sorted_bam
    ]
    run_and_profile_step(
        step_name="Genomic Mpileup Generation (BCFtools)",
        command=mpileup_cmd,
        input_files=[sorted_bam, REF_FASTA],
        output_files=[mpileup_bcf]
    )

    # --- STEP 5: Variant Calling to VCF using BCFtools ---
    # Calls actual SNVs and indels, outputting the final readable text format VCF
    call_cmd = [
        "bcftools", "call",
        "--threads", THREADS,
        "-v",  # Output variant sites only
        "-m",  # Use the multiallelic caller algorithm
        "-O", "v",  # Output format uncompressed VCF
        "-o", final_vcf,
        mpileup_bcf
    ]
    run_and_profile_step(
        step_name="Variant Calling (BCFtools Call to VCF)",
        command=call_cmd,
        input_files=[mpileup_bcf],
        output_files=[final_vcf]
    )

    total_duration = (time.perf_counter() - pipeline_start) / 60
    print(f"\n✅ Pipeline completed successfully in {total_duration:.2f} minutes!")
    print(f"Profile logging completed in: {os.path.abspath(LOG_EXCEL_FILE)}")


if __name__ == "__main__":
    main()
