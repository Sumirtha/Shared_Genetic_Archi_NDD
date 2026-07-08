import os
import subprocess
import time
import pandas as pd
from datetime import datetime

# =============================================================================
# PIPELINE CONFIGURATION
# =============================================================================
# Scale threads according to your server capabilities to process heavy data
THREADS = "16"

# Target ENA Accession provided by user
SRR_ACCESSION = "SRR1269322"

# Genomic Reference Configurations (Update paths to match your system)
REF_FASTA = "/home/ibab/WES_Project/ChrX_Pilot/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
STAR_INDEX_DIR = "/path/to/STAR_index_directory"

# Operational Paths
OUTPUT_DIR = "./rna_seq_pipeline_output"
LOG_EXCEL_FILE = "rna_seq_pipeline_profiling_report.xlsx"

os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# CORE METRIC LOGGING ENGINE
# =============================================================================
def run_and_profile_step(step_name, command, input_files, output_files):
    """
    Executes external command strings, monitors elapsed time,
    and appends metrics to an active Excel tracking workbook.
    """
    print(f"\n [LAUNCHING] {step_name}")
    print(f"Executing command: {' '.join(command)}")

    start_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    start_time = time.perf_counter()

    try:
        # Run process synchronously; collect runtime shell errors securely
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        status = "Success"
        notes = "Step completed without errors."
    except subprocess.CalledProcessError as e:
        status = "Failed"
        notes = e.stderr if e.stderr else str(e)
        print(f" [CRITICAL FAILURE] Step '{step_name}' hit an exception!\nDetails: {notes}")

    end_time = time.perf_counter()
    duration_seconds = end_time - start_time
    duration_minutes = duration_seconds / 60

    # Compile performance row
    log_entry = {
        "Step Name": step_name,
        "Start Timestamp": start_timestamp,
        "Input File(s)": ", ".join(input_files) if isinstance(input_files, list) else input_files,
        "Output File(s)": ", ".join(output_files) if isinstance(output_files, list) else output_files,
        "Duration (Seconds)": round(duration_seconds, 2),
        "Duration (Minutes)": round(duration_minutes, 2),
        "Status": status,
        "Command Executed": " ".join(command),
        "Error/Notes Summary": notes[:250]
    }

    # Update spreadsheet instantly to keep data safe from potential system crashes
    df_new = pd.DataFrame([log_entry])
    if os.path.exists(LOG_EXCEL_FILE):
        try:
            df_existing = pd.read_excel(LOG_EXCEL_FILE)
            df_combined = pd.concat([df_existing, df_new], ignore_index=True)
            df_combined.to_excel(LOG_EXCEL_FILE, index=False)
        except Exception as ex:
            print(f"Excel locked or corrupt, routing backup log to CSV: {ex}")
            df_new.to_csv("emergency_pipeline_log.csv", mode='a',
                          header=not os.path.exists("emergency_pipeline_log.csv"), index=False)
    else:
        df_new.to_excel(LOG_EXCEL_FILE, index=False)

    print(f" [COMPLETED] {step_name} | Elapsed Time: {duration_minutes:.2f} mins | Status: {status}")

    if status == "Failed":
        raise RuntimeError(f"Pipeline execution halted due to step failure: {step_name}")


# =============================================================================
# PIPELINE STEP COORDINATOR
# =============================================================================
def main():
    print(f"Starting pipeline sequence for dataset: {SRR_ACCESSION}")
    overall_start = time.perf_counter()

    # Define file path variables
    fastq_1 = os.path.join(OUTPUT_DIR, f"{SRR_ACCESSION}_1.fastq")
    fastq_2 = os.path.join(OUTPUT_DIR, f"{SRR_ACCESSION}_2.fastq")
    bam_prefix = os.path.join(OUTPUT_DIR, "Aligned.out.")
    aligned_bam = f"{bam_prefix}bam"
    sorted_bam = os.path.join(OUTPUT_DIR, "aligned_sorted.bam")
    mpileup_bcf = os.path.join(OUTPUT_DIR, "raw_calls.bcf")
    final_vcf = os.path.join(OUTPUT_DIR, f"{SRR_ACCESSION}_final_variants.vcf")

    # --- STEP 1: Automated Download from ENA/SRA ---
    # Downloads sequence files into the working destination
    download_cmd = [
        "fasterq-dump",
        "--outdir", OUTPUT_DIR,
        "--threads", THREADS,
        "--progress",
        SRR_ACCESSION
    ]
    run_and_profile_step(
        step_name="Data Download (fasterq-dump)",
        command=download_cmd,
        input_files=f"ENA Accession: {SRR_ACCESSION}",
        output_files=[fastq_1, fastq_2]
    )

    # --- STEP 2: Splice-Aware Alignment via STAR ---
    # Maps raw sequencing reads directly to the reference genome layout
    star_cmd = [
        "STAR",
        "--runThreadN", THREADS,
        "--genomeDir", STAR_INDEX_DIR,
        "--readFilesIn", fastq_1, fastq_2,
        "--outFileNamePrefix", bam_prefix,
        "--outSAMtype", "BAM", "Unsorted"
    ]
    run_and_profile_step(
        step_name="Sequence Read Mapping (STAR Alignment)",
        command=star_cmd,
        input_files=[fastq_1, fastq_2],
        output_files=[aligned_bam]
    )

    # --- STEP 3: Coordinate Sorting via Samtools ---
    # Organizes BAM records chronologically by genomic coordinates
    sort_cmd = [
        "samtools", "sort",
        "-@", THREADS,
        "-o", sorted_bam,
        aligned_bam
    ]
    run_and_profile_step(
        step_name="BAM File Sorting (Samtools Sort)",
        command=sort_cmd,
        input_files=[aligned_bam],
        output_files=[sorted_bam]
    )

    # --- STEP 4: Indexing Coordinates via Samtools ---
    # Yields an index lookup table for instant random access across heavy regions
    index_cmd = [
        "samtools", "index",
        "-@", THREADS,
        sorted_bam
    ]
    run_and_profile_step(
        step_name="BAM Index Generation (Samtools Index)",
        command=index_cmd,
        input_files=[sorted_bam],
        output_files=[f"{sorted_bam}.bai"]
    )

    # --- STEP 5: Genomic Mpileup Generation via BCFtools ---
    # Amasses positional match/mismatch ratios across covered target locations
    mpileup_cmd = [
        "bcftools", "mpileup",
        "--threads", THREADS,
        "-f", REF_FASTA,
        "-O", "b",  # Outputs highly compressed BCF streams to conserve space
        "-o", mpileup_bcf,
        sorted_bam
    ]
    run_and_profile_step(
        step_name="Mpileup Mathematical Accumulation (BCFtools Mpileup)",
        command=mpileup_cmd,
        input_files=[sorted_bam, REF_FASTA],
        output_files=[mpileup_bcf]
    )

    # --- STEP 6: Variant Extraction to VCF via BCFtools ---
    # Executes the statistical variant-calling algorithms to determine real mutations
    call_cmd = [
        "bcftools", "call",
        "--threads", THREADS,
        "-v",  # Restricts data to altered variant lines only
        "-m",  # Incorporates multiallelic framework mechanics
        "-O", "v",  # Compiles final text in standardized VCF format
        "-o", final_vcf,
        mpileup_bcf
    ]
    run_and_profile_step(
        step_name="Variant Calling Execution (BCFtools Call)",
        command=call_cmd,
        input_files=[mpileup_bcf],
        output_files=[final_vcf]
    )

    pipeline_duration = (time.perf_counter() - overall_start) / 60
    print(f"\n Success! Whole process concluded in {pipeline_duration:.2f} minutes.")
    print(f"Review your generated performance sheet here: {os.path.abspath(LOG_EXCEL_FILE)}")


if __name__ == "__main__":
    main()
