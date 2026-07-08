import os
import subprocess
import time
import pandas as pd
from datetime import datetime

# =============================================================================
# CONFIGURATION & REFS
# =============================================================================
THREADS = "6"  # Allocates 6 cores per step as used in your history

# Base Project Directories
BASE_DIR = os.path.expanduser("~/WES_Project/ChrX_Pilot")
REF_FASTA = os.path.join(BASE_DIR, "reference/GRCh38_primary.fa")
LOG_EXCEL_FILE = os.path.join(BASE_DIR, "logs/pipeline_profiling_report.xlsx")

# Sample Inventory Sheet (Updated with your filenames)
SAMPLES = {
    # =========================
    # ASD Samples
    # =========================
    "ASD_X_01": {
        "cohort": "ASD",
        "r1": os.path.join(BASE_DIR, "raw/ASD/ASD_X_01/SRR1269322_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ASD/ASD_X_01/SRR1269322_2.fastq.gz"),
        "rg": "@RG\\tID:ASD_X_01\\tSM:ASD_X_01\\tPL:ILLUMINA\\tLB:EZ2_WXS\\tPU:SRR1269322"
    },
    "ASD_X_02": {
        "cohort": "ASD",
        "r1": os.path.join(BASE_DIR, "raw/ASD/ASD_X_02/SRR1272237_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ASD/ASD_X_02/SRR1272237_2.fastq.gz"),
        "rg": "@RG\\tID:ASD_X_02\\tSM:ASD_X_02\\tPL:ILLUMINA\\tLB:EZ2_WXS\\tPU:SRR1272237"
    },
    "ASD_X_03": {
        "cohort": "ASD",
        "r1": os.path.join(BASE_DIR, "raw/ASD/ASD_X_03/SRR1272248_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ASD/ASD_X_03/SRR1272248_2.fastq.gz"),
        "rg": "@RG\\tID:ASD_X_03\\tSM:ASD_X_03\\tPL:ILLUMINA\\tLB:EZ2_WXS\\tPU:SRR1272248"
    },
    "ASD_X_04": {
        "cohort": "ASD",
        "r1": os.path.join(BASE_DIR, "raw/ASD/ASD_X_04/SRR1272295_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ASD/ASD_X_04/SRR1272295_2.fastq.gz"),
        "rg": "@RG\\tID:ASD_X_04\\tSM:ASD_X_04\\tPL:ILLUMINA\\tLB:EZ2_WXS\\tPU:SRR1272295"
    },
    "ASD_X_05": {
        "cohort": "ASD",
        "r1": os.path.join(BASE_DIR, "raw/ASD/ASD_X_05/SRR1301297_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ASD/ASD_X_05/SRR1301297_2.fastq.gz"),
        "rg": "@RG\\tID:ASD_X_05\\tSM:ASD_X_05\\tPL:ILLUMINA\\tLB:EZ2_WXS\\tPU:SRR1301297"
    },
    "ASD_X_06": {
        "cohort": "ASD",
        "r1": os.path.join(BASE_DIR, "raw/ASD/ASD_X_06/SRR1301340_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ASD/ASD_X_06/SRR1301340_2.fastq.gz"),
        "rg": "@RG\\tID:ASD_X_06\\tSM:ASD_X_06\\tPL:ILLUMINA\\tLB:EZ2_WXS\\tPU:SRR1301340"
    },
    "ASD_X_07": {
        "cohort": "ASD",
        "r1": os.path.join(BASE_DIR, "raw/ASD/ASD_X_07/SRR1301376_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ASD/ASD_X_07/SRR1301376_2.fastq.gz"),
        "rg": "@RG\\tID:ASD_X_07\\tSM:ASD_X_07\\tPL:ILLUMINA\\tLB:EZ2_WXS\\tPU:SRR1301376"
    },
    "ASD_X_08": {
        "cohort": "ASD",
        "r1": os.path.join(BASE_DIR, "raw/ASD/ASD_X_08/SRR1301424_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ASD/ASD_X_08/SRR1301424_2.fastq.gz"),
        "rg": "@RG\\tID:ASD_X_08\\tSM:ASD_X_08\\tPL:ILLUMINA\\tLB:EZ2_WXS\\tPU:SRR1301424"
    },
    "ASD_X_09": {
        "cohort": "ASD",
        "r1": os.path.join(BASE_DIR, "raw/ASD/ASD_X_09/SRR1301432_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ASD/ASD_X_09/SRR1301432_2.fastq.gz"),
        "rg": "@RG\\tID:ASD_X_09\\tSM:ASD_X_09\\tPL:ILLUMINA\\tLB:EZ2_WXS\\tPU:SRR1301432"
    },
    "ASD_X_10": {
        "cohort": "ASD",
        "r1": os.path.join(BASE_DIR, "raw/ASD/ASD_X_10/SRR1301440_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ASD/ASD_X_10/SRR1301440_2.fastq.gz"),
        "rg": "@RG\\tID:ASD_X_10\\tSM:ASD_X_10\\tPL:ILLUMINA\\tLB:EZ2_WXS\\tPU:SRR1301440"
    },

    # =========================
    # ADHD Samples
    # =========================
    "ADHD_X_01": {
        "cohort": "ADHD",
        "r1": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_01/SRR25993612_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_01/SRR25993612_2.fastq.gz"),
        "rg": "@RG\\tID:ADHD_X_01\\tSM:ADHD_X_01\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR25993612"
    },
    "ADHD_X_02": {
        "cohort": "ADHD",
        "r1": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_02/SRR25993759_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_02/SRR25993759_2.fastq.gz"),
        "rg": "@RG\\tID:ADHD_X_02\\tSM:ADHD_X_02\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR25993759"
    },
    "ADHD_X_03": {
        "cohort": "ADHD",
        "r1": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_03/SRR25993697_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_03/SRR25993697_2.fastq.gz"),
        "rg": "@RG\\tID:ADHD_X_03\\tSM:ADHD_X_03\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR25993697"
    },
    "ADHD_X_04": {
        "cohort": "ADHD",
        "r1": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_04/SRR25993647_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_04/SRR25993647_2.fastq.gz"),
        "rg": "@RG\\tID:ADHD_X_04\\tSM:ADHD_X_04\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR25993647"
    },
    "ADHD_X_05": {
        "cohort": "ADHD",
        "r1": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_05/SRR25993733_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_05/SRR25993733_2.fastq.gz"),
        "rg": "@RG\\tID:ADHD_X_05\\tSM:ADHD_X_05\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR25993733"
    },
    "ADHD_X_06": {
        "cohort": "ADHD",
        "r1": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_06/SRR25993815_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_06/SRR25993815_2.fastq.gz"),
        "rg": "@RG\\tID:ADHD_X_06\\tSM:ADHD_X_06\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR25993815"
    },
    "ADHD_X_07": {
        "cohort": "ADHD",
        "r1": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_07/SRR25993763_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_07/SRR25993763_2.fastq.gz"),
        "rg": "@RG\\tID:ADHD_X_07\\tSM:ADHD_X_07\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR25993763"
    },
    "ADHD_X_08": {
        "cohort": "ADHD",
        "r1": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_08/SRR25993653_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_08/SRR25993653_2.fastq.gz"),
        "rg": "@RG\\tID:ADHD_X_08\\tSM:ADHD_X_08\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR25993653"
    },
    "ADHD_X_09": {
        "cohort": "ADHD",
        "r1": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_09/SRR25993826_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_09/SRR25993826_2.fastq.gz"),
        "rg": "@RG\\tID:ADHD_X_09\\tSM:ADHD_X_09\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR25993826"
    },
    "ADHD_X_10": {
        "cohort": "ADHD",
        "r1": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_10/SRR25993758_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/ADHD/ADHD_X_10/SRR25993758_2.fastq.gz"),
        "rg": "@RG\\tID:ADHD_X_10\\tSM:ADHD_X_10\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR25993758"
    },

    # =========================
    # IDDD Samples
    # =========================
    "IDDD_X_01": {
        "cohort": "IDDD",
        "r1": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_01/SRR27706707_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_01/SRR27706707_2.fastq.gz"),
        "rg": "@RG\\tID:IDDD_X_01\\tSM:IDDD_X_01\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR27706707"
    },
    "IDDD_X_02": {
        "cohort": "IDDD",
        "r1": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_02/SRR27706184_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_02/SRR27706184_2.fastq.gz"),
        "rg": "@RG\\tID:IDDD_X_02\\tSM:IDDD_X_02\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR27706184"
    },
    "IDDD_X_03": {
        "cohort": "IDDD",
        "r1": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_03/SRR27706683_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_03/SRR27706683_2.fastq.gz"),
        "rg": "@RG\\tID:IDDD_X_03\\tSM:IDDD_X_03\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR27706683"
    },
    "IDDD_X_04": {
        "cohort": "IDDD",
        "r1": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_04/SRR27706219_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_04/SRR27706219_2.fastq.gz"),
        "rg": "@RG\\tID:IDDD_X_04\\tSM:IDDD_X_04\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR27706219"
    },
    "IDDD_X_05": {
        "cohort": "IDDD",
        "r1": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_05/SRR27706506_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_05/SRR27706506_2.fastq.gz"),
        "rg": "@RG\\tID:IDDD_X_05\\tSM:IDDD_X_05\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR27706506"
    },
    "IDDD_X_06": {
        "cohort": "IDDD",
        "r1": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_06/SRR27706345_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_06/SRR27706345_2.fastq.gz"),
        "rg": "@RG\\tID:IDDD_X_06\\tSM:IDDD_X_06\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR27706345"
    },
    "IDDD_X_07": {
        "cohort": "IDDD",
        "r1": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_07/SRR27706183_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_07/SRR27706183_2.fastq.gz"),
        "rg": "@RG\\tID:IDDD_X_07\\tSM:IDDD_X_07\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR27706183"
    },
    "IDDD_X_08": {
        "cohort": "IDDD",
        "r1": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_08/SRR27706469_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_08/SRR27706469_2.fastq.gz"),
        "rg": "@RG\\tID:IDDD_X_08\\tSM:IDDD_X_08\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR27706469"
    },
    "IDDD_X_09": {
        "cohort": "IDDD",
        "r1": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_09/SRR27706308_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_09/SRR27706308_2.fastq.gz"),
        "rg": "@RG\\tID:IDDD_X_09\\tSM:IDDD_X_09\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR27706308"
    },
    "IDDD_X_10": {
        "cohort": "IDDD",
        "r1": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_10/SRR27706638_1.fastq.gz"),
        "r2": os.path.join(BASE_DIR, "raw/IDDD/IDDD_X_10/SRR27706638_2.fastq.gz"),
        "rg": "@RG\\tID:IDDD_X_10\\tSM:IDDD_X_10\\tPL:ILLUMINA\\tLB:WXS_chrX\\tPU:SRR27706638"
    }
}

# Ensure log and output folders exist before launching
os.makedirs(os.path.dirname(LOG_EXCEL_FILE), exist_ok=True)


# =============================================================================
# REAL-TIME PROFILE LOG ENGINE
# =============================================================================
def run_and_profile(sample_id, cohort_name, step_name, command, inputs, outputs, use_shell=False):
    """
    Runs shell commands, calculates runtime down to millisecond precision,
    and logs metrics to an Excel file immediately to protect data from crashes.
    """
    print(f"\n[{cohort_name}][{sample_id}] Starting: {step_name}")
    start_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    start_time = time.perf_counter()

    cmd_str = command if use_shell else " ".join(command)
    status = "Success"
    error_msg = "None"

    try:
        # Executes commands safely; captures stderr pipeline logs
        subprocess.run(command,
                       check=True,
                       shell=use_shell,
                       executable="/bin/bash" if use_shell else None,
                       capture_output=True,
                       text=True)
    except subprocess.CalledProcessError as e:
        status = "Failed"
        error_msg = e.stderr if e.stderr else str(e)
        print(f" Error during step {step_name}:\n{error_msg}")

    duration_sec = time.perf_counter() - start_time

    # Structure logging metrics
    log_entry = {
        "Cohort": cohort_name,
        "Sample ID": sample_id,
        "Step Name": step_name,
        "Start Timestamp": start_timestamp,
        "Input File(s)": ", ".join(inputs) if isinstance(inputs, list) else inputs,
        "Output File(s)": ", ".join(outputs) if isinstance(outputs, list) else outputs,
        "Duration (Seconds)": round(duration_sec, 2),
        "Duration (Minutes)": round(duration_sec / 60, 2),
        "Status": status,
        "Command Run": cmd_str,
        "Notes": error_msg[:300]
    }

    # Save data directly to Excel sheet
    df_new = pd.DataFrame([log_entry])
    if os.path.exists(LOG_EXCEL_FILE):
        df_existing = pd.read_excel(LOG_EXCEL_FILE)
        pd.concat([df_existing, df_new], ignore_index=True).to_excel(LOG_EXCEL_FILE, index=False)
    else:
        df_new.to_excel(LOG_EXCEL_FILE, index=False)

    if status == "Failed":
        raise RuntimeError(f"Pipeline stopped: Critical error encountered at step: {step_name}")


# =============================================================================
# WES CORE WORKFLOW PIPELINE
# =============================================================================
def run_wes_pipeline(sample_id, metadata):
    cohort = metadata["cohort"]

    # Verify input files exist before starting heavy computations
    if not os.path.exists(metadata["r1"]) or not os.path.exists(metadata["r2"]):
        print( f" [WARNING] Fastq files missing for {sample_id}. Expected paths:\n R1: {metadata['r1']}\n R2: {metadata['r2']}\n Skipping this sample.")
        return

    # Create required directory targets dynamically
    fastqc_dir = os.path.join(BASE_DIR, f"fastqc/{sample_id}")
    clean_dir = os.path.join(BASE_DIR, f"clean/{cohort}/{sample_id}")
    fastp_dir = os.path.join(BASE_DIR, f"fastp/{sample_id}")
    align_dir = os.path.join(BASE_DIR, f"alignment/{sample_id}")
    chrX_dir = os.path.join(align_dir, "chrX")
    vcf_dir = os.path.join(BASE_DIR, f"vcf/{sample_id}")
    log_dir = os.path.join(BASE_DIR, "logs")

    for d in [fastqc_dir, clean_dir, fastp_dir, align_dir, chrX_dir, vcf_dir, log_dir]:
        os.makedirs(d, exist_ok=True)

    # 1. Raw Quality Control (FastQC)
    run_and_profile(sample_id, cohort, "Raw Quality Control (FastQC)",
                    ["fastqc", "-t", THREADS, "-o", fastqc_dir, metadata["r1"], metadata["r2"]],
                    inputs=[metadata["r1"], metadata["r2"]], outputs=[fastqc_dir])

    # 2. Adapter Trimming and Filtering (fastp)
    trimmed_r1 = os.path.join(clean_dir, f"{sample_id}_1.trimmed.fastq.gz")
    trimmed_r2 = os.path.join(clean_dir, f"{sample_id}_2.trimmed.fastq.gz")
    html_report = os.path.join(fastp_dir, f"{sample_id}_fastp.html")
    json_report = os.path.join(fastp_dir, f"{sample_id}_fastp.json")

    fastp_cmd = [
        "fastp", "-i", metadata["r1"], "-I", metadata["r2"], "-o", trimmed_r1, "-O", trimmed_r2,
        "--detect_adapter_for_pe", "--cut_tail", "--cut_window_size", "4", "--cut_mean_quality", "20",
        "--length_required", "50", "--thread", THREADS, "--html", html_report, "--json", json_report
    ]
    run_and_profile(sample_id, cohort, "Data Quality Trimming (fastp)", fastp_cmd,
                    inputs=[metadata["r1"], metadata["r2"]], outputs=[trimmed_r1, trimmed_r2])

    # 3. Read Mapping & Sorting Pipeline (BWA MEM to Sorted BAM)
    sorted_bam = os.path.join(align_dir, f"{sample_id}.sorted.bam")
    bwa_log = os.path.join(log_dir, f"{sample_id}_bwa.log")

    bwa_mem_pipe = (
        f"set -o pipefail; bwa mem -t {THREADS} -M -R '{metadata['rg']}' {REF_FASTA} "
        f"{trimmed_r1} {trimmed_r2} 2> {bwa_log} | samtools sort -@ {THREADS} -o {sorted_bam}"
    )
    run_and_profile(sample_id, cohort, "Read Mapping (BWA MEM -> samtools sort)", bwa_mem_pipe,
                    inputs=[trimmed_r1, trimmed_r2, REF_FASTA], outputs=[sorted_bam], use_shell=True)

    # 4. SAM/BAM Coordinate Indexing
    run_and_profile(sample_id, cohort, "BAM Coordinate Indexing", ["samtools", "index", "-@", THREADS, sorted_bam],
                    inputs=[sorted_bam], outputs=[f"{sorted_bam}.bai"])

    # 5. Alignment Quality Review (Flagstat)
    flagstat_out = os.path.join(align_dir, f"{sample_id}.flagstat.txt")
    run_and_profile(sample_id, cohort, "Mapping Stats Extraction (Flagstat)",
                    f"samtools flagstat {sorted_bam} > {flagstat_out}",
                    inputs=[sorted_bam], outputs=[flagstat_out], use_shell=True)

    # 6. Name Sorting for Pairing Layout Checks
    namesort_bam = os.path.join(align_dir, f"{sample_id}.namesort.bam")
    run_and_profile(sample_id, cohort, "BAM Name Sorting",
                    ["samtools", "sort", "-n", "-@", THREADS, "-o", namesort_bam, sorted_bam],
                    inputs=[sorted_bam], outputs=[namesort_bam])

    # 7. Mate Correction Setup (Fixmate)
    fixmate_bam = os.path.join(align_dir, f"{sample_id}.fixmate.bam")
    run_and_profile(sample_id, cohort, "Read Mate Correction (Fixmate)",
                    ["samtools", "fixmate", "-m", "-@", THREADS, namesort_bam, fixmate_bam],
                    inputs=[namesort_bam], outputs=[fixmate_bam])

    # 8. Coordinate Re-sorting Post Fixmate
    fixmate_sorted = os.path.join(align_dir, f"{sample_id}.fixmate.sorted.bam")
    run_and_profile(sample_id, cohort, "Fixmate BAM Coordinate Sorting",
                    ["samtools", "sort", "-@", THREADS, "-o", fixmate_sorted, fixmate_bam],
                    inputs=[fixmate_bam], outputs=[fixmate_sorted])

    # 9. Duplicate Read Marking (Markdup)
    markdup_bam = os.path.join(align_dir, f"{sample_id}.markdup.bam")
    markdup_log = os.path.join(log_dir, f"{sample_id}_markdup.log")
    run_and_profile(sample_id, cohort, "PCR Duplicate Marking (Markdup)",
                    f"samtools markdup -@ {THREADS} -s {fixmate_sorted} {markdup_bam} 2> {markdup_log}",
                    inputs=[fixmate_sorted], outputs=[markdup_bam], use_shell=True)

    # 10. Indexed Markdup Lookup Setup
    run_and_profile(sample_id, cohort, "Markdup BAM Indexing", ["samtools", "index", "-@", THREADS, markdup_bam],
                    inputs=[markdup_bam], outputs=[f"{markdup_bam}.bai"])

    # 11. Chromosome X Isolation
    chrX_bam = os.path.join(chrX_dir, f"{sample_id}.X.markdup.bam")
    run_and_profile(sample_id, cohort, "Chromosome X Slicing",
                    ["samtools", "view", "-@", THREADS, "-b", "-F", "0x904", markdup_bam, "X", "-o", chrX_bam],
                    inputs=[markdup_bam], outputs=[chrX_bam])

    # 12. Slice Index Generation
    run_and_profile(sample_id, cohort, "Chromosome X BAM Indexing", ["samtools", "index", "-@", THREADS, chrX_bam],
                    inputs=[chrX_bam], outputs=[f"{chrX_bam}.bai"])

    # 13. Metric Coverage Output Calculation
    coverage_out = os.path.join(chrX_dir, f"{sample_id}.coverage.txt")
    run_and_profile(sample_id, cohort, "ChrX Coverage Evaluation",
                    f"samtools coverage -G 0x400 {chrX_bam} > {coverage_out}",
                    inputs=[chrX_bam], outputs=[coverage_out], use_shell=True)

    # 14. Variant Calling (BCFtools Variant Call to Final VCF Format)
    mpileup_bcf = os.path.join(vcf_dir, f"{sample_id}_raw.bcf")
    final_vcf = os.path.join(vcf_dir, f"{sample_id}_final_variants.vcf")
    run_and_profile(sample_id, cohort, "Mpileup Generation (BCFtools)",
                    ["bcftools", "mpileup", "--threads", THREADS, "-f", REF_FASTA, "-O", "b", "-o", mpileup_bcf, chrX_bam],
                    inputs=[chrX_bam, REF_FASTA], outputs=[mpileup_bcf])
    run_and_profile(sample_id, cohort, "Genomic Variant Calling (BCFtools Call)",
                    ["bcftools", "call", "--threads", THREADS, "-v", "-m", "-O", "v", "-o", final_vcf, mpileup_bcf],
                    inputs=[mpileup_bcf], outputs=[final_vcf])


# =============================================================================
# RUN COORDINATOR
# =============================================================================
if __name__ == "__main__":
    print("Initializing WES Target Chromosome X Automation Run across ASD, ADHD, and IDDD cohorts...")

    for sample, data in SAMPLES.items():
        try:
            print(f"\n===============================\nPROCESSING SAMPLE: {sample}\n===============================")
            run_wes_pipeline(sample, data)
        except Exception as err:
            print(f"Skipping processing chain for sample {sample} due to an unexpected failure: {err}")
            continue

    print(f"\nAll tasks finished. Check your profiling performance report workbook at:\n--> {LOG_EXCEL_FILE}")
