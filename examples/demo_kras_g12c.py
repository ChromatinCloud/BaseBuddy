#!/usr/bin/env python3
"""
KRAS G12C Demonstration Script for BaseBuddy

This script demonstrates BaseBuddy's capabilities by:
1. Creating a small reference containing the KRAS gene region
2. Simulating short reads using ART
3. Spiking in a KRAS G12C mutation (chr12:25245350 G>C in GRCh38)
4. Running QC on the generated reads
5. Applying strand bias to demonstrate that feature

The KRAS G12C mutation is a common oncogenic mutation where the glycine at 
codon 12 is replaced by cysteine due to a G>C substitution.
"""

import subprocess
import sys
from pathlib import Path
import tempfile

def run_command(cmd, description):
    """Run a command and handle errors."""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Command: {' '.join(cmd)}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("SUCCESS!")
        if result.stdout:
            print(f"Output:\n{result.stdout}")
        return result
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Command failed with exit code {e.returncode}")
        if e.stdout:
            print(f"Stdout:\n{e.stdout}")
        if e.stderr:
            print(f"Stderr:\n{e.stderr}")
        sys.exit(1)

def create_kras_reference():
    """Create a small reference FASTA containing KRAS gene region."""
    # This is a simplified KRAS gene region sequence (not actual genomic sequence)
    # The actual KRAS G12C mutation occurs at chr12:25245350 in GRCh38
    # For demo purposes, we'll create a simple sequence with the mutation site
    
    kras_fasta = """
>chr12_KRAS_region
ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAA
TTCAGAATCATTTTGTGGACGAATATGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGA
TGGAGAAACCTGTCTCTTGGATATTCTCGACACAGCAGGTCAAGAGGAGTACAGTGCAATGAGGGACCAG
TACATGAGGACTGGGGAGGGCTTTCTTTGTGTATTTGCCATAAATAATACTAAATCATTTGAAGATATTC
ACCATTATAGAGAACAAATTAAAAGAGTTAAGGACTCTGAAGATGTACCTATGGTCCTAGTAGGAAATAA
ATGTGATTTGCCTTCTAGAACAGTAGACACAAAACAGGCTCAGGACTTAGCAAGAAGTTATGGAATTCCT
TTTATTGAAACATCAGCAAAGACAAGACAGGGTGTTGATGATGCCTTCTATACATTAGTTCGAGAAATTC
""".strip()
    
    ref_file = Path("kras_reference.fasta")
    with open(ref_file, "w") as f:
        f.write(kras_fasta)
    
    print(f"Created reference FASTA: {ref_file}")
    return ref_file

def create_kras_g12c_vcf():
    """Create a VCF file with the KRAS G12C mutation."""
    # In the simplified reference above, we'll place a G>C mutation
    # Position 35 contains a G that we'll change to C for demonstration
    
    vcf_content = """##fileformat=VCFv4.2
##reference=kras_reference.fasta
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##contig=<ID=chr12_KRAS_region,length=490>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr12_KRAS_region\t35\tKRAS_G12C\tG\tC\t100\tPASS\tAF=0.5
"""
    
    vcf_file = Path("kras_g12c.vcf")
    with open(vcf_file, "w") as f:
        f.write(vcf_content.strip())
    
    print(f"Created VCF file: {vcf_file}")
    return vcf_file

def main():
    """Run the KRAS G12C demonstration."""
    print("="*80)
    print("BaseBuddy KRAS G12C Demonstration")
    print("="*80)
    
    # Check if basebuddy is available
    try:
        subprocess.run(["basebuddy", "--help"], check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ERROR: basebuddy command not found. Please ensure BaseBuddy is installed.")
        print("Try: pip install -e .")
        sys.exit(1)
    
    # Create demo directory
    demo_dir = Path("demo_kras_g12c_output")
    demo_dir.mkdir(exist_ok=True)
    print(f"\nDemo output directory: {demo_dir}")
    
    # Step 1: Create reference FASTA
    print("\n" + "="*60)
    print("STEP 1: Creating KRAS reference sequence")
    print("="*60)
    ref_fasta = create_kras_reference()
    
    # Step 2: Simulate short reads
    print("\n" + "="*60)
    print("STEP 2: Simulating Illumina short reads")
    print("="*60)
    
    short_reads_dir = demo_dir / "short_reads"
    run_command([
        "basebuddy", "short",
        "--reference", str(ref_fasta),
        "--outdir", str(short_reads_dir),
        "--depth", "50",
        "--read-len", "150",
        "--profile", "HS25",
        "--mean-frag-len", "300",
        "--std-frag-len", "30",
        "--paired"
    ], "Simulating paired-end Illumina reads at 50x depth")
    
    # Find the generated FASTQ files
    fastq_files = list(short_reads_dir.rglob("*.fq"))
    if len(fastq_files) < 2:
        print("ERROR: Expected paired-end FASTQ files not found")
        sys.exit(1)
    
    print(f"\nGenerated FASTQ files:")
    for fq in fastq_files:
        print(f"  - {fq}")
    
    # Step 3: Create BAM from simulated reads
    print("\n" + "="*60)
    print("STEP 3: Creating BAM file from simulated reads")
    print("="*60)
    
    # This would normally be done with a proper aligner like BWA or bowtie2
    # For demo purposes, we'll use samtools to create a simple BAM
    print("Note: In a real workflow, you would align reads with BWA or similar")
    
    # Create a dummy BAM file for spiking demonstration
    dummy_sam = demo_dir / "dummy.sam"
    with open(dummy_sam, "w") as f:
        # SAM header
        f.write("@HD\tVN:1.0\tSO:coordinate\n")
        f.write("@SQ\tSN:chr12_KRAS_region\tLN:490\n")
        # Add some dummy aligned reads around position 35
        for i in range(100):
            f.write(f"read{i}\t0\tchr12_KRAS_region\t{20+i%30}\t60\t150M\t*\t0\t0\t{'A'*150}\t{'I'*150}\n")
    
    dummy_bam = demo_dir / "dummy.bam"
    run_command([
        "samtools", "view", "-b", str(dummy_sam), "-o", str(dummy_bam)
    ], "Converting SAM to BAM")
    
    run_command([
        "samtools", "sort", str(dummy_bam), "-o", str(dummy_bam.with_stem("dummy_sorted"))
    ], "Sorting BAM")
    
    sorted_bam = dummy_bam.with_stem("dummy_sorted")
    
    run_command([
        "samtools", "index", str(sorted_bam)
    ], "Indexing BAM")
    
    # Step 4: Spike in KRAS G12C mutation
    print("\n" + "="*60)
    print("STEP 4: Spiking KRAS G12C mutation into BAM")
    print("="*60)
    
    kras_vcf = create_kras_g12c_vcf()
    spiked_output_dir = demo_dir / "spiked_bams"
    
    # Note: This requires picard.jar to be available
    print("\nNote: BAMSurgeon requires picard.jar. Setting BAMSURGEON_PICARD_JAR if available...")
    
    run_command([
        "basebuddy", "spike",
        "--reference", str(ref_fasta),
        "--input-bam", str(sorted_bam),
        "--snp-vcf", str(kras_vcf),
        "--output-prefix", str(spiked_output_dir / "kras_g12c_spiked"),
        "--vaf", "0.3",  # 30% variant allele frequency
        "--seed", "42"
    ], "Spiking KRAS G12C mutation at 30% VAF")
    
    # Step 5: Run QC on original FASTQ files
    print("\n" + "="*60)
    print("STEP 5: Running FastQC on simulated reads")
    print("="*60)
    
    qc_output_dir = demo_dir / "fastqc_results"
    fastq_args = [str(f) for f in fastq_files]
    
    run_command([
        "basebuddy", "qc",
        *fastq_args,
        "--output-dir", str(qc_output_dir),
        "--threads", "2"
    ], "Running FastQC quality control")
    
    # Step 6: Demonstrate strand bias
    print("\n" + "="*60)
    print("STEP 6: Introducing strand bias")
    print("="*60)
    
    biased_bam = demo_dir / "strand_biased.bam"
    
    run_command([
        "basebuddy", "strand-bias",
        "--input-bam", str(sorted_bam),
        "--output-bam", str(biased_bam),
        "--forward-fraction", "0.8",  # Keep 80% of forward strand reads
        "--seed", "123"
    ], "Creating strand-biased BAM (80% forward strand)")
    
    # Summary
    print("\n" + "="*80)
    print("DEMONSTRATION COMPLETE!")
    print("="*80)
    print(f"\nAll outputs are in: {demo_dir.absolute()}")
    print("\nKey outputs:")
    print(f"  - Reference FASTA: {ref_fasta}")
    print(f"  - Simulated reads: {short_reads_dir}")
    print(f"  - KRAS G12C VCF: {kras_vcf}")
    print(f"  - Spiked BAMs: {spiked_output_dir}")
    print(f"  - FastQC results: {qc_output_dir}")
    print(f"  - Strand-biased BAM: {biased_bam}")
    
    print("\nThis demonstration showed:")
    print("  1. Short read simulation with ART")
    print("  2. Variant spiking with BAMSurgeon (KRAS G12C)")
    print("  3. Quality control with FastQC")
    print("  4. Strand bias introduction")
    
    print("\nNote: The KRAS G12C mutation (chr12:25245350 G>C) is a common")
    print("oncogenic mutation found in various cancers, particularly lung cancer.")

if __name__ == "__main__":
    main()