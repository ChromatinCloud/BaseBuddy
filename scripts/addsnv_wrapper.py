#!/Users/lauferva/miniconda3/envs/basebuddy/bin/python
"""
Wrapper script that mimics addsnv.py interface for BaseBuddy
Uses pysam to spike variants without needing BAMSurgeon
"""

import argparse
import pysam
import random
import sys
import os
from pathlib import Path

def spike_snv(args):
    """Main function to spike SNV into BAM file"""
    
    # Set random seed if provided
    if args.randseed:
        random.seed(args.randseed)
    
    # Open reference if provided (for validation)
    ref_fasta = None
    if args.reffile:
        ref_fasta = pysam.FastaFile(args.reffile)
    
    # Parse variants - either from string or BED file
    variants = []
    
    if os.path.exists(args.varstring) and os.path.isfile(args.varstring):
        # It's a BED file
        with open(args.varstring, 'r') as bed_file:
            for line in bed_file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 5:
                    # BED format: chrom, start (0-based), end, vaf, alt_base
                    chrom = parts[0]
                    pos = int(parts[1])  # Already 0-based in BED
                    vaf = float(parts[3])
                    alt_base = parts[4]
                    
                    # Try to get ref base from reference if available
                    ref_base = 'N'
                    if ref_fasta and chrom in ref_fasta.references:
                        ref_base = ref_fasta.fetch(chrom, pos, pos + 1).upper()
                    
                    variants.append({
                        'chrom': chrom,
                        'pos': pos,
                        'ref_base': ref_base,
                        'alt_base': alt_base,
                        'vaf': vaf
                    })
    else:
        # It's a variant string
        var_parts = args.varstring.split(',')
        if len(var_parts) >= 4:
            chrom = var_parts[0]
            pos = int(var_parts[1]) - 1  # Convert to 0-based
            ref_base = var_parts[2]
            alt_base = var_parts[3]
            variants.append({
                'chrom': chrom,
                'pos': pos,
                'ref_base': ref_base,
                'alt_base': alt_base,
                'vaf': args.mutfrac
            })
    
    # Open BAM files
    inbam = pysam.AlignmentFile(args.bamfile, "rb")
    outbam = pysam.AlignmentFile(args.outbamfile, "wb", template=inbam)
    
    # Track statistics
    total_reads = 0
    overlapping_reads = {}  # Track per variant
    modified_reads = {}  # Track per variant
    forward_modified = {}  # Track strand-specific modifications
    reverse_modified = {}
    
    for i, var in enumerate(variants):
        overlapping_reads[i] = 0
        modified_reads[i] = 0
        forward_modified[i] = 0
        reverse_modified[i] = 0
    
    # Process all reads
    for read in inbam:
        total_reads += 1
        read_modified = False
        
        # Check each variant
        for i, var in enumerate(variants):
            chrom = var['chrom']
            pos = var['pos']
            ref_base = var['ref_base']
            alt_base = var['alt_base']
            vaf = var['vaf']
            
            # Check if read overlaps this variant position
            if read.reference_name == str(chrom) and read.reference_start <= pos < read.reference_end:
                overlapping_reads[i] += 1
                
                # Apply strand bias
                # For paired reads: R1 strand = is_reverse, R2 effective strand = not is_reverse
                if read.is_paired and read.is_read2:
                    # For R2, the effective strand for variant calling is opposite
                    is_forward_strand = read.is_reverse
                else:
                    # For R1 or single-end reads
                    is_forward_strand = not read.is_reverse
                
                # Decide whether to modify based on strand bias
                strand_pass = False
                if is_forward_strand:
                    strand_pass = random.random() < args.strand_bias
                else:
                    strand_pass = random.random() < (1 - args.strand_bias)
                
                # Decide whether to modify this read based on VAF AND strand bias
                if not read_modified and random.random() < vaf and strand_pass:
                    # Find the position in the read sequence
                    pairs = read.get_aligned_pairs(matches_only=False)
                    
                    for read_pos, ref_pos in pairs:
                        if ref_pos == pos and read_pos is not None:
                            # Get the read sequence as a list
                            seq_list = list(read.query_sequence)
                            qual_list = list(read.query_qualities)
                            
                            # Modify the base
                            current_base = seq_list[read_pos]
                            if ref_base == 'N' or current_base == ref_base:
                                seq_list[read_pos] = alt_base
                                read.query_sequence = ''.join(seq_list)
                                read.query_qualities = qual_list
                                modified_reads[i] += 1
                                if is_forward_strand:
                                    forward_modified[i] += 1
                                else:
                                    reverse_modified[i] += 1
                                read_modified = True
                            break
        
        outbam.write(read)
    
    inbam.close()
    outbam.close()
    
    # Index the output BAM
    pysam.index(args.outbamfile)
    
    # Print statistics
    if not args.quiet:
        print(f"Total reads processed: {total_reads}")
        print(f"Total variants processed: {len(variants)}")
        
        total_overlapping = sum(overlapping_reads.values())
        total_modified = sum(modified_reads.values())
        
        print(f"Total reads overlapping any variant: {total_overlapping}")
        print(f"Total reads modified: {total_modified}")
        
        if args.strand_bias != 0.5:
            print(f"Strand bias applied: {args.strand_bias:.3f} (forward strand fraction)")
        
        if len(variants) <= 10:  # Show per-variant stats for small numbers
            for i, var in enumerate(variants):
                if overlapping_reads[i] > 0:
                    actual_vaf = modified_reads[i] / overlapping_reads[i]
                    fwd_pct = (forward_modified[i] / modified_reads[i] * 100) if modified_reads[i] > 0 else 0
                    rev_pct = (reverse_modified[i] / modified_reads[i] * 100) if modified_reads[i] > 0 else 0
                    print(f"  Variant {var['chrom']}:{var['pos']+1} {var['ref_base']}->{var['alt_base']}: "
                          f"{overlapping_reads[i]} overlapping, {modified_reads[i]} modified "
                          f"(Target VAF: {var['vaf']:.3f}, Actual VAF: {actual_vaf:.3f}, "
                          f"Fwd: {forward_modified[i]} ({fwd_pct:.1f}%), Rev: {reverse_modified[i]} ({rev_pct:.1f}%))")
    
    if ref_fasta:
        ref_fasta.close()

def main():
    parser = argparse.ArgumentParser(description='addsnv.py wrapper - spike SNVs into BAM files')
    
    # Required arguments (matching addsnv.py interface)
    parser.add_argument('-v', '--varstring', required=True,
                        help='Variant string in format: chrom,pos,ref,alt OR path to BED file')
    parser.add_argument('-f', '--bamfile', required=True,
                        help='Input BAM file')
    parser.add_argument('-o', '--outbamfile', required=True,
                        help='Output BAM file')
    parser.add_argument('-r', '--reffile', required=True,
                        help='Reference FASTA file')
    
    # Optional arguments
    parser.add_argument('-p', '--mutfrac', type=float, default=0.1,
                        help='Variant allele frequency (default: 0.1)')
    parser.add_argument('-s', '--randseed', type=int,
                        help='Random seed for reproducibility')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress output messages')
    parser.add_argument('--strand-bias', type=float, default=0.5,
                        help='Fraction of variants on forward strand (0.5 = no bias, 0 = all reverse, 1 = all forward)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.bamfile):
        print(f"Error: Input BAM file not found: {args.bamfile}")
        sys.exit(1)
    
    if not os.path.exists(args.reffile):
        print(f"Error: Reference file not found: {args.reffile}")
        sys.exit(1)
    
    try:
        spike_snv(args)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()