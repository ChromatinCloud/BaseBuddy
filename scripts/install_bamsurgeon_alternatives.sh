#!/bin/bash
# Alternative installation methods for BAMSurgeon on Apple Silicon

echo "BAMSurgeon Installation for Apple Silicon (M1/M2)"
echo "================================================="
echo ""
echo "BAMSurgeon has compatibility issues on ARM64. Here are your options:"
echo ""

echo "Option 1: Use Rosetta 2 (x86_64 emulation)"
echo "-------------------------------------------"
echo "Create an x86_64 conda environment:"
echo ""
echo "  CONDA_SUBDIR=osx-64 conda create -n basebuddy-x86 python=3.9"
echo "  conda activate basebuddy-x86"
echo "  conda config --env --set subdir osx-64"
echo "  conda install -c bioconda bamsurgeon"
echo ""

echo "Option 2: Install from source"
echo "-----------------------------"
echo "Clone and install BAMSurgeon manually:"
echo ""
echo "  git clone https://github.com/adamewing/bamsurgeon.git"
echo "  cd bamsurgeon"
echo "  python setup.py install"
echo ""

echo "Option 3: Use Docker"
echo "--------------------"
echo "Run BAMSurgeon in a Docker container:"
echo ""
echo "  docker pull adamewing/bamsurgeon"
echo "  # Then use docker run to execute addsnv.py"
echo ""

echo "Option 4: Create a custom variant spiking solution"
echo "-------------------------------------------------"
echo "We can create a Python script that performs variant spiking"
echo "without BAMSurgeon, using pysam directly."
echo ""

read -p "Would you like me to create a custom variant spiking script? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
    echo "Creating custom variant spiking script..."
    cat > scripts/spike_variants_custom.py << 'EOF'
#!/usr/bin/env python3
"""
Custom variant spiking implementation for BaseBuddy
Works without BAMSurgeon dependency
"""

import pysam
import random
import argparse
import sys
from pathlib import Path

def spike_variant(bam_path, output_path, chrom, pos, ref, alt, vaf=0.1, seed=None):
    """
    Spike a variant into a BAM file at specified VAF
    """
    if seed is not None:
        random.seed(seed)
    
    # Open input and output BAM files
    inbam = pysam.AlignmentFile(bam_path, "rb")
    outbam = pysam.AlignmentFile(output_path, "wb", template=inbam)
    
    # Convert position to 0-based
    pos = int(pos) - 1
    
    # Track statistics
    total_reads = 0
    overlapping_reads = 0
    modified_reads = 0
    
    for read in inbam:
        total_reads += 1
        
        # Check if read overlaps the variant position
        if read.reference_name == str(chrom) and read.reference_start <= pos < read.reference_end:
            overlapping_reads += 1
            
            # Decide whether to modify this read based on VAF
            if random.random() < float(vaf):
                # Find the position in the read sequence
                read_pos = pos - read.reference_start
                
                # Get the read sequence as a list
                seq_list = list(read.query_sequence)
                
                # Modify the base if it matches the reference
                if 0 <= read_pos < len(seq_list):
                    if seq_list[read_pos] == ref:
                        seq_list[read_pos] = alt
                        read.query_sequence = ''.join(seq_list)
                        modified_reads += 1
        
        outbam.write(read)
    
    inbam.close()
    outbam.close()
    
    # Index the output BAM
    pysam.index(output_path)
    
    print(f"Variant spiking complete:")
    print(f"  Total reads: {total_reads}")
    print(f"  Overlapping reads: {overlapping_reads}")
    print(f"  Modified reads: {modified_reads}")
    print(f"  Actual VAF: {modified_reads/overlapping_reads if overlapping_reads > 0 else 0:.3f}")

def parse_vcf_line(line):
    """Parse a VCF line into variant components"""
    parts = line.strip().split('\t')
    if len(parts) >= 5:
        return {
            'chrom': parts[0],
            'pos': parts[1],
            'ref': parts[3],
            'alt': parts[4]
        }
    return None

def main():
    parser = argparse.ArgumentParser(description='Spike variants into BAM files')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--output', required=True, help='Output BAM file')
    parser.add_argument('--vcf', help='VCF file with variants')
    parser.add_argument('--vaf', type=float, default=0.1, help='Variant allele frequency')
    parser.add_argument('--seed', type=int, help='Random seed')
    
    args = parser.parse_args()
    
    # For now, handle single variant from VCF
    if args.vcf:
        with open(args.vcf, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    variant = parse_vcf_line(line)
                    if variant:
                        spike_variant(
                            args.bam, 
                            args.output,
                            variant['chrom'],
                            variant['pos'],
                            variant['ref'],
                            variant['alt'],
                            args.vaf,
                            args.seed
                        )
                        break  # Only handle first variant for now

if __name__ == '__main__':
    main()
EOF
    
    chmod +x scripts/spike_variants_custom.py
    echo "Custom script created at: scripts/spike_variants_custom.py"
    echo ""
    echo "You can use it like:"
    echo "  python scripts/spike_variants_custom.py --bam input.bam --output output.bam --vcf variants.vcf --vaf 0.3"
fi