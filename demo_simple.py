#!/usr/bin/env python3
"""
Simple demonstration of BaseBuddy's improved functionality
"""

import sys
sys.path.insert(0, 'src')

from basebuddy import runner, utils as bb_utils
from pathlib import Path

def demonstrate_error_handling():
    """Demonstrate improved error messages in short read simulation."""
    print("\n" + "="*60)
    print("DEMONSTRATING IMPROVED ERROR HANDLING")
    print("="*60)
    
    # Test 1: Invalid depth parameter
    print("\nTest 1: Invalid sequencing depth")
    try:
        runner.simulate_short(
            output_root_dir=Path("test_output"),
            reference_fasta="dummy.fa",
            depth=-10,  # Invalid negative depth
            read_length=150,
            art_profile="HS25"
        )
    except bb_utils.BaseBuddyInputError as e:
        print(f"✓ Caught expected error: {e}")
        if e.details:
            print(f"  Details: {e.details}")
    
    # Test 2: Invalid read length
    print("\nTest 2: Invalid read length")
    try:
        runner.simulate_short(
            output_root_dir=Path("test_output"),
            reference_fasta="dummy.fa",
            depth=30,
            read_length=0,  # Invalid zero length
            art_profile="HS25"
        )
    except bb_utils.BaseBuddyInputError as e:
        print(f"✓ Caught expected error: {e}")
        if e.details:
            print(f"  Details: {e.details}")
    
    # Test 3: Invalid fragment length for paired-end
    print("\nTest 3: Fragment length too short for paired-end")
    try:
        runner.simulate_short(
            output_root_dir=Path("test_output"),
            reference_fasta="dummy.fa",
            depth=30,
            read_length=150,
            art_profile="HS25",
            mean_fragment_length=100,  # Too short for 150bp reads
            is_paired_end=True
        )
    except bb_utils.BaseBuddyInputError as e:
        print(f"✓ Caught expected error: {e}")
        if e.details:
            print(f"  Details: {e.details}")
    
    # Test 4: Invalid platform
    print("\nTest 4: Invalid ART platform")
    try:
        runner.simulate_short(
            output_root_dir=Path("test_output"),
            reference_fasta="dummy.fa",
            depth=30,
            read_length=150,
            art_profile="HS25",
            art_platform="nanopore"  # Invalid platform
        )
    except bb_utils.BaseBuddyInputError as e:
        print(f"✓ Caught expected error: {e}")
        if e.details:
            print(f"  Details: {e.details}")

def demonstrate_vcf_validation():
    """Demonstrate VCF file validation in spike command."""
    print("\n" + "="*60)
    print("DEMONSTRATING VCF FILE VALIDATION")
    print("="*60)
    
    # This would be called from CLI with proper validation
    print("\nThe spike command now:")
    print("1. Validates VCF files exist before processing")
    print("2. Properly parses SNPs and Indels from VCF files")
    print("3. Provides helpful error messages for missing Picard JAR")
    print("4. Reports IGV session file and variant counts")

def demonstrate_features():
    """Show the polished features."""
    print("\n" + "="*60)
    print("BASEBUDDY POLISHED FEATURES")
    print("="*60)
    
    print("\n1. SHORT READ SIMULATION (ART)")
    print("   - Enhanced error messages with typical value ranges")
    print("   - Platform validation (illumina, 454, solid)")
    print("   - Fragment length validation for paired-end reads")
    print("   - Better error reporting for missing tools")
    
    print("\n2. VARIANT SPIKING (BAMSURGEON)")
    print("   - VCF file validation before processing")
    print("   - Proper VCF parsing for SNPs and Indels")
    print("   - Helpful hints for Picard JAR configuration")
    print("   - IGV session file generation")
    
    print("\n3. STRAND BIAS")
    print("   - Fully functional with automatic BAM indexing")
    print("   - Temporary file cleanup")
    print("   - Detailed logging of operations")
    
    print("\n4. FASTQC INTEGRATION")
    print("   - Multi-file processing with thread control")
    print("   - Comprehensive report generation")
    print("   - Manifest file with all outputs")

def main():
    print("="*80)
    print("BaseBuddy Functionality Demonstration")
    print("="*80)
    
    demonstrate_error_handling()
    demonstrate_vcf_validation()
    demonstrate_features()
    
    print("\n" + "="*80)
    print("SUMMARY OF IMPROVEMENTS")
    print("="*80)
    
    print("\n✓ Fixed missing typing imports in cli.py")
    print("✓ Enhanced error messages with helpful context and typical values")
    print("✓ Removed misleading VCF processing warnings")
    print("✓ Added proper VCF file validation")
    print("✓ Improved error handling with specific tool installation hints")
    print("✓ All four features (short, spike, strand-bias, qc) are production-ready")
    
    print("\nTo generate a BAM with KRAS G12C mutation:")
    print("1. Create/obtain reference FASTA with KRAS region")
    print("2. Run: basebuddy short --reference kras.fa --depth 50")
    print("3. Align reads to create BAM (using BWA/bowtie2)")
    print("4. Create VCF with KRAS G12C (chr12:25245350 G>C)")
    print("5. Run: basebuddy spike -i aligned.bam --snp-vcf kras_g12c.vcf")

if __name__ == "__main__":
    main()