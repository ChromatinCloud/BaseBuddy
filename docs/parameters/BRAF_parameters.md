# BRAF V600K MNV Vignette Parameters

## Input Parameters for BaseBuddy

### Common Parameters (Both Scenarios)
- **Reference Region**: chr7:140753300-140924900 (BRAF gene)
- **Target Region**: chr7:140753336-140753337 (codon 600, positions 1-2)
- **Mutation**: V600K (Valine → Lysine)
- **Nucleotide Change**: GT → AA (dinucleotide substitution)
- **Read Length**: 150bp
- **Insert Size**: 300bp ± 40bp
- **Sequencing Profile**: HS25 (HiSeq 2500)
- **Coverage**: 150x depth
- **Read Type**: Paired-end

### Scenario 2.1: Correct MNV Representation
**Biological Parameters:**
- **Variant Type**: Multi-nucleotide variant (MNV/dinucleotide)
- **Change**: GT>AA at chr7:140753336-140753337
- **VAF**: 30% (typical for melanoma driver)
- **Phasing**: Both mutations on same DNA molecule
- **Clinical Context**: UV-induced melanoma mutation

**BaseBuddy Commands:**
```bash
# 1. Simulate normal reads
basebuddy short \
  --reference braf_reference.fa \
  --depth 150 \
  --read-len 150 \
  --profile HS25 \
  --mean-frag-len 300 \
  --std-frag-len 40 \
  --paired

# 2. Spike MNV (both changes together)
basebuddy spike \
  --reference braf_reference.fa \
  --input-bam normal.bam \
  --snp-vcf braf_v600k_mnv.vcf \
  --vaf 0.30 \
  --seed 42
```

**VCF Format for MNV:**
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr7	140753336	BRAF_V600K_MNV	GT	AA	100	PASS	AF=0.30;TYPE=MNV
```

### Scenario 2.2: Incorrect Separate SNV Representation
**Mis-calling Parameters:**
- **Variant Type**: Two separate SNVs
- **Change 1**: G>A at chr7:140753336
- **Change 2**: T>A at chr7:140753337
- **VAF**: 30% each (but on different molecules)
- **Phasing**: Mutations on DIFFERENT DNA molecules
- **Issue**: Misses functional MNV, incorrect interpretation

**BaseBuddy Workflow:**
```bash
# 1. Split reads into two subsets
samtools view -h -s 0.5 normal.bam > subset1.bam
samtools view -h -s 0.5 normal.bam > subset2.bam

# 2. Spike first SNV into subset1
basebuddy spike \
  --reference braf_reference.fa \
  --input-bam subset1.bam \
  --snp-vcf braf_snv1.vcf \
  --vaf 0.60 \
  --seed 123

# 3. Spike second SNV into subset2
basebuddy spike \
  --reference braf_reference.fa \
  --input-bam subset2.bam \
  --snp-vcf braf_snv2.vcf \
  --vaf 0.60 \
  --seed 456

# 4. Merge BAMs
samtools merge final.bam snv1.bam snv2.bam
```

## Key Differences Between Scenarios

| Feature | Correct (MNV) | Incorrect (2 SNVs) |
|---------|--------------|-------------------|
| Variant Type | Single MNV | Two separate SNVs |
| VCF Entries | 1 entry (GT>AA) | 2 entries (G>A, T>A) |
| Read Pattern | Both changes on same read | Changes on different reads |
| Phasing | Cis (linked) | Trans (unlinked) |
| Functional Impact | Correctly identifies V600K | May miss V600K |
| IGV Display | Reads show "AA" | Reads show "A-" or "-A" |

## Biological Significance

### BRAF V600K Biology:
- **Gene**: BRAF (B-Raf proto-oncogene)
- **Location**: Chromosome 7, minus strand
- **Codon 600**: GTG (Valine) → AAG (Lysine)
- **Function**: Activating mutation in MAPK pathway
- **Cancer Type**: ~20% of melanomas, some colorectal

### UV Signature:
- Dinucleotide substitutions are characteristic of UV damage
- CC>TT is the classic UV signature
- GT>AA (complement of CA>TT) also UV-related

## IGV Visualization Guide

### What to Look For:

**Scenario 2.1 (Correct MNV):**
- Reads spanning both positions show "AA"
- No reads with partial changes ("AT" or "GA")
- Clean substitution pattern

**Scenario 2.2 (Incorrect SNVs):**
- Some reads show "AT" (only first change)
- Other reads show "GA" (only second change)
- No reads with both changes "AA"
- Demonstrates phasing problem

## Clinical Implications

1. **Correct MNV calling**:
   - Identifies V600K mutation
   - Enables targeted therapy (BRAF inhibitors)
   - Accurate prognostic information

2. **Incorrect separate SNV calling**:
   - May miss V600K if not properly phased
   - Could lead to missed therapeutic opportunity
   - Highlights importance of MNV-aware variant callers

## Technical Notes

### MNV Detection Requirements:
1. Variant caller must support MNV detection
2. Sufficient read length to span both positions
3. Proper read phasing algorithms
4. Annotation tools that recognize MNVs

### Recommended Variant Callers:
- Mutect2 (GATK) - MNV-aware
- Strelka2 - Can detect MNVs
- VarDict - Good MNV support
- FreeBayes - Haplotype-based, handles MNVs

### Quality Metrics:
- Check phasing quality scores
- Verify read support for complete MNV
- Look for strand bias (UV signature)
- Confirm both changes in cis