# KRAS G12C Vignette Parameters

## Input Parameters for BaseBuddy

### Common Parameters (Both Scenarios)
- **Reference Region**: chr12:25245027-25245627 (600bp covering KRAS)
- **Mutation Position**: chr12:25245327 (G>T change)
- **Read Length**: 150bp
- **Insert Size**: 350bp ± 50bp
- **Sequencing Profile**: HS25 (HiSeq 2500)
- **Coverage**: 200x depth
- **Read Type**: Paired-end

### Scenario 1.1: True Positive KRAS G12C (NSCLC)
**Biological Parameters:**
- **VAF**: 45% (high, typical for driver mutation)
- **Strand Bias**: Balanced (50/50 forward/reverse)
- **Base Quality**: High quality for variant bases
- **Read Distribution**: Normal, properly paired

**BaseBuddy Commands:**
```bash
# 1. Simulate normal reads
basebuddy short \
  --reference kras_reference.fa \
  --depth 200 \
  --read-len 150 \
  --profile HS25 \
  --mean-frag-len 350 \
  --std-frag-len 50 \
  --paired

# 2. Spike mutation at 45% VAF
basebuddy spike \
  --reference kras_reference.fa \
  --input-bam normal.bam \
  --snp-vcf kras_g12c_true.vcf \
  --vaf 0.45 \
  --seed 42
```

### Scenario 1.2: Artifactual KRAS G12C
**Artifact Parameters:**
- **VAF**: 1.5% (very low, typical for artifacts)
- **Strand Bias**: Extreme (95% forward strand)
- **Base Quality**: Lower for variant bases
- **Read Distribution**: Soft-clipped reads near variant

**BaseBuddy Commands:**
```bash
# 1. Simulate normal reads (same as above)
basebuddy short \
  --reference kras_reference.fa \
  --depth 200 \
  --read-len 150 \
  --profile HS25 \
  --mean-frag-len 350 \
  --std-frag-len 50 \
  --paired

# 2. Spike mutation at low VAF
basebuddy spike \
  --reference kras_reference.fa \
  --input-bam normal.bam \
  --snp-vcf kras_g12c_artifact.vcf \
  --vaf 0.015 \
  --seed 123

# 3. Introduce strand bias
basebuddy strand-bias \
  --input-bam low_vaf.bam \
  --output-bam artifact_final.bam \
  --forward-fraction 0.95 \
  --seed 456
```

## Key Differences Between Scenarios

| Parameter | True Positive (1.1) | Artifact (1.2) |
|-----------|-------------------|----------------|
| VAF | 45% | 1.5% |
| Strand Bias | Balanced (50/50) | Extreme (95/5) |
| Base Quality | High | Lower for variants |
| Read Pattern | Normal | Soft-clipped |
| Clinical Context | NSCLC tumor | Normal tissue |

## Additional Considerations

### For More Realistic Artifacts:
1. **Position-specific artifacts**: Variants appearing at read ends
2. **Oxidative damage**: C>A changes with strand bias
3. **FFPE artifacts**: C>T changes with specific patterns

### Quality Metrics to Track:
- Strand bias (SB)
- Base quality rank sum
- Mapping quality rank sum
- Read position rank sum
- Clipping rank sum

### VCF Format for Mutations:
```
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr12	25245327	KRAS_G12C	G	T	100	PASS	AF=0.45;SOMATIC
```

## Notes on KRAS Biology:
- KRAS is on the minus strand
- G12C: Glycine to Cysteine at codon 12
- Reference: G on plus strand → T mutation
- Transcript: C>A change (reverse complement)
- Common in NSCLC, rare in other tissues