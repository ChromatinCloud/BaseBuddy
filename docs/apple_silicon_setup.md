# BaseBuddy on Apple Silicon: Hybrid Environment Setup

This guide addresses the architecture compatibility challenges when running BaseBuddy on Apple Silicon (M1/M2/M3) Macs.

## The Challenge

Several critical bioinformatics tools in the BaseBuddy pipeline are only available as x86_64 binaries in Bioconda:
- **exonerate** - No ARM64 build
- **velvet** - No ARM64 build  
- **wgsim** - No ARM64 build
- **nanosim** - No ARM64 build (due to htseq/genometools dependencies)
- **BAMSurgeon** - Depends on above tools

## The Solution: Hybrid Environment

We use a two-environment approach:
1. **ARM64 environment** (`basebuddy-arm`) - For native performance with compatible tools
2. **x86_64 environment** (`basebuddy-x86`) - For legacy tools via Rosetta 2

## Quick Setup

```bash
# Run the automated setup
./scripts/setup_hybrid_env.sh

# Activate native ARM environment (default)
source activate_basebuddy.sh

# Or activate x86 environment for legacy tools
source activate_basebuddy.sh x86
```

## Environment Contents

### ARM64 Environment (`basebuddy-arm`)
Native ARM64 builds for optimal performance:
- ✅ Python 3.10+
- ✅ BWA
- ✅ samtools, bcftools
- ✅ minimap2
- ✅ ART (Illumina simulator)
- ✅ All Python packages (pysam, numpy, pandas, etc.)
- ✅ Picard (Java-based)

### x86_64 Environment (`basebuddy-x86`)
Legacy tools running under Rosetta 2:
- ⚠️ BAMSurgeon (addsnv.py, addindel.py)
- ⚠️ exonerate
- ⚠️ velvet
- ⚠️ wgsim
- ⚠️ NanoSim
- ⚠️ LAST aligner
- ⚠️ htseq
- ⚠️ genometools

## Usage Patterns

### Pattern 1: Manual Environment Switching

```bash
# For short read simulation (native ARM64)
source activate_basebuddy.sh
python run_gui.py  # or basebuddy simulate-short-reads

# For variant spiking (needs BAMSurgeon)
source activate_basebuddy.sh x86
python run_gui.py  # Use variant spiking tab
```

### Pattern 2: Smart Tool Wrapper

The setup creates a `basebuddy-tool` wrapper that automatically switches environments:

```bash
# Automatically uses ARM64 env
./basebuddy-tool bwa mem ref.fa reads.fq > output.sam

# Automatically switches to x86 env
./basebuddy-tool addsnv.py -v chr1,100,A,T -f input.bam -o output.bam
```

### Pattern 3: Homebrew Alternatives

Some tools have ARM64 builds via Homebrew:
```bash
# Install native alternatives
brew install exonerate
brew install blast

# Then use system version instead of conda version
/opt/homebrew/bin/exonerate --version
```

## GUI Usage

The BaseBuddy GUI works in both environments, but functionality varies:

| Tab | ARM64 Environment | x86_64 Environment |
|-----|-------------------|-------------------|
| Short Read Sim | ✅ Full support | ✅ Full support |
| Long Read Sim | ❌ No NanoSim | ✅ Full support |
| Variant Spiking | ❌ No BAMSurgeon | ✅ Full support |
| Apply Signature | ✅ Full support | ✅ Full support |
| Germline Sim | ✅ Full support* | ✅ Full support |

*SimuG is installed separately in both environments

## Performance Considerations

### Native ARM64 (basebuddy-arm)
- **2-3x faster** for compute-intensive tasks
- Lower memory usage
- Better battery life
- Recommended for most workflows

### Rosetta 2 (basebuddy-x86)
- **20-40% slower** than native
- Higher memory overhead
- Required for certain tools
- Transparent operation (no user intervention needed)

## Troubleshooting

### "Tool not found" Errors

1. Check which environment is active:
   ```bash
   echo $CONDA_DEFAULT_ENV
   ```

2. Verify tool installation:
   ```bash
   conda list | grep <tool_name>
   ```

3. Check architecture:
   ```bash
   file $(which <tool_name>)
   ```

### Rosetta 2 Installation

If you get Rosetta errors:
```bash
softwareupdate --install-rosetta --agree-to-license
```

### Environment Activation Issues

If `source activate_basebuddy.sh` fails:
```bash
# Manually activate
conda activate basebuddy-arm  # or basebuddy-x86

# Add SimuG to PATH if needed
export PATH="$PWD/external/simuG:$PATH"
```

## Advanced: Custom Tool Installation

### Building from Source

For tools without ARM64 binaries, you can compile from source:

```bash
# Example: Building exonerate from source
git clone https://github.com/nathanweeks/exonerate.git
cd exonerate
./configure --prefix=$CONDA_PREFIX
make && make install
```

### Using Docker

Alternative approach using x86_64 containers:
```bash
# Run x86_64 tools in Docker
docker run --platform linux/amd64 -v $PWD:/data biocontainers/bamsurgeon:latest \
  addsnv.py -v chr1,100,A,T -f /data/input.bam -o /data/output.bam
```

## Recommendations by Workflow

### Primarily Short Read Analysis
- Use ARM64 environment exclusively
- Install any missing tools via Homebrew

### Mixed Short/Long Read Analysis
- Keep both environments
- Use the activation script to switch as needed

### Heavy Variant Manipulation
- Consider keeping a Linux VM or Docker setup
- Or use the x86 environment despite performance penalty

## Future Outlook

As more bioinformatics tools get ARM64 builds, the need for x86_64 environments will diminish. Check periodically for updates:

```bash
# Check for ARM64 builds
conda search -c bioconda <package_name> --platform osx-arm64
```

Tools likely to get ARM64 support soon:
- minimap2 (already available)
- STAR aligner (in progress)
- More to come as Apple Silicon adoption increases