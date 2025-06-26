# BaseBuddy macOS App Packaging Guide

This guide explains how to package BaseBuddy as a standalone macOS application that users can simply drag and drop to install, without needing Python, Conda, or any other dependencies.

## Overview

We create a self-contained `.app` bundle that includes:
- Python runtime
- All Python packages (pysam, numpy, customtkinter, etc.)
- All bioinformatics tools (BWA, BAMSurgeon, etc.)
- GUI and all resources

The app is built as x86_64 to ensure compatibility with both Intel and Apple Silicon Macs (via Rosetta 2).

## Quick Build

```bash
# One command to build everything
./scripts/build_app.sh
```

This creates:
- `dist/BaseBuddy.app` - The application bundle
- `dist/BaseBuddy.dmg` - Installer disk image

## Detailed Steps

### 1. Prerequisites

- macOS (Intel or Apple Silicon)
- Conda/Miniforge installed
- Xcode Command Line Tools (`xcode-select --install`)

### 2. Create Packaging Environment

The build script automatically creates an x86_64 conda environment with all dependencies:

```bash
export CONDA_SUBDIR=osx-64
conda create -n bb-pack python=3.10 [all packages...] -c conda-forge -c bioconda
```

### 3. Prepare External Tools

Before building, we copy all CLI tools into a `tools/` directory:

```bash
conda activate bb-pack
./scripts/prepare_tools.sh
```

This collects:
- BAMSurgeon tools (addsnv.py, etc.)
- Aligners (BWA, samtools)
- Simulators (ART, wgsim, NanoSim)
- Other dependencies

### 4. Build the App

PyInstaller bundles everything into a macOS app:

```bash
arch -x86_64 pyinstaller basebuddy.spec
```

The spec file:
- Bundles all Python code
- Includes external tools in `Contents/tools/`
- Sets up proper paths
- Creates a native macOS app structure

### 5. Create Installer

The build script automatically creates a DMG:

```bash
hdiutil create -volname "BaseBuddy" -srcfolder dist/BaseBuddy.app -format UDZO dist/BaseBuddy.dmg
```

## Distribution

### For Users

1. Download `BaseBuddy.dmg`
2. Open the DMG
3. Drag BaseBuddy to Applications
4. Double-click to run

No installation of Python, Conda, or command-line tools required!

### Code Signing (Optional)

For distribution without security warnings:

```bash
# Sign the app
codesign --deep --force --sign "Developer ID Application: Your Name" dist/BaseBuddy.app

# Notarize for Gatekeeper
xcrun altool --notarize-app --primary-bundle-id com.basebuddy.gui \
  --username "your@email.com" --password "@keychain:AC_PASSWORD" \
  --file dist/BaseBuddy.dmg
```

## Architecture Details

### Why x86_64?

- Several bioinformatics tools (exonerate, velvet) only have x86_64 builds
- Building for x86_64 ensures all tools work
- Rosetta 2 provides transparent translation on Apple Silicon
- Single build works on all Macs

### Bundle Structure

```
BaseBuddy.app/
├── Contents/
│   ├── Info.plist          # App metadata
│   ├── MacOS/
│   │   └── BaseBuddy       # Main executable
│   ├── Resources/
│   │   ├── basebuddy.icns  # App icon
│   │   └── data/           # Signature files, etc.
│   └── tools/              # All CLI tools
│       ├── addsnv.py
│       ├── bwa
│       ├── samtools
│       └── ...
```

### Path Management

The app automatically:
1. Detects if running as bundled app
2. Adds `tools/` directory to PATH
3. Sets environment variables
4. Ensures all tools are executable

## Customization

### App Icon

Replace `resources/basebuddy.icns` with your custom icon:
- 1024x1024 PNG recommended
- Convert with: `sips -s format icns icon.png --out basebuddy.icns`

### Version Info

Edit `basebuddy.spec` to update:
- Bundle identifier
- Version numbers
- App name
- Copyright info

### Including Additional Files

Add to the `datas` list in `basebuddy.spec`:
```python
datas=[
    ('path/to/file', 'destination/in/bundle'),
    ('path/to/directory/', 'destination/directory/'),
]
```

## Troubleshooting

### "App is damaged" Error

Right-click the app and select "Open" to bypass Gatekeeper on first run.

### Tools Not Found

Check that tools are properly bundled:
```bash
ls -la /Applications/BaseBuddy.app/Contents/tools/
```

### Python Import Errors

Ensure all hidden imports are listed in `basebuddy.spec`:
```python
hiddenimports=['module_name', ...]
```

### Large App Size

The app may be 500MB-1GB due to bundled tools. Consider:
- Excluding unused tools
- Using UPX compression
- Creating separate "lite" version

## Platform Support

### macOS
- Full support via this guide
- Universal app possible with additional work

### Windows
- Use same approach with Windows-specific tools
- Create `.exe` instead of `.app`
- Bundle as installer with NSIS/Inno Setup

### Linux
- Create AppImage or Flatpak
- Bundle tools in similar way
- Distribute as `.tar.gz` or package

## Maintenance

### Updating Tools

1. Update conda environment
2. Re-run `prepare_tools.sh`
3. Rebuild app

### Updating Python Code

Just rebuild - PyInstaller will pick up changes automatically.

### Testing

Always test the packaged app on:
- Clean macOS system
- Both Intel and Apple Silicon
- Different macOS versions (10.13+)