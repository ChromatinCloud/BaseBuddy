#!/usr/bin/env python3
"""
Launch BaseBuddy GUI on macOS with conda Python 3.13 tkinter workaround

This script provides a workaround for the tkinter issue with conda Python 3.13 on macOS.
"""

import sys
import os
import subprocess

print("BaseBuddy GUI Launcher for macOS")
print("================================\n")

# Check if we're on macOS
if sys.platform != 'darwin':
    print("This launcher is specifically for macOS. Use run_gui.py directly on other platforms.")
    sys.exit(1)

# Try to detect available Python versions
python_options = []

# Check for homebrew Python
if os.path.exists('/opt/homebrew/bin/python3'):
    python_options.append(('/opt/homebrew/bin/python3', 'Homebrew Python'))
elif os.path.exists('/usr/local/bin/python3'):
    python_options.append(('/usr/local/bin/python3', 'Homebrew Python (Intel Mac)'))

# Check for pyenv
pyenv_root = os.environ.get('PYENV_ROOT', os.path.expanduser('~/.pyenv'))
if os.path.exists(pyenv_root):
    versions_dir = os.path.join(pyenv_root, 'versions')
    if os.path.exists(versions_dir):
        for version in sorted(os.listdir(versions_dir), reverse=True):
            if version.startswith('3.') and not version.startswith('3.13'):
                py_path = os.path.join(versions_dir, version, 'bin', 'python3')
                if os.path.exists(py_path):
                    python_options.append((py_path, f'pyenv Python {version}'))

# Check for Python.org installation
if os.path.exists('/Library/Frameworks/Python.framework/Versions'):
    versions_dir = '/Library/Frameworks/Python.framework/Versions'
    for version in sorted(os.listdir(versions_dir), reverse=True):
        if version.startswith('3.') and not version.startswith('3.13'):
            py_path = os.path.join(versions_dir, version, 'bin', 'python3')
            if os.path.exists(py_path):
                python_options.append((py_path, f'Python.org {version}'))

print("The current conda Python 3.13 has tkinter compatibility issues.")
print("Looking for alternative Python installations...\n")

if not python_options:
    print("No alternative Python installations found.")
    print("\nTo run the GUI, you need to install Python with tkinter support:")
    print("  Option 1: Install Python from python.org")
    print("  Option 2: Install Python via Homebrew: brew install python@3.12")
    print("  Option 3: Use pyenv to install Python 3.12: pyenv install 3.12")
    print("\nAlternatively, you can use the CLI:")
    print("  basebuddy --help")
    sys.exit(1)

print("Found Python installations:")
for i, (path, desc) in enumerate(python_options, 1):
    print(f"  {i}. {desc} - {path}")

# Try each Python until one works
for path, desc in python_options:
    print(f"\nTrying {desc}...")
    
    # First check if it has tkinter
    result = subprocess.run([path, '-c', 'import tkinter'], capture_output=True)
    if result.returncode != 0:
        print(f"  ❌ No tkinter support")
        continue
    
    # Check for customtkinter
    result = subprocess.run([path, '-c', 'import customtkinter'], capture_output=True)
    if result.returncode != 0:
        print(f"  ⚠️  Installing customtkinter...")
        subprocess.run([path, '-m', 'pip', 'install', 'customtkinter', '--user'])
    
    # Check for other dependencies
    for dep in ['pysam', 'typer']:
        result = subprocess.run([path, '-c', f'import {dep}'], capture_output=True)
        if result.returncode != 0:
            print(f"  ⚠️  Installing {dep}...")
            subprocess.run([path, '-m', 'pip', 'install', dep, '--user'])
    
    # Try launching the GUI
    print(f"  ✅ Launching GUI with {desc}...")
    subprocess.run([path, 'run_gui.py'])
    break
else:
    print("\n❌ Could not find a working Python installation with tkinter support.")
    print("\nPlease install Python 3.12 or earlier with tkinter support.")