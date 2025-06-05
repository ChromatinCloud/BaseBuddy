# File: README.md

# BaseBuddy

**BaseBuddy** is a command-line toolkit for simulating sequencing reads, spiking in variants, generating mutational signatures, and introducing strand bias—all in one package. You can install and run it locally in a Python virtual environment, or build and run it inside Docker.

---

## Table of Contents

1. [Installation (Local)](#installation-local)  
2. [Usage (Local)](#usage-local)  
   - [Simulating Short Reads (ART)](#simulating-short-reads-art)  
   - [Simulating Long Reads (NanoSim-h)](#simulating-long-reads-nanosim-h)  
   - [Spiking Variants (BAMSurgeon)](#spiking-variants-bamsurgeon)  
   - [Simulating Signatures (SigProfilerSimulator)](#simulating-signatures-sigprofilersimulator)  
   - [Introducing Strand Bias (samtools)](#introducing-strand-bias-samtools)  
3. [Docker Setup & Usage](#docker-setup--usage)  
4. [Running Tests](#running-tests)  
5. [Examples](#examples)  
6. [License](#license)

---

## Installation (Local)

1. **Clone the repo** and create a virtual environment:
   ```bash
   git clone https://github.com/yourusername/BaseBuddy.git
   cd BaseBuddy
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip

