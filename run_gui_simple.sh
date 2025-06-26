#!/bin/bash
# Simple launcher when already in basebuddy environment

# Ensure the basebuddy environment's bin directory is in PATH
export PATH="/Users/lauferva/miniconda3/envs/basebuddy/bin:$PATH"

# Use the specific Python from the basebuddy environment
/Users/lauferva/miniconda3/envs/basebuddy/bin/python run_gui.py