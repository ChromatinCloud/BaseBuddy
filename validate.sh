#!/bin/bash

# =================================================================
# BaseBuddy: End-to-End Validation Script
#
# This script runs the validation pipeline. It can be run from the
# beginning or started at a specific step using the --start-at flag.
#
# Usage:
#   ./validate.sh                (Runs all steps from the beginning)
#   ./validate.sh --start-at build
#   ./validate.sh --start-at smoke
#   ./validate.sh --start-at integration
#
# =================================================================

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration ---
IMAGE_TAG="basebuddy-validation:latest"
START_STEP="build" # Default starting step

# --- Argument Parsing ---
# Check if a command-line flag was provided.
if [[ "$1" == "--start-at" ]]; then
  if [[ -n "$2" ]]; then
    START_STEP="$2"
  else
    echo "Error: --start-at flag requires an argument (build, smoke, or integration)."
    exit 1
  fi
fi

# --- Script Start ---
echo "🔵 Starting BaseBuddy End-to-End Validation..."
echo "   Image to be used: $IMAGE_TAG"
echo "   Starting at step: $START_STEP"
echo "-------------------------------------------------"

# STEP 1: Build the Docker Image
if [[ "$START_STEP" == "build" ]]; then
  echo "STEP 1: Building Docker Image..."
  docker build -t "$IMAGE_TAG" .
  echo "✅ Docker image built successfully."
  START_STEP="smoke" # Automatically proceed to the next step
else
  echo "SKIPPING STEP 1: Docker Build"
fi
echo "-------------------------------------------------"


# STEP 2: Run the Smoke Test
if [[ "$START_STEP" == "smoke" ]]; then
  echo "STEP 2: Running Smoke Test..."
  tests/smoke_test.sh
  echo "✅ Smoke test passed."
  START_STEP="integration" # Automatically proceed to the next step
else
  echo "SKIPPING STEP 2: Smoke Test"
fi
echo "-------------------------------------------------"


# STEP 3: Run the Integration Test Suite
if [[ "$START_STEP" == "integration" ]]; then
  echo "STEP 3: Running Integration Test Suite (inside Docker)..."
  docker run --rm "$IMAGE_TAG" pytest tests/test_integration.py
  echo "✅ Integration test suite passed."
  START_STEP="gui" # Automatically proceed to the next step
else
  echo "SKIPPING STEP 3: Integration Tests"
fi
echo "-------------------------------------------------"


# STEP 4: Final Manual GUI Test
if [[ "$START_STEP" == "gui" ]]; then
  echo "🎉 All automated tests passed successfully!"
  echo ""
  echo "The final validation step is to test the GUI manually."
  echo "Please run the following command in your terminal:"
  echo ""
  echo "    python -m src.basebuddy.gui.main_app"
  echo ""
else
  echo "Validation run complete."
fi