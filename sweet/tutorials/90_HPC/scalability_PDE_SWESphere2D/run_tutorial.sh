#! /bin/bash

# Stop on first error
set -e

# First, cleanup things
mule.benchmark.cleanup_all

# Create all job directories
./1_benchmark_create_jobs.py

# Compile stuff
./compile_platform_*.sh

# Run benchmarks
mule.benchmark.jobs_run_directly job_*

# Plot nice looking pictures
./2_postprocessing_plot.py
