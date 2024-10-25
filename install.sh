#!/bin/bash

# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate dr_sasa_env

# Install package
pip install -e .

# Run tests (optional)
pytest tests/