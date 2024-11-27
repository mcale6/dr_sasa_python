#!/bin/bash
# Installation script for dr_sasa_python using a virtual environment
set -e  # Exit on any error

echo "Starting dr_sasa_python installation..."

# 1. System updates and dependencies (if on a fresh VM)
echo "Installing system dependencies..."
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    git \
    python3 \
    python3-dev \
    python3-venv \
    python3-full \
    ocl-icd-opencl-dev \
    #nvidia-cuda-toolkit \
    #gcc-offload-nvptx \
    #libgomp-plugin-nvptx1 \
    #nvidia-cuda-dev \
    #gcc-11-offload-nvptx \
    #libgomp1

# 2. Create and activate virtual environment
echo "Setting up virtual environment..."
python3 -m venv ~/dr_sasa_venv
source ~/dr_sasa_venv/bin/activate

# 3. Install Python package management tools in virtual environment
echo "Installing Python package management tools..."
pip install --upgrade pip setuptools wheel

# 4. Install Python dependencies from pyproject.toml not needed
echo "Installing Python dependencies..."
pip install "pybind11[global]" numpy  pandas pytest

# 5. Clone and build the repository
echo "Cloning and building dr_sasa_python..."
git clone --recursive https://github.com/mcale6/dr_sasa_python.git
cd dr_sasa_python

# 6. Build the project
echo "Building the project..."
mkdir -p build
cd build
cmake ..
make -j4

# 7. Set up Python path and virtual environment activation
echo "Setting up environment..."
INSTALL_PATH=$(pwd)/build/lib
echo "export PYTHONPATH=\$PYTHONPATH:$INSTALL_PATH" >> ~/.bashrc

# Add virtual environment activation to bashrc 
echo "# DR-SASA Python Virtual Environment" >> ~/.bashrc
echo "source ~/dr_sasa_venv/bin/activate" >> ~/.bashrc

# 8 Test 
pytest -v

