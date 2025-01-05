#!/bin/bash
# Installation script for dr_sasa_python

# Exit on any error and enable command printing
set -e
set -x  # Print commands being executed

# Script variables
REPO_PATH="$HOME/dr_sasa_python"
BUILD_PATH="$REPO_PATH/build"

echo "Starting dr_sasa_python installation..."

# 1. System updates and dependencies
echo "Installing system dependencies..."
sudo apt-get update && sudo apt-get install -y \
    build-essential \
    cmake \
    git \
    python3 \
    python3-dev \
    ocl-icd-opencl-dev

# 2. Install Python package management tools
echo "Installing Python package management tools..."
pip install --upgrade pip setuptools wheel build scikit-build

# 3. Install build dependencies
echo "Installing build dependencies..."
pip install "pybind11[global]" numpy pandas pytest pytest-cov

# 4. Clone or update repository
if [ ! -d "$REPO_PATH" ]; then
    echo "Cloning dr_sasa_python repository..."
    git clone --recursive https://github.com/mcale6/dr_sasa_python.git "$REPO_PATH"
else
    echo "Updating existing repository..."
    cd "$REPO_PATH"
    git pull
    git submodule update --init --recursive
fi

# 5. Install the package in development mode
echo "Installing dr_sasa_python in development mode..."
cd "$REPO_PATH"
pip install -e .

# 6. Build the project
echo "Building the project..."
mkdir -p "$BUILD_PATH"
cd "$BUILD_PATH"
cmake ..
make -j$(nproc)  # Use all available CPU cores

# 7. Set up environment variables
echo "Setting up environment variables..."
# Remove any existing dr_sasa entries from .bashrc
sed -i '/# DR-SASA Python/d' "$HOME/.bashrc"
sed -i '/export PYTHONPATH.*dr_sasa_python/d' "$HOME/.bashrc"

# Add new environment settings
{
    echo "# DR-SASA Python Path"
    echo "export PYTHONPATH=\$PYTHONPATH:$BUILD_PATH/lib:$BUILD_PATH:$REPO_PATH"
} >> "$HOME/.bashrc"

# 8. Display installation summary
echo "
Installation Summary:
--------------------
Repository: $REPO_PATH
Build Directory: $BUILD_PATH
Python Path: ${PYTHONPATH}

To start using dr_sasa_python:
1. Close and reopen your terminal, or run: source ~/.bashrc
2. Run 'pytest -v' to verify the installation
"

# Optional: Run tests if specified
if [ "$1" = "--test" ]; then
    echo "Running tests..."
    pytest -v
fi