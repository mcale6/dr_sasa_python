#!/bin/bash
# Installation script for dr_sasa_python using a virtual environment

# Exit on any error and enable command printing
set -e
set -x  # Print commands being executed

# Script variables
VENV_PATH="$HOME/dr_sasa_venv"
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
    python3-venv \
    python3-full \
    ocl-icd-opencl-dev

# 2. Create and activate virtual environment
echo "Setting up virtual environment..."
if [ ! -d "$VENV_PATH" ]; then
    python3 -m venv "$VENV_PATH"

fi
source "$VENV_PATH/bin/activate"

# 3. Install Python package management tools
echo "Installing Python package management tools..."
pip install --upgrade pip setuptools wheel build scikit-build

# 4. Install build dependencies
echo "Installing build dependencies..."
pip install "pybind11[global]" numpy pandas pytest pytest-cov

# 5. Clone or update repository
if [ ! -d "$REPO_PATH" ]; then
    echo "Cloning dr_sasa_python repository..."
    git clone --recursive https://github.com/mcale6/dr_sasa_python.git "$REPO_PATH"

else
    echo "Updating existing repository..."
    cd "$REPO_PATH"
    git pull
    git submodule update --init --recursive

fi

# 6. Install the package in development mode
echo "Installing dr_sasa_python in development mode..."
cd "$REPO_PATH"
pip install -e .

# 7. Build the project
echo "Building the project..."
mkdir -p "$BUILD_PATH"
cd "$BUILD_PATH"
cmake ..
make -j$(nproc)  # Use all available CPU cores

# 8. Set up environment variables
echo "Setting up environment variables..."
# Remove any existing dr_sasa entries from .bashrc
sed -i '/# DR-SASA Python/d' "$HOME/.bashrc"
sed -i '/source.*dr_sasa_venv\/bin\/activate/d' "$HOME/.bashrc"
sed -i '/export PYTHONPATH.*dr_sasa_python/d' "$HOME/.bashrc"

# Add new environment settings
{
    echo "# DR-SASA Python Virtual Environment and Path"
    echo "source $VENV_PATH/bin/activate"
    echo "export PYTHONPATH=\$PYTHONPATH:$BUILD_PATH/lib:$BUILD_PATH:$REPO_PATH"
} >> "$HOME/.bashrc"

# 9. Display installation summary
echo "
Installation Summary:
--------------------
Virtual Environment: $VENV_PATH
Repository: $REPO_PATH
Build Directory: $BUILD_PATH
Python Path: $BUILD_PATH/lib:$BUILD_PATH:$REPO_PATH

To start using dr_sasa_python:
1. Close and reopen your terminal, or run: source ~/.bashrc
2. The virtual environment will be activated automatically
3. Run 'pytest -v' to verify the installation
"

# Optional: Run tests if specified
if [ "$1" = "--test" ]; then
    echo "Running tests..."
    pytest -v
fi