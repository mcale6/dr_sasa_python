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
    #python3-dev \
    #python3-venv \
    #python3-full \
    #ocl-icd-opencl-dev

# 2. Create and activate virtual environment
echo "Setting up virtual environment..."
python3 -m venv ~/dr_sasa_venv
source ~/dr_sasa_venv/bin/activate

# 3. Install Python package management tools in virtual environment
echo "Installing Python package management tools..."
pip install --upgrade pip setuptools wheel build scikit-build

# 4. Install build dependencies
echo "Installing build dependencies..."
pip install "pybind11[global]" numpy pandas pytest pytest-cov

# 5. Clone the repository (if not already in it)
if [ ! -d "dr_sasa_python" ] && [ ! -f "setup.py" ]; then
    echo "Cloning dr_sasa_python repository..."
    git clone --recursive https://github.com/mcale6/dr_sasa_python.git
    cd dr_sasa_python
fi

# 6. Install the package in development mode
echo "Installing dr_sasa_python in development mode..."
pip install -e ".[test]"

# 7. Run tests to verify installation
echo "Running tests..."
pytest -v

# 8. Set up virtual environment activation (optional)
read -p "Would you like to automatically activate the virtual environment in your bashrc? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
    echo "# DR-SASA Python Virtual Environment" >> ~/.bashrc
    echo "source ~/dr_sasa_venv/bin/activate" >> ~/.bashrc
    echo "Virtual environment activation added to ~/.bashrc"
fi

echo "Installation complete!"
echo "You can now import dr_sasa_python in Python"
echo "Remember to activate the virtual environment with: source ~/dr_sasa_venv/activate"