#!/bin/bash
# Installation script for dr_sasa_python using a virtual environment
# Exit on any error and enable command printing
set -e
set -x

# Detect Python version
echo "Detecting Python version..."
# Get system Python version first
SYS_PYTHON_VERSION=$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
echo "System Python version: $SYS_PYTHON_VERSION"

# Try to use system Python version first, then fall back to Python 3.10
if command -v python$SYS_PYTHON_VERSION &> /dev/null; then
    PYTHON_CMD="python$SYS_PYTHON_VERSION"
    PYTHON_VERSION="$SYS_PYTHON_VERSION"
elif command -v python3.10 &> /dev/null; then
    PYTHON_CMD="python3.10"
    PYTHON_VERSION="3.10"
elif command -v python3 &> /dev/null; then
    PYTHON_CMD="python3"
    PYTHON_VERSION=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
else
    echo "Error: Python 3.x not found"
    exit 1
fi

echo "Using Python version: $PYTHON_VERSION (command: $PYTHON_CMD)"

# Script variables
REPO_PATH="$HOME/dr_sasa_python"
VENV_PATH="$REPO_PATH/.venv"
BUILD_PATH="$REPO_PATH/build"

echo "Starting dr_sasa_python installation..."

# 1. System updates and dependencies
echo "Installing system dependencies..."
# Construct package names based on detected version
PYTHON_PACKAGES="python${PYTHON_VERSION} python${PYTHON_VERSION}-dev python${PYTHON_VERSION}-venv"

# Try to update and install packages, continue if fails (useful for non-Ubuntu systems)
sudo apt-get update || true
sudo apt-get install -y \
    build-essential \
    cmake \
    git \
    ${PYTHON_PACKAGES} \
    python3-full \
    ocl-icd-opencl-dev || true

# 2. Clone or update repository first
echo "Setting up repository..."
if [ ! -d "$REPO_PATH" ]; then
    git clone --recursive https://github.com/mcale6/dr_sasa_python.git "$REPO_PATH"
else
    cd "$REPO_PATH"
    git pull
    git submodule update --init --recursive
fi

# 3. Create and activate virtual environment
echo "Setting up virtual environment..."
if [ ! -d "$VENV_PATH" ]; then
    cd "$REPO_PATH"
    $PYTHON_CMD -m venv .venv
fi
source "$VENV_PATH/bin/activate"

# 4. Install Python package management tools
echo "Installing Python package management tools..."
pip install --upgrade pip setuptools wheel build scikit-build

# 5. Install build dependencies
echo "Installing build dependencies..."
pip install "pybind11[global]" numpy pandas pytest pytest-cov

# 6. Install the package in development mode
echo "Installing dr_sasa_python in development mode..."
cd "$REPO_PATH"
pip install -e .

# 7. Build the project
echo "Building the project..."
mkdir -p "$BUILD_PATH"
cd "$BUILD_PATH"
cmake ..
make -j4

# 8. Set up environment variables
echo "Setting up environment variables..."
# Remove any existing entries
sed -i '/# DR-SASA Python/d' "$HOME/.bashrc"
sed -i '/source.*\.venv\/bin\/activate/d' "$HOME/.bashrc"
sed -i '/export PYTHONPATH.*dr_sasa_python/d' "$HOME/.bashrc"

# Add new environment settings
{
    echo "# DR-SASA Python Virtual Environment and Path"
    echo "source $REPO_PATH/.venv/bin/activate"
    echo "export PYTHONPATH=$BUILD_PATH/lib:$BUILD_PATH:$REPO_PATH:\$PYTHONPATH"
} >> "$HOME/.bashrc"

# 9. Display installation summary
echo "
Installation Summary:
--------------------
Python Version: $PYTHON_VERSION
Python Command: $PYTHON_CMD
Repository: $REPO_PATH
Virtual Environment: $REPO_PATH/.venv
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
