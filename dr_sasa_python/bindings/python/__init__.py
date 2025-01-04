"""Python bindings for DR-SASA."""
import os
import sys

# Add current directory to Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# Import the compiled module
from .dr_sasa_py import *

# Clean up namespace
del os, sys, current_dir