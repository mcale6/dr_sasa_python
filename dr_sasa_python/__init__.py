# dr_sasa_python/__init__.py
"""DR-SASA Python package."""
import os
import sys
from pathlib import Path

# Add compiled module directory to path
_pkg_dir = Path(__file__).parent
_module_dir = _pkg_dir / "bindings" / "python"
if str(_module_dir) not in sys.path:
    sys.path.insert(0, str(_module_dir))

# Import version info
from .version import __version__, __version_info__

# Import the compiled module (dr_sasa_py.so)
import dr_sasa_py

# Make everything available at package level
from dr_sasa_py import *

# Import utilities
from . import utils

# Clean up namespace
del os, sys, Path, _pkg_dir, _module_dir

__all__ = [
    '__version__',
    '__version_info__',
    'SimpleSASA',
    'GenericSASA',
    'DecoupledSASA',
    'AtomStruct',
    'utils'
]