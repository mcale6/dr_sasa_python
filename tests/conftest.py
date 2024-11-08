import os
import sys
import platform
import sysconfig
from pathlib import Path

def get_lib_paths():
    """Dynamically generate possible library paths based on current platform and Python version"""
    root_dir = Path(__file__).parent.parent
    
    # Get Python implementation and version info
    implementation = platform.python_implementation().lower()
    version = sysconfig.get_config_var('py_version_nodot')
    full_version = sysconfig.get_config_var('py_version_nodot_plat')
    
    # Get platform-specific info
    system = platform.system().lower()
    machine = platform.machine().lower()
    
    # Base paths that always should be checked
    paths = [
        root_dir / "build" / "lib",
        root_dir / "build"
    ]
    
    # Platform specific path construction
    if system == 'linux':
        platform_tag = f"linux-{machine}"
    elif system == 'darwin':
        platform_tag = f"macosx-{machine}"
    elif system == 'windows':
        platform_tag = f"win-{machine}"
    else:
        platform_tag = f"{system}-{machine}"
    
    # Add platform-specific paths
    lib_patterns = [
        f"lib.{platform_tag}-{implementation}-{version}",
        f"lib.{platform_tag}-{implementation}{version}",
        f"lib.{platform_tag}-{full_version}",
        f"lib.{platform_tag}",
        f"lib.{system}-{machine}",
    ]
    
    # Add all possible combinations
    for pattern in lib_patterns:
        paths.append(root_dir / "build" / pattern)
    
    # Print debug info if needed
    if os.environ.get('DEBUG'):
        print(f"Python Implementation: {implementation}")
        print(f"Python Version: {version}")
        print(f"Platform: {platform_tag}")
        print("Checking paths:")
        for path in paths:
            print(f"  {path}")
    
    return paths

def pytest_sessionstart(session):
    """Add library paths to sys.path at the start of testing"""
    paths = get_lib_paths()
    
    # Add all existing build paths to sys.path
    added_paths = []
    for build_path in paths:
        if build_path.exists():
            build_path_str = str(build_path.absolute())
            if build_path_str not in sys.path:
                sys.path.insert(0, build_path_str)
                added_paths.append(build_path_str)
    
    if added_paths:
        print("\nAdded build paths:")
        for path in added_paths:
            print(f"  ✓ {path}")
    else:
        print("\nWarning: No valid build paths found!")
        print("Searched in:")
        for path in paths:
            print(f"  ✗ {path}")