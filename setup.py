from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os
import subprocess
import platform

class CMakeExtension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        project_root = os.path.dirname(os.path.abspath(__file__))
        
        # Create build temp directory
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
            
        # Set architecture flags for macOS
        if platform.system() == "Darwin":
            # Check if we're on Apple Silicon
            if platform.machine() == 'arm64':
                arch_flag = "-DCMAKE_OSX_ARCHITECTURES=arm64"
            else:
                arch_flag = "-DCMAKE_OSX_ARCHITECTURES=x86_64"
        else:
            arch_flag = ""
            
        # Configure cmake arguments
        cmake_args = [
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={self.build_temp}',
            f'-DPYTHON_EXECUTABLE={sys.executable}',
            f'-DCMAKE_BUILD_TYPE=Debug'
        ]
        
        if arch_flag:
            cmake_args.append(arch_flag)
        
        print(f"Project root: {project_root}")
        print(f"Build temp: {self.build_temp}")
        print(f"CMake args: {cmake_args}")
        print(f"Machine architecture: {platform.machine()}")
        
        try:
            # Run cmake
            print("\nRunning CMake...")
            subprocess.check_call(
                ['cmake', project_root] + cmake_args,
                cwd=self.build_temp
            )
            
            # Run make
            print("\nRunning make...")
            subprocess.check_call(
                ['make', 'VERBOSE=1'],
                cwd=self.build_temp
            )
            
            # Find the built extension file
            ext_path = None
            for root, _, files in os.walk(self.build_temp):
                for file in files:
                    if file.endswith('.so') and 'dr_sasa_py' in file:
                        ext_path = os.path.join(root, file)
                        break
                if ext_path:
                    break
                    
            if not ext_path:
                raise RuntimeError("Built extension file not found!")
                
            # Create the final directory if it doesn't exist
            os.makedirs(os.path.dirname(self.get_ext_fullpath(ext.name)), exist_ok=True)
            
            # Copy the built extension to its final location
            print(f"\nCopying {ext_path} to {self.get_ext_fullpath(ext.name)}")
            self.copy_file(ext_path, self.get_ext_fullpath(ext.name))
            
        except subprocess.CalledProcessError as e:
            print(f"Build failed with error code {e.returncode}")
            print("Command output:")
            print(e.output if hasattr(e, 'output') else "No output available")
            raise
        except Exception as e:
            print(f"Error during build: {str(e)}")
            raise

setup(
    name='dr_sasa',
    version='0.1',
    author='Alessio DAddio',
    description='Python bindings for dr_sasa',
    ext_modules=[CMakeExtension('dr_sasa_py')],
    cmdclass={'build_ext': CMakeBuild},
    install_requires=[
        'numpy',
        'pytest'
    ],
)