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
            
        # Get OpenMP settings for macOS
        if platform.system() == "Darwin":
            try:
                libomp_path = subprocess.check_output(
                    ['brew', '--prefix', 'libomp'], 
                    universal_newlines=True
                ).strip()
                print(f"Found libomp at: {libomp_path}")
                
                # Set environment variables for OpenMP
                os.environ['CC'] = 'clang'
                os.environ['CXX'] = 'clang++'
                os.environ['CPPFLAGS'] = f"-Xpreprocessor -fopenmp -I{libomp_path}/include"
                os.environ['CFLAGS'] = f"-Xpreprocessor -fopenmp -I{libomp_path}/include"
                os.environ['CXXFLAGS'] = f"-Xpreprocessor -fopenmp -I{libomp_path}/include"
                os.environ['LDFLAGS'] = f"-Wl,-rpath,{libomp_path}/lib {libomp_path}/lib/libomp.dylib"
                
                omp_flags = [
                    f'-DOpenMP_C_FLAGS=-Xpreprocessor -fopenmp -I{libomp_path}/include',
                    f'-DOpenMP_CXX_FLAGS=-Xpreprocessor -fopenmp -I{libomp_path}/include',
                    f'-DOpenMP_C_LIB_NAMES=omp',
                    f'-DOpenMP_CXX_LIB_NAMES=omp',
                    f'-DOpenMP_omp_LIBRARY={libomp_path}/lib/libomp.dylib',
                    f'-DOpenMP_CXX_LIB_NAMES=omp',
                    f'-DOpenMP_CXX_LIBRARIES={libomp_path}/lib/libomp.dylib'
                ]
            except subprocess.CalledProcessError:
                print("Warning: libomp not found, please install with: brew install libomp")
                omp_flags = []
                
            # Set architecture
            if platform.machine() == 'arm64':
                arch_flag = "-DCMAKE_OSX_ARCHITECTURES=arm64"
            else:
                arch_flag = "-DCMAKE_OSX_ARCHITECTURES=x86_64"
        else:
            omp_flags = []
            arch_flag = ""
            
        # Configure cmake arguments
        cmake_args = [
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={self.build_temp}',
            f'-DPYTHON_EXECUTABLE={sys.executable}',
            f'-DCMAKE_BUILD_TYPE=Debug'
        ]
        
        if arch_flag:
            cmake_args.append(arch_flag)
        cmake_args.extend(omp_flags)
        
        print(f"Project root: {project_root}")
        print(f"Build temp: {self.build_temp}")
        print(f"CMake args: {cmake_args}")
        print(f"Environment:")
        for key in ['CC', 'CXX', 'CPPFLAGS', 'CFLAGS', 'CXXFLAGS', 'LDFLAGS']:
            print(f"  {key}: {os.environ.get(key, 'not set')}")
        
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