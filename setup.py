from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os
import subprocess
import platform

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # Add OpenMP flags for macOS
        if platform.system() == "Darwin":
            try:
                libomp_path = subprocess.check_output(['brew', '--prefix', 'libomp'], 
                                                    universal_newlines=True).strip()
                print(f"Found libomp at: {libomp_path}")
                openmp_flags = [
                    f'-DOpenMP_C_FLAGS=-Xpreprocessor -fopenmp -I{libomp_path}/include',
                    f'-DOpenMP_CXX_FLAGS=-Xpreprocessor -fopenmp -I{libomp_path}/include',
                    f'-DOpenMP_C_LIB_NAMES=omp',
                    f'-DOpenMP_CXX_LIB_NAMES=omp',
                    f'-DOpenMP_omp_LIBRARY={libomp_path}/lib/libomp.dylib',
                ]
            except subprocess.CalledProcessError:
                print("Warning: libomp not found. Please install with 'brew install libomp'")
                openmp_flags = []
        else:
            openmp_flags = []

        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DPYTHON_EXECUTABLE=' + sys.executable,
        ] + openmp_flags

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
        build_args += ['--', '-j2']

        env = os.environ.copy()
        if platform.system() == "Darwin":
            try:
                env['LDFLAGS'] = f"-L{libomp_path}/lib"
                env['CPPFLAGS'] = f"-I{libomp_path}/include"
            except:
                pass

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        print("CMake args:", cmake_args)
        print("Build args:", build_args)
        print("Environment:", env)
            
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, 
                            cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                            cwd=self.build_temp)

setup(
    name='dr_sasa',
    version='0.1',
    author='Alessio DAddio',
    description='Python bindings for dr_sasa',
    long_description='',
    ext_modules=[CMakeExtension('dr_sasa_py')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'pytest'
    ],
)