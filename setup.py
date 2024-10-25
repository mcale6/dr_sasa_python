from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os
import subprocess

class CMakeExtension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])

class CMakeBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        # Where to put the library
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # Configure cmake
        cmake_args = [
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}',
            f'-DPYTHON_EXECUTABLE={sys.executable}',
            '-DCMAKE_BUILD_TYPE=Debug'
        ]

        # Build
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
            
        subprocess.check_call(['cmake', '..'] + cmake_args, cwd=self.build_temp)
        subprocess.check_call(['cmake', '--build', '.'], cwd=self.build_temp)

setup(
    name='dr_sasa',
    version='0.1',
    author='Alessio DAddio',
    description='Python bindings for dr_sasa',
    ext_modules=[CMakeExtension('dr_sasa_py')],
    cmdclass=dict(build_ext=CMakeBuild),
    install_requires=['numpy', 'pytest'],
)