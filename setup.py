from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
import sys
import os
import subprocess

class CMakeExtension:
    def __init__(self, name, sourcedir=""):
        self.name = name
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # Required for auto-detection & inclusion of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead of Python3_EXECUTABLE
        # to avoid CMake warning
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
        ]
        build_args = []

        # Default to ninja-build if available
        if not cmake_generator and os.system("which ninja") == 0:
            cmake_args += ["-GNinja"]

        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control parallel build
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            build_args += ["-j"]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)

setup(
    name="dr_sasa",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Python bindings for dr_sasa library",
    long_description="",
    packages=find_packages(),
    ext_modules=[CMakeExtension("dr_sasa_py")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
    ],
    extras_require={
        "test": ["pytest>=6.0"],
        "dev": ["pytest>=6.0", "pytest-cov", "black", "flake8"],
    },
)