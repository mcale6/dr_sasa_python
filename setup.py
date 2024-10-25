from setuptools import setup
from setuptools.command.build_ext import build_ext
from pathlib import Path
import platform
import subprocess
import sys

try:
    from pybind11.setup_helpers import Pybind11Extension, build_ext
except ImportError:
    from setuptools import Extension as Pybind11Extension


class CMakeExtension(Pybind11Extension):
    def __init__(self, name):
        super().__init__(name, [])


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        # Convert paths to Path objects for better handling
        project_root = Path(__file__).parent.absolute()
        build_temp = Path(self.build_temp)
        ext_fullpath = Path(self.get_ext_fullpath(ext.name))

        # Create build directory
        build_temp.mkdir(parents=True, exist_ok=True)
        ext_fullpath.parent.mkdir(parents=True, exist_ok=True)

        # Basic CMake configuration
        config = 'Debug' if self.debug else 'Release'
        cmake_args = [
            f'-DCMAKE_BUILD_TYPE={config}',
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={build_temp}',
            f'-DPYTHON_EXECUTABLE={sys.executable}',
        ]

        # Platform-specific configurations
        if platform.system() == "Linux":
            cmake_args.extend([
                '-DCMAKE_C_COMPILER=gcc',
                '-DCMAKE_CXX_COMPILER=g++',
            ])

        try:
            print("Configuring CMake...")
            subprocess.run(
                ['cmake', str(project_root)] + cmake_args,
                cwd=build_temp,
                check=True,
                capture_output=True,
                text=True
            )

            print("Building extension...")
            subprocess.run(
                ['cmake', '--build', '.', '--config', config],
                cwd=build_temp,
                check=True,
                capture_output=True,
                text=True
            )

            # Find and copy the built extension
            built_ext = list(build_temp.rglob(f'*dr_sasa_py*.so'))[0]
            self.copy_file(str(built_ext), str(ext_fullpath))

        except subprocess.CalledProcessError as e:
            print(f"Build failed!")
            print(f"Command output:\n{e.output}")
            raise
        except Exception as e:
            print(f"Error: {e}")
            raise

setup(
    name='dr_sasa',
    version='0.1',
    author='mcale6',
    description='Python bindings for dr_sasa',
    long_description='',
    ext_modules=[CMakeExtension('dr_sasa_py')],
    cmdclass={'build_ext': CMakeBuild},
    zip_safe=False,
    python_requires='>=3.7',
    install_requires=[
        'numpy>=1.20.0',
        'pytest>=6.0.0',
    ],
    extras_require={
        'dev': [
            'pytest',
            'pytest-cov',
            'flake8',
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
    ],
)