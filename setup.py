from skbuild import setup

if __name__ == "__main__":
    setup(
        packages=["dr_sasa_python"],
        cmake_install_dir="dr_sasa_python",
        cmake_args=[
            "-DCMAKE_BUILD_TYPE=Release",
            "-DBUILD_SHARED_LIBS=ON",
        ],
        cmake_source_dir=".",
        cmake_languages=("C", "CXX"),
    )