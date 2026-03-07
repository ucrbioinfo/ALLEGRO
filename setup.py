from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import subprocess, pathlib, shutil, os, re, stat, sys

root = pathlib.Path(__file__).resolve().parent
long_description = (root / "README.md").read_text() if (root / "README.md").exists() else ""

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        src_dir    = root / "cpp"
        build_dir  = pathlib.Path(self.build_temp) / "cmake-build"
        build_dir.mkdir(parents=True, exist_ok=True)

        # Standardizing OR-Tools build parameters
        ortools_args = [
            "-DBUILD_DEPS=ON",
            "-DUSE_COINOR=OFF", "-DUSE_HIGHS=OFF", "-DUSE_SCIP=OFF",
            "-DUSE_BOP=OFF", "-DUSE_MATH_OPT=OFF",
            "-DUSE_PDLP=OFF", "-DBUILD_TESTING=OFF",
            "-DBUILD_SAMPLES=OFF", "-DBUILD_EXAMPLES=OFF",
            "-DUSE_DOTNET_8=OFF",
        ]
        
        cmake_args = [
            "-DCMAKE_BUILD_TYPE=Release",
            f"-DPython3_EXECUTABLE={sys.executable}",
        ]

        print(f"[BUILD] Configuring in {build_dir}")
        subprocess.check_call(["cmake", "-S", str(src_dir), "-B", str(build_dir)] + ortools_args + cmake_args)
        
        print("[BUILD] Compiling C++ core and executable...")
        jobs = os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL", "4")
        subprocess.check_call([
            "cmake", "--build", str(build_dir), 
            "--target", "_kirschtorte_core", "uCRISPR_scorer", 
            "--parallel", jobs
        ])

        # 1) Install the C++ Engine (_kirschtorte_core)
        ext_path = pathlib.Path(self.get_ext_fullpath(ext.name))
        ext_path.parent.mkdir(parents=True, exist_ok=True)
        
        so_candidates = list(build_dir.glob("_kirschtorte_core*.so"))
        if not so_candidates:
            so_candidates = list(build_dir.glob("**/_kirschtorte_core*.so"))
        
        if not so_candidates:
            raise RuntimeError(f"Could not find _kirschtorte_core*.so in {build_dir}")
        
        shutil.copyfile(so_candidates[0], ext_path)
        os.chmod(ext_path, os.stat(ext_path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        # 2) Stage uCRISPR_scorer into allegro/bin
        pkg_stage = pathlib.Path(self.build_lib) / "allegro"
        bin_out   = pkg_stage / "bin"
        bin_out.mkdir(parents=True, exist_ok=True)

        exe_candidates = list(build_dir.glob("uCRISPR_scorer"))
        if not exe_candidates:
            exe_candidates = list(build_dir.glob("**/uCRISPR_scorer"))

        if not exe_candidates:
            raise RuntimeError(f"uCRISPR_scorer executable not found in {build_dir}")
        
        target_exe = bin_out / "uCRISPR_scorer"
        shutil.copy2(exe_candidates[0], target_exe)

        # RPATH fixing for Linux portability
        patchelf = shutil.which("patchelf")
        if patchelf:
            subprocess.check_call([patchelf, "--set-rpath", "$ORIGIN/../lib:$ORIGIN", str(target_exe)])
        
        # Optimization: Strip unneeded symbols
        strip = shutil.which("strip")
        if strip:
            for path in [ext_path, target_exe]:
                try: 
                    subprocess.check_call([strip, "--strip-unneeded", str(path)])
                except Exception: 
                    pass

setup(
    name="allegro-bio",
    version="1.0.0",
    description="ALLEGRO: CRISPR Guide Design Tool",
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(where="python"),
    package_dir={"": "python"},
    include_package_data=True,
    package_data={
        "allegro": [
            "bin/uCRISPR_scorer",
            "lib/*.so*",
            "_kirschtorte.py",
            "scorers/uCRISPR/RNAstructure/data_tables/**/*"
        ],
    },
    cmdclass={"build_ext": CMakeBuild},
    ext_modules=[Extension("allegro._kirschtorte_core", sources=[])],
    entry_points={
        "console_scripts": [
            "allegro = allegro.main:main",
        ]
    },
    install_requires=["numpy", "pandas", "pyyaml", "setproctitle", "scipy", "biopython"],
    python_requires=">=3.10",
)