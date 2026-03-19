from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import subprocess, pathlib, shutil, os, sys

root = pathlib.Path(__file__).resolve().parent
long_description = (root / "README.md").read_text() if (root / "README.md").exists() else ""

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        # 1. Paths setup
        src_dir    = root / "cpp" 
        build_dir  = pathlib.Path(self.build_temp) / "cmake-build"
        build_dir.mkdir(parents=True, exist_ok=True)

        # 2. Configure and Build
        cmake_args = [
            "-DCMAKE_BUILD_TYPE=Release",
            f"-DPython3_EXECUTABLE={sys.executable}",
        ]

        print(f"[BUILD] Configuring in {build_dir}")
        subprocess.check_call(["cmake", "-S", str(src_dir), "-B", str(build_dir)] + cmake_args)
        
        print("[BUILD] Compiling C++ core and executable...")
        jobs = os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL", "4")
        subprocess.check_call([
            "cmake", "--build", str(build_dir), 
            "--target", "_kirschtorte_core", "uCRISPR_scorer", 
            "--parallel", jobs
        ])

        # 3. Install the Python Extension (_kirschtorte_core)
        ext_path = pathlib.Path(self.get_ext_fullpath(ext.name))
        ext_path.parent.mkdir(parents=True, exist_ok=True)
        
        so_candidates = list(build_dir.glob("**/_kirschtorte_core*"))
        if not so_candidates:
            raise RuntimeError(f"Could not find _kirschtorte_core extension in {build_dir}")
        
        shutil.copyfile(so_candidates[0], ext_path)

        # 4. Stage uCRISPR_scorer into 'allegro/bin'
        pkg_stage = pathlib.Path(self.build_lib) / "allegro"
        bin_out   = pkg_stage / "bin"
        bin_out.mkdir(parents=True, exist_ok=True)

        exe_candidates = list(build_dir.glob("**/uCRISPR_scorer"))
        if not exe_candidates:
             raise RuntimeError("uCRISPR_scorer executable not found")
        
        target_exe = bin_out / "uCRISPR_scorer"
        shutil.copy2(exe_candidates[0], target_exe)

        # 5. Optimization: Strip symbols
        strip = shutil.which("strip")
        if strip:
            for path in [ext_path, target_exe]:
                try: 
                    subprocess.check_call([strip, "--strip-unneeded", str(path)])
                except Exception: pass

setup(
    name="allegro-bio",
    version="1.0.0",
    description="ALLEGRO: CRISPR guide design tool with ILP/CP-SAT core",
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(where="python"),
    package_dir={"": "python"},
    include_package_data=True,
    package_data={
        "allegro": [
            "bin/uCRISPR_scorer",
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