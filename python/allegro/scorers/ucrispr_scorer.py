import os
import pickle
import subprocess
import scipy.stats

from pathlib import Path
from importlib import resources as ir

from allegro.utils.shell_colors import bcolors
from allegro.scorers.scorer_base import Scorer

class uCRISPR_scorer(Scorer):
    def __init__(self) -> None:
        self.use_secondary_memory: bool = True
        self.context_toward_five_prime: int = 4
        self.context_toward_three_prime: int = 3

        # Resolve internal paths using modern importlib.resources
        # self._pkg_root is a Traversable object
        self._pkg_root = ir.files("allegro") 
        self.cpp_program_path = self._pkg_root / "bin" / "uCRISPR_scorer"
        self._datapath = self._pkg_root / "scorers" / "uCRISPR" / "RNAstructure" / "data_tables"

        if self.use_secondary_memory:
            if not os.path.exists('allegro_cache/'):
                os.makedirs('allegro_cache/', exist_ok=True)

    def _run_ucrispr(self, payload: bytes) -> list[float]:
        # Convert Traversable to string for subprocess and OS calls
        exe = str(self.cpp_program_path)
        
        if not os.path.exists(exe):
            print(f"{bcolors.RED}> Dev error{bcolors.RESET}: uCRISPR executable not found at {exe}")
            sys.exit(1)

        # Handle permissions (common in some wheel installations)
        if not os.access(exe, os.X_OK):
            try:
                os.chmod(exe, os.stat(exe).st_mode | 0o111)
            except Exception:
                pass

        if not os.path.exists(str(self._datapath)):
            print(f"{bcolors.RED}> Dev error{bcolors.RESET}: RNAstructure data_tables not found.")
            sys.exit(1)

        env = os.environ.copy()
        env["DATAPATH"] = str(self._datapath)

        # Auditwheel repairs often move libs to allegro.libs or a custom lib folder
        # Adding our internal lib folder to LD_LIBRARY_PATH ensures the binary works
        libdir = str(self._pkg_root / "lib")
        env["LD_LIBRARY_PATH"] = f"{libdir}:{env.get('LD_LIBRARY_PATH', '')}".rstrip(":")

        try:
            process = subprocess.Popen(
                [exe], 
                stdin=subprocess.PIPE, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                env=env
            )
        except PermissionError:
            print(f"{bcolors.RED}> Error{bcolors.RESET}: Permission denied for {exe}")
            sys.exit(1)

        output, stderr = process.communicate(payload)

        if stderr:
            print(f"{bcolors.RED}>{bcolors.RESET} uCRISPR_scorer error: {stderr.decode(errors='replace')}")
            sys.exit(1)

        lines = output.decode().strip().split("\n") if output else []
        return [float(s) for s in lines if s]

    def score_sequence(self, guides_context_list: list[str]) -> list[float]:
        scores: list = list()
        saved_guides: dict[str, float] = dict()
        indices_of_uncached_guides: list[int] = list()

        if self.use_secondary_memory:
            cache_path = f"allegro_cache/cached_guides.pickle"
            if os.path.exists(cache_path):
                with open(cache_path, "rb") as f:
                    saved_guides = pickle.load(f)

        for idx, guide in enumerate(guides_context_list):
            score = saved_guides.get(guide)
            scores.append(score)
            if score is None:
                indices_of_uncached_guides.append(idx)

        if len(indices_of_uncached_guides) > 0:
            for i in range(0, len(indices_of_uncached_guides), 1000):
                payload = b""
                for j in indices_of_uncached_guides[i:i+1000]:
                    payload += guides_context_list[j].encode() + b"\n"

                raw_scores = self._run_ucrispr(payload)

                uncached_scores = [
                    scipy.stats.norm.cdf(float(s), loc=11.92658, scale=0.2803797) * 100
                    for s in raw_scores
                ]

                for idx, s in enumerate(uncached_scores):
                    scores[indices_of_uncached_guides[idx+i]] = s
                    saved_guides[guides_context_list[indices_of_uncached_guides[idx+i]]] = s

            if self.use_secondary_memory:
                cache_path = f"allegro_cache/cached_guides.pickle"
                with open(cache_path, "wb") as f:
                    pickle.dump(saved_guides, f)

        return scores
