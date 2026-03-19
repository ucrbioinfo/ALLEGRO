import os
import re
import sys
import time
import yaml
import pandas
import signal
import shutil
import argparse
from datetime import timedelta

from allegro.utils.iupac import iupac_dict
from allegro.utils.shell_colors import bcolors
from allegro.scorers.scorer_factory import scorer_names

def normalize_scorer(scorer: str) -> str:
    normalized_scorer = scorer.strip().lower()
    if normalized_scorer not in scorer_names:
        print(f'{bcolors.RED}> Error{bcolors.RESET}: Scorer {scorer} is not supported. Exiting.')
        sys.exit(1)
    return normalized_scorer

def parse_pam(v: str) -> str:
    if v is None:
        return ''

    pam = str(v).strip().upper()

    if len(pam) == 0:
        print(f'{bcolors.RED}> Error{bcolors.RESET}: PAM must be specified. Exiting.')
        sys.exit(1)

    for c in pam:
        if c not in iupac_dict:
            print(f'{bcolors.RED}> Error{bcolors.RESET}: Invalid IUPAC code "{c}" in PAM "{pam}". Exiting.')
            sys.exit(1)

    return pam

def require_on_path(binary_name: str, hint: str = "") -> None:
    if shutil.which(binary_name) is not None:
        return

    print(f"{bcolors.RED}> Error{bcolors.RESET}: Required executable '{binary_name}' was not found on your PATH. Exiting.")
    if hint:
        print(f"{bcolors.BLUE}>{bcolors.RESET} {hint}")
    sys.exit(1)

def coerce_bool(v) -> bool:
    if isinstance(v, bool):
        return v

    if v is None:
        print(f"{bcolors.RED}> Error{bcolors.RESET}: Expected boolean (True/False), got None. Exiting.")
        sys.exit(1)

    s = str(v).strip().lower()
    if s in ("true", "1", "yes", "y", "t"):
        return True
    if s in ("false", "0", "no", "n", "f"):
        return False

    print(f"{bcolors.RED}> Error{bcolors.RESET}: Expected boolean (True/False), got: {v}. Exiting.")
    print(f"{bcolors.BLUE}>{bcolors.RESET} Accepted values: True/False, 1/0, yes/no, y/n, t/f")
    sys.exit(1)

def sanitize_filename(filename, max_length=255) -> str:
    if not filename or filename == '':
        return 'ALLEGRO_TEST_RUN'

    # Remove disallowed characters (e.g., /, \0, *, ?, ", <, >, |)
    sanitized = re.sub(r'[\/\0\*\?"<>\|]', '', filename)
    
    # Replace whitespace characters with underscores
    sanitized = re.sub(r'\s+', '_', sanitized)
    
    # Remove leading periods to avoid hidden files, unless it's a special case (e.g., ".", "..")
    if sanitized.startswith('.') and sanitized not in ['.', '..']:
        sanitized = sanitized.lstrip('.')

    # Böse böse...
    if sanitized == '.':
        print(f'{bcolors.ORANGE}> Böse böse...{bcolors.RESET}')
        sanitized = '._'
    elif sanitized == '..':
        print(f'{bcolors.ORANGE}> Böse böse...{bcolors.RESET}')
        sanitized = '.._'
    
    # Shorten the filename to comply with filesystem limits
    if len(sanitized.encode('utf-8')) > max_length:
        # Shorten while trying to preserve file extension
        name, dot, extension = sanitized.rpartition('.')
        if dot and extension and len(extension) <= max_length - 1:
            # Shorten name part, leave room for dot and extension
            name = name[:max_length - len(extension) - 2]
            sanitized = f"{name}.{extension}"
        else:
            # If there's no extension or it's too long, just truncate the name
            sanitized = sanitized[:max_length]
    
    return sanitized

def signal_handler(sig, frame) -> None:
    print(f'\n{bcolors.BLUE}>{bcolors.RESET} Interrupted {bcolors.ORANGE}ALLEGRO{bcolors.RESET}. Ciao.')
    sys.exit(0)

class Configurator:
    def __init__(self) -> None:
        self.start_time = time.time()
        self.species_column_name = 'species_name';
        
        signal.signal(signal.SIGINT, signal_handler)

    def begruessung(self) -> None:
        print(f'{bcolors.BLUE}>{bcolors.RESET} Welcome to {bcolors.ORANGE}ALLEGRO{bcolors.RESET}.')

    def parse_configurations(self) -> None:
        # ---------------------------
        # 1. Pre-parse only --config
        # ---------------------------
        config_parser = argparse.ArgumentParser(add_help=False)

        config_parser.add_argument(
            "-c",
            "--config",
            type=str,
            default="config.yaml",
            help="Config file to use."
        )

        config_args, remaining_args = config_parser.parse_known_args()

        # ---------------------------
        # 2. Load YAML if present
        # ---------------------------
        config_defaults = {}

        if os.path.exists(config_args.config):
            try:
                with open(config_args.config, "r") as f:
                    config_defaults = yaml.safe_load(f) or {}
            except yaml.YAMLError as e:
                m = re.search(r"line (\d+)", str(e))
                line_number = m.group(1) if m else "?"
                print(f"{bcolors.RED}> Error{bcolors.RESET}: Problem reading {config_args.config} (line {line_number})")
                print(e)
                sys.exit(1)

        parser = argparse.ArgumentParser(
            prog="ALLEGRO",
            description="ALLEGRO: Find a small set of guide RNAs that covers all targets.",
            parents=[config_parser],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )

        parser.add_argument(
            "--help-examples",
            action="store_true",
            help="Show usage examples and exit."
        )

        g_general = parser.add_argument_group("General")
        g_inputs  = parser.add_argument_group("Inputs")
        g_design  = parser.add_argument_group("Design objective")
        g_scoring = parser.add_argument_group("Scoring")
        g_ot      = parser.add_argument_group("Off-target reporting")
        g_cluster = parser.add_argument_group("Clustering")
        g_output  = parser.add_argument_group("Output")
        g_adv     = parser.add_argument_group("Advanced")

        g_general.add_argument("--soundcheck", action="store_true", help="Exit after basic sanity checks.")
        g_general.add_argument("-n", "--experiment_name", type=sanitize_filename, help="Name of the experiment.")

        g_inputs.add_argument("-id", "--input_directory", type=str, default="data/input/", help="Directory containing input .fna files.")
        g_inputs.add_argument("-isp", "--input_species_path", type=str, default="data/input/manifest.csv", help="CSV mapping species to genome/CDS files.")
        g_inputs.add_argument("-ispc", "--input_species_path_column", type=str, default="file_name", help="Column name for file selection.")

        g_design.add_argument("-t", "--track", type=str, default="e", help="Track: e/track_e or a/track_a.")
        g_design.add_argument("-m", "--multiplicity", type=int, default=1, help="Require each target to be hit at least this many times.")
        g_design.add_argument("-b", "--beta", type=int, default=0, help="Soft target for max guides; 0 disables beta objective.")

        g_scoring.add_argument("-s", "--scorer", type=normalize_scorer, default="dummy", help="Scorer: dummy, ucrispr.")
        g_scoring.add_argument("-sthresh", "--guide_score_threshold", type=int, default=0, help="Discard guides below this score.")

        g_scoring.add_argument("-gc", "--filter_by_gc", type=coerce_bool, default=True, help="Filter guides by GC range.")
        g_scoring.add_argument("-gcmin", "--gc_min", type=float, default=0.4, help="Minimum GC fraction.")
        g_scoring.add_argument("-gcmax", "--gc_max", type=float, default=0.6, help="Maximum GC fraction.")
        g_scoring.add_argument("-pte", "--patterns_to_exclude", type=lambda s: s.split(','), default=[], help="Comma-separated IUPAC patterns to exclude.")

        g_ot.add_argument("-off", "--output_offtargets", type=coerce_bool, default=False, help="Generate off-target report. Requires Bowtie v1.")
        g_ot.add_argument("-reportmm", "--report_up_to_n_mismatches", type=int, default=3, help="Max mismatches after seed region (0-3). Requires output_offtargets.")
        g_ot.add_argument("-seed", "--seed_region_is_n_upstream_of_pam", type=int, default=12, help="Seed length upstream of PAM. Requires output_offtargets.")
        g_ot.add_argument("-isod", "--input_species_offtarget_dir", type=str, default="", help="Dir with background FASTAs. Requires output_offtargets.")
        g_ot.add_argument("-isoc", "--input_species_offtarget_column", type=str, default="", help="CSV column naming background FASTA. Requires output_offtargets.")

        g_cluster.add_argument("-prec", "--preclustering", type=coerce_bool, default=False, help="Enable preclustering (performance tradeoff). Requires Bowtie v1.")
        g_cluster.add_argument("-postc", "--postclustering", type=coerce_bool, default=False, help="Enable postclustering (compress output). Requires Bowtie v1.")
        g_cluster.add_argument("-mmafterseed", "--mismatches_allowed_after_seed_region", type=int, default=2, help="Allowed mismatches after seed for clustering. Requires a clustering flag.")

        g_output.add_argument("-od", "--output_directory", type=str, default="data/output/", help="Output base directory.")
        g_adv.add_argument("-esp", "--early_stopping_patience", type=int, default=60, help="ILP early-stopping patience (seconds).")
        g_adv.add_argument("-esd", "--enable_solver_diagnostics", type=coerce_bool, default=True, help="Try to diagnose infeasible instances.")
        g_adv.add_argument("--pam", type=parse_pam, default="NGG", help='PAM. Accepts "NGG".')
        g_adv.add_argument("-pl", "--protospacer_length", type=int, default=20, help='Protospacer length. Defaults to 20.')
        g_adv.add_argument("--mp_threshold", type=int, default=0, help="Preselect guides that hit <= this many targets.")
        g_adv.add_argument("--align_solution_to_input", type=coerce_bool, default=True, help="Align output guides back to where they came from for CSV report. Requires Bowtie v1.")
        
        parser.set_defaults(**config_defaults)
        
        args = parser.parse_args(remaining_args)

        if args.help_examples:
            print(f"{bcolors.BLUE}>{bcolors.RESET} Examples")
            print("  allegro -c config.yaml")
            print("  allegro -t e -isp species.csv -id data/input/")
            sys.exit(0)

        self.args = args

    def check_and_fix_configurations(self) -> argparse.Namespace:
        if self.args.soundcheck:
            print(f'{bcolors.BLUE}> {bcolors.RESET}{bcolors.ORANGE}ALLEGRO{bcolors.RESET} is ready to shred.')
            sys.exit(0)

        print(f'{bcolors.BLUE}>{bcolors.RESET} All unspecified command-line arguments default to the values in {self.args.config}.')
        
        species_df = None

        # ------------------------------------------------------------------------------
        #   input_directory
        # ------------------------------------------------------------------------------
        if not os.path.exists(self.args.input_directory):
            print(f'{bcolors.RED}> Error{bcolors.RESET}: Cannot find directory {self.args.input_directory}. Does it exist, and did you spell its name (input_directory) correctly? Exiting.')
            sys.exit(1)
        
        self.blocklist_path = os.path.join(self.args.input_directory, "_the_blocklist_.txt")

        # ------------------------------------------------------------------------------
        #   input_species_path
        # ------------------------------------------------------------------------------
        try:
            species_df = pandas.read_csv(self.args.input_species_path)
        except pandas.errors.EmptyDataError:
            print(f'{bcolors.RED}> Error{bcolors.RESET}: File {self.args.input_species_path} is empty. Exiting.')
            sys.exit(1)
        except FileNotFoundError:
            print(f'{bcolors.RED}> Error{bcolors.RESET}: Cannot find file {self.args.input_species_path}. Does it exist, and did you spell the path/file name (input_species_path) correctly? Exiting.')
            sys.exit(1)
        
        # ------------------------------------------------------------------------------
        #   input_species_path_column
        # ------------------------------------------------------------------------------
        if self.args.input_species_path_column not in species_df.columns:
            print(f'{bcolors.RED}> Error{bcolors.RESET}: Cannot find column the specified "{self.args.input_species_path_column}" in {self.args.input_species_path}. Does it exist, and did you spell the column name (input_species_path_column) correctly? Exiting.')
            sys.exit(1)
        
        # ------------------------------------------------------------------------------
        #   species_name
        # ------------------------------------------------------------------------------
        if "species_name" not in species_df.columns:
            print(f'{bcolors.RED}> Error{bcolors.RESET}: Cannot find the required column "species_name" in {self.args.input_species_path}. Does it exist, and did you format the CSV file correctly? Exiting.')
            sys.exit(1)

        # ------------------------------------------------------------------------------
        #   align_solution_to_input or preclustering or postclustering or output_offtargets 
        #   -> requires bowtie on PATH
        # ------------------------------------------------------------------------------
        if self.args.align_solution_to_input or self.args.preclustering or self.args.postclustering or self.args.output_offtargets:
            require_on_path(
                "bowtie",
                hint="Any of parameters align_solution_to_input, preclustering, output_offtargets, or postclustering require Bowtie v1. Either turn them off (set to False), or install Bowtie v1 and ensure 'bowtie' is available. Example: conda install -c bioconda bowtie"
            )
            # also require bowtie-build:
            require_on_path(
                "bowtie-build",
                hint="Any of parameters align_solution_to_input, preclustering, output_offtargets, or postclustering require Bowtie v1. Either turn them off (set to False), or install Bowtie v1 and ensure 'bowtie-build' is available. Example: conda install -c bioconda bowtie"
            )
        
        # ------------------------------------------------------------------------------
        #   patterns_to_exclude
        # ------------------------------------------------------------------------------
        pte = self.args.patterns_to_exclude

        # Normalize YAML + CLI shapes into a flat list[str]
        if not pte:
            self.args.patterns_to_exclude = []
        elif isinstance(pte, str):
            self.args.patterns_to_exclude = [pte]
        else:
            flat: list[str] = []
            for x in pte:
                if isinstance(x, list):
                    flat.extend(x)
                else:
                    flat.append(x)
            self.args.patterns_to_exclude = flat

        # Strip whitespace, uppercase, drop empties
        self.args.patterns_to_exclude = [
            str(p).strip().upper() for p in self.args.patterns_to_exclude
            if str(p).strip() != ""
        ]

        # Validate patterns
        bad_patterns: set[str] = set()

        for pattern_to_exclude in self.args.patterns_to_exclude:
            # Do not allow more than 5 non-ACTG characters in exclusion patterns
            if sum(c not in ['A', 'C', 'T', 'G'] for c in pattern_to_exclude) > 5:
                print(f'{bcolors.RED}> Error{bcolors.RESET}: Cannot have more than 5 non-ACTG IUPAC codes in pattern "{pattern_to_exclude}" due to the exponential number of combinations. Remove it and retry. Exiting.')
                sys.exit(1)
            
            if len(pattern_to_exclude) > self.args.protospacer_length:
                print(f'{bcolors.ORANGE}> Warning{bcolors.RESET}: Pattern "{pattern_to_exclude}" is longer than the length of your protospacer ({self.protospacer_pattern} – {len(pattern_to_exclude)} vs. {protospacer_length}). Ignoring this pattern.')
                bad_patterns.add(pattern_to_exclude)
                continue

            for character in pattern_to_exclude:
                if character not in iupac_dict:
                    print(f'{bcolors.ORANGE}> Warning{bcolors.RESET}: patterns_to_exclude pattern "{pattern_to_exclude}" contains "{character}" which is not an IUPAC code. Ignoring this pattern.')
                    bad_patterns.add(pattern_to_exclude)
                    break

        # Remove invalid patterns
        if bad_patterns:
            self.args.patterns_to_exclude = [
                p for p in self.args.patterns_to_exclude if p not in bad_patterns
            ]

        # ------------------------------------------------------------------------------
        #   output_offtargets
        # ------------------------------------------------------------------------------
        if self.args.output_offtargets:
            try:
                species_df[self.args.input_species_offtarget_column]
            except KeyError:
                print(f'{bcolors.RED}> Error{bcolors.RESET}: Cannot find column "{self.args.input_species_offtarget_column}" in {self.args.input_species_path}. Did you format the CSV file correctly? Exiting.')
                sys.exit(1)
        
            if self.args.report_up_to_n_mismatches > 3 or self.args.report_up_to_n_mismatches < 0:
                print(f'{bcolors.RED}> Error{bcolors.RESET}: The value for report_up_to_n_mismatches may be 0, 1, 2, or 3. Exiting.')
                sys.exit(1)

            if self.args.seed_region_is_n_upstream_of_pam < 0:
                print(f'{bcolors.RED}> Warning{bcolors.RESET}: seed_region_is_n_upstream_of_pam ({self.args.seed_region_is_n_upstream_of_pam}) is set to a negative value. Auto adjusting to 0.')
                self.args.seed_region_is_n_upstream_of_pam = 0

            if not os.path.exists(self.args.input_species_offtarget_dir):
                print(f'{bcolors.RED}> Error{bcolors.RESET}: Path {self.args.input_species_offtarget_dir} does not exist. Exiting.')
                sys.exit(1)
            elif len(os.listdir(self.args.input_species_offtarget_dir)) == 0:
                print(f'{bcolors.RED}> Error{bcolors.RESET}: Directory {self.args.input_species_offtarget_dir} exists but is empty. Exiting.')
                sys.exit(1)

            base_path = os.getcwd()
            for file_name in species_df[self.args.input_species_offtarget_column]:
                offtarget_file_full_path = os.path.join(base_path, self.args.input_species_offtarget_dir, file_name)
                if not os.path.exists(offtarget_file_full_path):
                    print(f'{bcolors.RED}> Error{bcolors.RESET}: Using the given configuration in {self.args.config} under "output_offtargets: True," path {base_path}/{self.args.input_species_offtarget_dir}/{file_name} does not exist. Exiting.')
                    sys.exit(1)
                elif os.stat(offtarget_file_full_path).st_size == 0:
                    print(f'{bcolors.RED}> Error{bcolors.RESET}: Using the given configuration in {self.args.config} under "output_offtargets: True," path {base_path}/{self.args.input_species_offtarget_dir}/{file_name} exists, but is empty. Exiting.')
                    sys.exit(1)

        # ------------------------------------------------------------------------------
        #   gc_max
        # ------------------------------------------------------------------------------
        if self.args.gc_max < self.args.gc_min:
            print(f'{bcolors.RED}> Error{bcolors.RESET}: gc_max ({self.args.gc_max}) is set to be lower than gc_min ({self.args.gc_min}). Fix these values in your config. Exiting.')
            sys.exit(1)

        # ------------------------------------------------------------------------------
        #   gc_min
        # ------------------------------------------------------------------------------
        if self.args.gc_min < 0:
            print(f'{bcolors.RED}> Warning{bcolors.RESET}: gc_min ({self.args.gc_min}) is set to a negative value. Auto adjusting to 0.')
            self.args.gc_min = 0

        # ------------------------------------------------------------------------------
        #   report_up_to_n_mismatches
        # ------------------------------------------------------------------------------
        self.args.report_up_to_n_mismatches = int(self.args.report_up_to_n_mismatches)

        if (self.args.report_up_to_n_mismatches > 3) and (self.args.output_offtargets == True):
            print(f'{bcolors.RED}> Error{bcolors.RESET}: report_up_to_n_mismatches ({self.args.report_up_to_n_mismatches}) is set to a value higher than 3. It may only be in range [0-3]. Adjust this value in the config and try again. Exiting.')
            sys.exit(1)

        # ------------------------------------------------------------------------------
        #   mismatches_allowed_after_seed_region
        # ------------------------------------------------------------------------------
        self.args.mismatches_allowed_after_seed_region = int(self.args.mismatches_allowed_after_seed_region)

        if self.args.mismatches_allowed_after_seed_region < 0:
            if self.args.preclustering == False:
                print(f'{bcolors.RED}> Warning{bcolors.RESET}: mismatches_allowed_after_seed_region is set to a negative value ({self.args.mismatches_allowed_after_seed_region}). Auto adjusting to 0.')
                self.args.mismatches_allowed_after_seed_region = 0
            else:
                print(f'{bcolors.RED}> Warning{bcolors.RESET}: mismatches_allowed_after_seed_region is set to a negative value ({self.args.mismatches_allowed_after_seed_region}) with preclustering enabled. Auto adjusting to 1.')
                self.args.mismatches_allowed_after_seed_region = 1
        elif (self.args.mismatches_allowed_after_seed_region > 3) and (self.args.preclustering):
            print(f'{bcolors.RED}> Warning{bcolors.RESET}: preclustering is enabled and mismatches_allowed_after_seed_region is set to a value higher than 3 ({self.args.mismatches_allowed_after_seed_region}). Beware that Bowtie can only report alignments with 3 mismatches. ALLEGRO will allow mismatches in the output guides, but in the targets.csv report, they may not be aligned to their mismatching intended targets.')

        # ------------------------------------------------------------------------------
        #   track
        # ------------------------------------------------------------------------------
        if not self.args.track:
            print(f'{bcolors.RED}> Error{bcolors.RESET}: Track must be set (track_a/a, or track_e/e) – Exiting.')
            sys.exit(1)
        self.args.track = self.args.track.lower()

        if self.args.track not in ['track_a', 'a', 'track_e', 'e']:
            print(f'{bcolors.RED}> Error{bcolors.RESET}: Unknown track "{self.args.track}" selected in {self.args.config}. Options are: "track_a"/"a" and "track_e"/"e" – Exiting.')
            sys.exit(1)
        
        if self.args.track == 'a':
            self.args.track = 'track_a'
        elif self.args.track == 'e':
            self.args.track = 'track_e'

        # ------------------------------------------------------------------------------
        #   multiplicity
        # ------------------------------------------------------------------------------
        self.args.multiplicity = int(self.args.multiplicity)

        if self.args.multiplicity < 1:
            print(f'{bcolors.RED}> Warning{bcolors.RESET}: Multiplicity is set to {self.args.multiplicity}, a value smaller than 1. Re-adjusting multiplicity to 1.')
            self.args.multiplicity = 1

        # ------------------------------------------------------------------------------
        #   mp_threshold
        # ------------------------------------------------------------------------------
        self.args.mp_threshold = int(self.args.mp_threshold)

        if self.args.mp_threshold <= 0:
            print(f'{bcolors.BLUE}>{bcolors.RESET} mp_threshold is set to {self.args.mp_threshold} and thus disabled.')
            self.args.mp_threshold = 0

        # ------------------------------------------------------------------------------
        #   mp_threshold & multiplicity
        # ------------------------------------------------------------------------------
        if self.args.mp_threshold > 0 and self.args.mp_threshold < self.args.multiplicity:
            print(f'{bcolors.RED}> Warning{bcolors.RESET}: mp_threshold is set to {self.args.mp_threshold}, a positive value smaller than the multiplicity {self.args.multiplicity}.')
            print(f'{bcolors.BLUE}>{bcolors.RESET} {bcolors.ORANGE}ALLEGRO{bcolors.RESET} cannot remove all but {self.args.mp_threshold} guides from each record and still ensure each record is targeted at least {self.args.multiplicity} times.')
            print(f'{bcolors.BLUE}>{bcolors.RESET} Re-adjusting mp_threshold to be equal to multiplicity {self.args.multiplicity}. You may also set mp_threshold to 0 to disable this feature. Refer to the manual for more details.')
            self.args.mp_threshold = self.args.multiplicity

        # ------------------------------------------------------------------------------
        #   beta
        # ------------------------------------------------------------------------------
        self.args.beta = int(self.args.beta)

        if self.args.beta <= 0:
            print(f'{bcolors.BLUE}>{bcolors.RESET} Beta is disabled. {bcolors.ORANGE}ALLEGRO{bcolors.RESET} will minimize the set size.')
            self.args.beta = 0

            if self.args.scorer != 'dummy':
                print(f'{bcolors.BLUE}>{bcolors.RESET} Scorer is set to {self.args.scorer} and beta is disabled. {bcolors.ORANGE}ALLEGRO{bcolors.RESET} will score the guides for information only and will not use them in calculations.')
        
        if self.args.scorer == 'dummy' and self.args.beta > 0:
            # No feasible solutions if there are fewer guides than beta
            # Say there are 5 species, 5 guides total, and beta is set to 1. Say that none of the species share any guides.
            # This will ask ALLEGRO to find 1 guide out of 5 to cover all 5 species. There is no solution.
            if self.args.track == 'track_e':
                print(f'{bcolors.BLUE}>{bcolors.RESET} The scorer is set to dummy and beta to {self.args.beta}. {bcolors.ORANGE}ALLEGRO{bcolors.RESET} will try to find {self.args.beta} guides to cover all genes.')
            if self.args.track == 'track_a':
                print(f'{bcolors.BLUE}>{bcolors.RESET} The scorer is set to dummy and beta to {self.args.beta}. {bcolors.ORANGE}ALLEGRO{bcolors.RESET} will try to find {self.args.beta} guides to cover all species.')
            
            if self.args.enable_solver_diagnostics:
                print(f'{bcolors.BLUE}>{bcolors.RESET} {bcolors.ORANGE}ALLEGRO{bcolors.RESET} may find that there is no feasible solution if the number of shared guides is fewer than beta in which case it will find the smallest beta for you.')
            else:
                print(f'{bcolors.RED}> Warning{bcolors.RESET}: {bcolors.ORANGE}ALLEGRO{bcolors.RESET} may find that there is no feasible solution if the number of shared guides is fewer than beta. Enabling solver diagnostics in {self.args.config} will find the smallest beta for you.')

        if self.args.beta > 0 and self.args.beta < self.args.multiplicity:
            print(f'{bcolors.RED}> Warning{bcolors.RESET}: Beta is set to {self.args.beta}, a positive value smaller than the multiplicity {self.args.multiplicity}')
            print(f'{bcolors.BLUE}>{bcolors.RESET} {bcolors.ORANGE}ALLEGRO{bcolors.RESET} cannot find a total of {self.args.beta} guides while each guide container is required to be targeted at least {self.args.multiplicity} times.')
            print(f'{bcolors.BLUE}>{bcolors.RESET} Re-adjusting beta to be equal to the multiplicity {self.args.multiplicity}. You may also set beta to 0. Refer to the manual for more details.')
            self.args.beta = self.args.multiplicity

        # ------------------------------------------------------------------------------
        #   early_stopping_patience
        # ------------------------------------------------------------------------------
        self.args.early_stopping_patience = int(self.args.early_stopping_patience)

        if self.args.early_stopping_patience < 10:
            print(f'{bcolors.RED}> Warning{bcolors.RESET}: early_stopping_patience is {self.args.early_stopping_patience} < 10 seconds. Auto adjusting to 10. ALLEGRO may be able to find a smaller sized solution with a higher patience.')
            self.args.early_stopping_patience = 10
        
        # ------------------------------------------------------------------------------
        #   protopacer length
        # ------------------------------------------------------------------------------
        if self.args.protospacer_length < 1:
            print(f'{bcolors.RED}> Error{bcolors.RESET}: Protospacer length is too short (< 1). Exiting.')
            sys.exit(1)

        if self.args.scorer == "ucrispr":
            if self.args.protospacer_length != 20 or self.args.pam != "NGG":
                print(f'{bcolors.RED}> Error{bcolors.RESET}: uCRISPR scorer only works with protospacer_length of 20 and the NGG PAM. Use the dummy scorer for guides other than Cas9. Exiting.')
                sys.exit(1)

        # ------------------------------------------------------------------------------

        # Create the output folder using the output directory and experiment name
        self.args.output_directory = self.create_output_directory(self.args.output_directory, self.args.experiment_name)

        return self.args

    def log_args(self) -> None:
        output_txt_path = os.path.join(self.args.output_directory, self.args.experiment_name + '_config_used.txt')

        with open(output_txt_path, 'w') as f:
            f.write(f'Config used for experiment {self.args.experiment_name}\n')
            for key, value in vars(self.args).items():
                if key != 'soundcheck':
                    f.writelines(f'{key}: {value}\n')

        self.output_txt_path = output_txt_path
        return output_txt_path

    def log_time(self, total_time_elapsed=0):
        end_time = time.time()

        # Calculate the elapsed time in seconds
        elapsed_seconds = end_time - self.start_time + total_time_elapsed

        # Convert elapsed_seconds to a timedelta object
        time_elapsed = timedelta(seconds=elapsed_seconds)

        # Extract the time components (hours, minutes, seconds)
        hours = time_elapsed.seconds // 3600
        minutes = (time_elapsed.seconds // 60) % 60
        seconds = time_elapsed.total_seconds() % 60 

        if hours == 0 and minutes == 0:
            print(f'{bcolors.BLUE}> {bcolors.RESET}{bcolors.ORANGE}ALLEGRO{bcolors.RESET} experiment took {seconds:.2f} seconds.')
        elif hours == 0:
            print(f'{bcolors.BLUE}> {bcolors.RESET}{bcolors.ORANGE}ALLEGRO{bcolors.RESET} experiment took {minutes} minutes, {seconds:.2f} seconds.')
        else:
            print(f'{bcolors.BLUE}> {bcolors.RESET}{bcolors.ORANGE}ALLEGRO{bcolors.RESET} experiment took {hours} hours, {minutes} minutes, {seconds:.2f} seconds.')
        
        with open(self.output_txt_path, 'a') as f:
            f.write(f'ALLEGRO experiment took {hours} hours, {minutes} minutes, {seconds:.2f} seconds.')

    def create_output_directory(self, output_directory: str, experiment_name: str) -> str:
        dir_name = os.path.join(output_directory, experiment_name)

        # Check if the directory already exists
        if os.path.exists(dir_name):
            # If it exists, append a number to the directory name
            i = 1

            while os.path.exists(f'{dir_name}_{i}'):
                i += 1

            dir_name = f'{dir_name}_{i}'

        print(f'{bcolors.BLUE}>{bcolors.RESET} Creating directory {dir_name}')
        os.makedirs(dir_name)

        return dir_name
