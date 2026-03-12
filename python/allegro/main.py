import sys

from allegro.utils.shell_colors import bcolors
from allegro.utils.guide_finder import GuideFinder
from allegro.utils.configurator import Configurator
import allegro.utils.postprocessing as postprocessing
import allegro.utils.write_solution_to_file as solution_writer
from allegro.scorers.scorer_factory import ScorerFactory

try:
    import setproctitle
    setproctitle.setproctitle("ALLEGRO")
except ModuleNotFoundError:
    pass

try:
    from allegro._kirschtorte import KirschtorteCython
except ImportError as e:
    print(f"{bcolors.RED}>{bcolors.RESET} Dev error. Failed to load native extension: {e}")
    raise

def main() -> int:
    config = Configurator()
    if not any(x in sys.argv for x in ("-h", "--help", "--soundcheck")):
        config.begruessung()  # Guten Tag
    
    # Read command line arguments and config.yaml.
    # Some arguments are not shown in config.yaml and are defaulted in parse_configurations().
    # Arguments specified on the command line have priority and will replace those in config.yaml.
    config.parse_configurations()

    # This is the only function that changes the arguments in the config file.
    # Warnings and info to help with user error. Also hard-codes some paths and arguments.
    args = config.check_and_fix_configurations()

    # Write the current configuration to a log file in the output folder.
    config.log_args()

    # Instantiate the appropriate scorer.
    scorer_factory = ScorerFactory()
    guide_scorer = scorer_factory.get_scorer(scorer_name=args.scorer)

    # Initialize which guides should be kept from container sequences.
    GuideFinder.initialize(
        pam_list=[args.pam],
        blocklist_path=config.blocklist_path,
        protospacer_length=args.protospacer_length,
        filter_by_gc=args.filter_by_gc,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        patterns_to_exclude=args.patterns_to_exclude,
        guide_scorer=guide_scorer,
        guide_score_threshold=args.guide_score_threshold)

    cython_engine = KirschtorteCython(
        beta=args.beta,
        track=args.track,
        multiplicity=args.multiplicity,
        monophonic_threshold=args.mp_threshold,
        preclustering=args.preclustering,
        seed_length=args.seed_region_is_n_upstream_of_pam,
        mismatched_allowed_after_seed=args.mismatches_allowed_after_seed_region,
        early_stopping_patience=args.early_stopping_patience,
        input_directory=args.input_directory,
        output_directory=args.output_directory,
        input_species_csv_file_path=args.input_species_path,
        input_species_path_column=args.input_species_path_column,
        enable_solver_diagnostics=args.enable_solver_diagnostics)

    solution_writer.output_solution_to_text(
        species_names=cython_engine.species_names,
        solution=cython_engine.solution,
        experiment_name=args.experiment_name,
        output_directory=args.output_directory)

    solution_writer.output_solution_to_library(
        solution=cython_engine.solution,
        experiment_name=args.experiment_name,
        output_directory=args.output_directory)

    if args.align_solution_to_input:
        solution_writer.align_solution_to_input_bowtie(
            pam=args.pam,
            solution=cython_engine.solution,
            experiment_name=args.experiment_name,
            input_csv_path=args.input_species_path,
            input_directory=args.input_directory,
            output_directory=args.output_directory,
            paths_csv_column_name=args.input_species_path_column,
            species_names_csv_column_name=config.species_column_name)

    if args.postclustering:
        postprocessing.cluster_solution(
            output_directory=args.output_directory,
            experiment_name=args.experiment_name,
            req_match_len=args.seed_region_is_n_upstream_of_pam,
            mm_allowed=args.mismatches_allowed_after_seed_region)

    if args.output_offtargets:
        postprocessing.report_offtargets(
            input_species_path=args.input_species_path,
            output_directory=args.output_directory,
            experiment_name=args.experiment_name,
            input_species_offtarget_dir=args.input_species_offtarget_dir,
            input_species_offtarget_column=args.input_species_offtarget_column,
            num_mismatches=args.report_up_to_n_mismatches,
            seed_region_is_n_upstream_of_pam=args.seed_region_is_n_upstream_of_pam,
            protospacer_length=args.protospacer_length,
            pam=args.pam)
    
    config.log_time()

    return 0

if __name__ == '__main__':
    sys.exit(main())
