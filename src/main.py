import sys

from utils.configurator import Configurator
import utils.postprocessing as postprocessing
import utils.write_solution_to_file as write_solution
from cython_libs.kirschtorte import KirschtorteCython as coverset  # type: ignore


def main() -> int:
    configurator = Configurator()
    configurator.begruessung()  # Guten Tag!
    
    # Read command line arguments and config.yaml.
    # Some arguments are not shown in config.yaml and are defaulted in parse_configurations().
    # Arguments specified on the command line have priority over and will replace those in config.yaml.
    args = configurator.parse_configurations()

    # This is the only function that changes the arguments in the config file.
    # Warnings and info to help with user error. Also hard-codes some paths and arguments.
    args, scorer_settings = configurator.check_and_fix_configurations()
    
    # Write the current configuration to a log file in the output folder.
    configurator.log_args()

    coversets_obj = coverset(
        beta=args.beta,
        track=args.track,
        early_stopping_patience=args.early_stopping_patience,
        scorer_name=args.scorer,
        cas_variant='cas9',
        guide_length=20,
        scorer_settings=scorer_settings,
        input_directory=args.input_directory,
        output_directory=args.output_directory,
        input_species_csv_file_path=args.input_species_path,
        input_species_path_column=args.input_species_path_column,
        cut_multiplicity=args.multiplicity,
        monophonic_threshold=args.mp_threshold
    )

    # TODO Use Bowtie to align these. Detect changes in last modified of the CDS files and rerun bowtie-build if file has been changed.
    output_csv_path = write_solution.write_solution_to_file_bowtie(
        pam='NGG',
        species_names=coversets_obj.species_names,
        solution=coversets_obj.solution,
        experiment_name=args.experiment_name,
        input_csv_path=args.input_species_path,
        input_directory=args.input_directory,
        paths_csv_column_name=args.input_species_path_column,
        species_names_csv_column_name='species_name',
        output_directory=args.output_directory
    )

    if args.cluster_guides:
        postprocessing.cluster_solution(
            solution_path=output_csv_path,
            req_match_len=args.seed_region_is_n_upstream_of_pam,
            mm_allowed=args.mismatches_allowed_after_seed_region
        )

    if args.output_offtargets:
        postprocessing.report_offtargets(
            input_species_path=args.input_species_path,
            solution_path=output_csv_path,
            output_dir=args.output_directory,
            input_species_offtarget_dir=args.input_species_offtarget_dir,
            input_species_offtarget_column=args.input_species_offtarget_column,
            experiment_name=args.experiment_name,
            num_mismatches=args.report_up_to_n_mismatches,
            seed_region_is_n_upstream_of_pam=args.seed_region_is_n_upstream_of_pam,
            pam_length=3
        )
        
    configurator.log_time()

    return 0


if __name__ == '__main__':
    sys.exit(main())