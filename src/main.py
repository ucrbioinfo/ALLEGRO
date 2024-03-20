import sys

from utils.configurator import Configurator
import utils.postprocessing as postprocessing
import utils.write_solution_to_file as solution_writer
from cython_libs.kirschtorte import KirschtorteCython as coverset  # type: ignore
from cython_libs.kirschtorte import EinfacherModusCython as einfacher_modus  # type: ignore


def main() -> int:
    configurator = Configurator()
    configurator.begruessung()  # Guten Tag
    
    # Read command line arguments and config.yaml.
    # Some arguments are not shown in config.yaml and are defaulted in parse_configurations().
    # Arguments specified on the command line have priority and will replace those in config.yaml.
    args = configurator.parse_configurations()

    # This is the only function that changes the arguments in the config file.
    # Warnings and info to help with user error. Also hard-codes some paths and arguments.
    args, scorer_settings = configurator.check_and_fix_configurations()
    
    # Write the current configuration to a log file in the output folder.
    configurator.log_args()

    if configurator.using_easy_mode:
        coversets_obj = einfacher_modus(
            beta=args.beta,
            track=args.track,
            early_stopping_patience=args.early_stopping_patience,
            cas_variant=args.cas,
            guide_length=20,
            output_directory=args.output_directory,
            multiplicity=args.multiplicity,
            monophonic_threshold=args.mp_threshold,
            input_csv_path_with_guides=args.input_csv_path_with_guides,
            enable_solver_diagnostics=args.enable_solver_diagnostics
        )

        solution_writer.output_solution_to_text(
        species_names=coversets_obj.species_names,
        solution=coversets_obj.solution,
        experiment_name=args.experiment_name,
        output_directory=args.output_directory
        )

        solution_writer.output_solution_to_library(
            solution=coversets_obj.solution,
            experiment_name=args.experiment_name,
            output_directory=args.output_directory
        )

        solution_writer.cross_reference_solution_to_input(
            solution=coversets_obj.solution,
            experiment_name=args.experiment_name,
            input_csv_path=args.input_csv_path_with_guides,
            output_directory=args.output_directory
        )
    else:
        coversets_obj = coverset(
            beta=args.beta,
            track=args.track,
            precluster=args.preclustering,
            seed_length=args.seed_region_is_n_upstream_of_pam,
            mismatched_allowed_after_seed=args.mismatches_allowed_after_seed_region,
            early_stopping_patience=args.early_stopping_patience,
            scorer_name=args.scorer,
            cas_variant=args.cas,
            guide_length=20,
            scorer_settings=scorer_settings,
            input_directory=args.input_directory,
            output_directory=args.output_directory,
            input_species_csv_file_path=args.input_species_path,
            input_species_path_column=args.input_species_path_column,
            multiplicity=args.multiplicity,
            monophonic_threshold=args.mp_threshold,
            enable_solver_diagnostics=args.enable_solver_diagnostics
        )

        solution_writer.output_solution_to_text(
            species_names=coversets_obj.species_names,
            solution=coversets_obj.solution,
            experiment_name=args.experiment_name,
            output_directory=args.output_directory
        )

        solution_writer.output_solution_to_library(
            solution=coversets_obj.solution,
            experiment_name=args.experiment_name,
            output_directory=args.output_directory
        )

        solution_writer.align_solution_to_input_bowtie(
            pam=args.pam,
            solution=coversets_obj.solution,
            experiment_name=args.experiment_name,
            input_csv_path=args.input_species_path,
            input_directory=args.input_directory,
            output_directory=args.output_directory,
            paths_csv_column_name=args.input_species_path_column,
            species_names_csv_column_name='species_name'
        )

    if args.postclustering:
        postprocessing.cluster_solution(
            output_directory=args.output_directory,
            experiment_name=args.experiment_name,
            req_match_len=args.seed_region_is_n_upstream_of_pam,
            mm_allowed=args.mismatches_allowed_after_seed_region
        )

    if args.output_offtargets:
        postprocessing.report_offtargets(
            input_species_path=args.input_species_path,
            output_directory=args.output_directory,
            experiment_name=args.experiment_name,
            input_species_offtarget_dir=args.input_species_offtarget_dir,
            input_species_offtarget_column=args.input_species_offtarget_column,
            num_mismatches=args.report_up_to_n_mismatches,
            seed_region_is_n_upstream_of_pam=args.seed_region_is_n_upstream_of_pam,
            pam_length=3
        )
    
    configurator.log_time()

    return 0


if __name__ == '__main__':
    sys.exit(main())