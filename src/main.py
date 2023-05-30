import sys

import utils.configurator as configurator
import utils.cluster_solution as cluster_solution
import utils.write_solution_to_file as write_solution
from cython_libs.rammingstone import RammingstoneCython as coverset  # type: ignore


def main() -> int:
    print('Welcome to ALLEGRO. All unspecified command-line arguments default to the values in config.yaml')
    
    args = configurator.parse_configurations()
    args = configurator.check_and_fix_configurations(args)
    args.output_directory = configurator.create_output_directory(args.output_directory, args.experiment_name)
    scorer_settings = configurator.configure_scorer_settings(args)
    mode_settings = configurator.configure_mode_settings(args)
    configurator.log_args(args)

    coversets_obj = coverset(
        beta=args.beta,
        monophonic_threshold=args.mp_threshold,
        cut_multiplicity=args.multiplicity,
        container_names=args.container_names,
        track=args.track,
        num_trials=args.num_trials,
        cas_variant=args.cas,
        guide_source=args.mode,
        guide_length=args.protospacer_length,
        scorer_name=args.scorer,
        scorer_settings=scorer_settings,
        output_directory=args.output_directory,
        input_cds_directory=args.input_cds_directory,
        input_species_csv_file_path=args.input_species_path,
        input_genome_directory=args.input_genomes_directory,
    )

    output_path: str = ''
    if args.mode == 'from_genome':
        output_path = write_solution.write_solution_to_file(
            species_names=coversets_obj.species_names,
            solution=coversets_obj.solution,
            experiment_name=args.experiment_name,
            input_csv_path=args.input_species_path,
            input_sequence_directory=mode_settings['input_sequence_directory'],
            paths_csv_column_name=mode_settings['paths_csv_column_name'],
            species_names_csv_column_name='species_name',
            output_directory=args.output_directory
        )
    elif args.mode == 'from_cds':
        output_path = write_solution.write_cds_solution_to_file(
            container_names=mode_settings['container_names'],
            solution=coversets_obj.solution,
            experiment_name=args.experiment_name,
            input_csv_path=args.input_species_path,
            input_sequence_directory=mode_settings['input_sequence_directory'],
            paths_csv_column_name=mode_settings['paths_csv_column_name'],
            species_names_csv_column_name='species_name',
            output_directory=args.output_directory
        )

    if args.cluster_guides:
        cluster_solution.cluster_solution(
            output_path,
            args.seed_region_is_n_from_pam,
            args.mismatches_allowed_after_seed_region)

    return 0


if __name__ == '__main__':
    sys.exit(main())