import os


def make_species_bowtie_index( 
    species_name: str, 
    absolute_path_to_bowtie_build: str,
    absolute_path_to_species_genome: str) -> None:
    output_directory = '../../data/output/bowtie_indices/'

    print('Running bowtie-build for', species_name)
    # conda run -n chopchop path_to/bowtie/bowtie-build ../../data/input/genomes/hpolymorpha_genomic.fna hpoly
    os.system('conda run -n chopchop' + ' ' + absolute_path_to_bowtie_build + ' ' + absolute_path_to_species_genome + ' ' + species_name + ' ' + '-q')
    os.system('mv *.ebwt' + ' ' + output_directory)


def run_chopchop_for_species_genome(
    species_name: str,
    absolute_path_to_chopchop: str,
    absolute_path_to_species_genome: str,
    scoring_method: str = 'DOENCH_2016') -> None:
    
    output_path = os.path.join('../../data/output/chopchop_scores/', species_name, 'output.csv')
    output_directory = os.path.join('../../data/output/chopchop_scores/', species_name)

    if os.path.exists(output_path):
        print('Path', output_path, 'already exists for species', species_name + '. Aborting running ChopChop.')
    else:
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        print('Running ChopChop for', species_name)

        os.system('conda run -n chopchop ' + 
        absolute_path_to_chopchop + ' -F ' + 
        ' -Target ' + absolute_path_to_species_genome + 
        ' -o ' + output_directory + 
        ' -G ' + species_name + 
        ' --scoringMethod ' + scoring_method + 
        ' > ' + output_path)

        os.system('rm -rf ' + output_directory + ' *.offtargets')