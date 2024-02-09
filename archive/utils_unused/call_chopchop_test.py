import chopchop_wrapper

species_name = 'hpolymorpha'
absolute_path_to_chopchop = '/home/amohs002/projects/research/chopchop/chopchop.py'
absolute_path_to_bowtie_build = '/home/amohs002/projects/research/chopchop/bowtie/bowtie-build'
absolute_path_to_species_genome = '/home/amohs002/projects/research/application_v2/data/input/genomes/hpolymorpha_genomic.fna'


chopchop_wrapper.make_species_bowtie_index(
    species_name=species_name,
    absolute_path_to_bowtie_build=absolute_path_to_bowtie_build,
    absolute_path_to_species_genome=absolute_path_to_species_genome,
)


chopchop_wrapper.run_chopchop_for_species_genome(
    species_name=species_name,
    absolute_path_to_chopchop=absolute_path_to_chopchop,
    absolute_path_to_species_genome=absolute_path_to_species_genome,
    scoring_method='DOENCH_2016'
)