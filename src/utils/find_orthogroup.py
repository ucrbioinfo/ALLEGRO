import os
import re
import sys
import yaml
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def hash_gene_to_prot_id(
    cds_path: str, 
    gene_re=r'\[gene=(.*?)\]', 
    prot_id_re=r'\[protein_id=(.*?)\]'
    ) -> dict:

    hash = dict()
    cds = list(SeqIO.parse(open(cds_path), 'fasta'))

    for record in cds:
        gene_match = re.search(gene_re, record.description)

        if gene_match:
            prot_id = re.search(prot_id_re, record.description)
            
            if prot_id:
                hash[gene_match.group(1)] = prot_id.group(1)
            else:
                hash[gene_match.group(1)] = ''
    
    return hash


def hash_locus_tag_to_prot_id(
    cds_path: str, 
    tag_re=r'\[locus_tag=(.*?)\]', 
    prot_id_re=r'\[protein_id=(.*?)\]'
    ) -> dict:

    hash = dict()
    cds = list(SeqIO.parse(open(cds_path), 'fasta'))

    for record in cds:
        tag_match = re.search(tag_re, record.description)

        if tag_match:
            prot_id = re.search(prot_id_re, record.description)
            
            if prot_id:
                hash[tag_match.group(1)] = prot_id.group(1)
            else:
                hash[tag_match.group(1)] = ''
    
    return hash


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--config',
        type=argparse.FileType(mode='r'),
        default='src/utils/find_orthogroup_config.yaml', 
        help='The config file to use. Must be placed in the root folder.',
    )

    args = parser.parse_args()
    if args.config:
        data = yaml.load(args.config, Loader=yaml.FullLoader)
        arg_dict = args.__dict__

        for key, value in data.items():
            arg_dict[key] = value

    return args


def main() -> int:
    args = parse_arguments()

    output_directory = os.path.join(args.cds_directory, 'orthogroups', '')
    ortho_path = 'data/input/orthogroups/Orthogroups.tsv'

    print('Reading Orthogroups.tsv from {path}'.format(
        path=ortho_path,
    ))
    df = pd.read_csv(ortho_path, sep='\t', low_memory=False)
    
    print('Read {n} rows.'.format(n=len(df)))
    # print('Only keeping rows where all species have an entry.')
    # shared = df.dropna()

    # print('Dropped {n} rows where not all species had an orthologous entry.'.format(
    #     n=len(df)-len(shared),
    # ))
    shared = df.copy()
    species_names = shared.columns[1:]
    n_species = len(species_names)

    print('{n} species loaded.'.format(n=n_species))
    print('Species names {name}'.format(name=species_names.values.tolist()))

    file_paths = [args.cds_directory + name for name in os.listdir(args.cds_directory) if name.split('_')[0] in species_names]

    path_dict = dict()
    for name in species_names:
        for path in file_paths:
            if name in path:
                path_dict[name] = path

    print('File paths for each species:')
    for name, path in path_dict.items():
        print(name + ':', path)
    print()

    name_to_prot_id = dict()
    prot_id_to_name = dict()

    gene_to_prot_id_hash = hash_gene_to_prot_id(path_dict[args.species])
    for auxo in args.auxotroph_names:
        if auxo in gene_to_prot_id_hash:
            id = gene_to_prot_id_hash[auxo]

            if len(shared[shared[args.species] == id]) > 0:
                name_to_prot_id[auxo] = id
                prot_id_to_name[id] = auxo
            else:
                print('Auxo gene', auxo, 'has no orthogroup across all', n_species, 'species.')
                print(df[df[args.species] == id])

        else:
            print('Auxo gene name', auxo, 'not found in', path_dict[args.species])


    tag_to_prot_id_hash = hash_locus_tag_to_prot_id(path_dict[args.species])
    for auxo in args.auxotroph_tags:
        if auxo in tag_to_prot_id_hash:
            id = tag_to_prot_id_hash[auxo]

            if len(shared[shared[args.species] == id]) > 0:
                name_to_prot_id[auxo] = id
                prot_id_to_name[id] = auxo
            else:
                print('Auxo tag', auxo, 'has no orthogroup across all', n_species, '.')
                print(df[df[args.species] == id])
        else:
            print('Auxo tag', auxo, 'not found in', path_dict[args.species])


    print('{n} out of {k} input genes (auxotroph_name + auxotroph_tags) for {species} have orthogroups across all {all} species.'.format(
        n=len(name_to_prot_id),
        k=len(args.auxotroph_names) + len(args.auxotroph_tags),
        species=args.species,
        all=len(species_names)
        )
    )

    prot_id_to_orthogroup = dict()
    orthogroup_to_prot_id = dict()

    for id in name_to_prot_id.values():
        group = shared[shared[args.species] == id]['Orthogroup'].values[0]
        prot_id_to_orthogroup[id] = group
        orthogroup_to_prot_id[group] = id

    for gene_name in args.gene_names:
        print('Selecting for gene {name}'.format(name=gene_name))
        print('Associated protein_id is {id}'.format(id=name_to_prot_id[gene_name]))
        print('Found in orthogroup {group}'.format(group=prot_id_to_orthogroup[name_to_prot_id[gene_name]]))
        print(df[df[args.species] == name_to_prot_id[gene_name]])
        print()

    name_to_selected_ids = {name: name_to_prot_id[name] for name in args.gene_names}

    species_name_to_prot_ids = shared[shared[args.species].isin(name_to_selected_ids.values())].iloc[:, 1:].to_dict(orient='list')  # type: ignore

    dir = os.listdir(output_directory)
    # Checking if the list is empty or not
    if len(dir) != 0:
        answer = input('There are already files in {dir}. REMOVE fasta files? (yes/no): '.format(dir=output_directory)) 
        if answer == 'yes' or answer == 'y':
            for f in dir:
                if any([f.endswith('.fna'), f.endswith('.fa'), f.endswith('.fasta')]):
                    os.remove(os.path.join(output_directory, f))

    found_id_for_these = list()
    for name in species_names:
        prot_ids = species_name_to_prot_ids[name]
        cds_path = path_dict[name]

        records = list(SeqIO.parse(open(cds_path), 'fasta'))

        for record in records:
            prot_id_in_record = re.search(args.prot_id_regex, record.description)

            if prot_id_in_record:
                if prot_id_in_record.group(1) in prot_ids:
                    # print('Found prot_id {id} for species {name} in {path}'.format(
                    #     id=prot_id_in_record.group(1), 
                    #     name=name, path=cds_path,
                    #     )
                    # )

                    sequence = SeqRecord(
                        record.seq, 
                        id=prot_id_in_record.group(1), 
                        description=record.description,
                    )
                    output_path = os.path.join(output_directory, name + '_cds.fna')

                    with open(output_path, 'a') as f:
                        SeqIO.write(sequence, f, 'fasta')

                    found_id_for_these.append(name)

    print()
    print('Found genes in the following {n}/{m} species:'.format(
        n=len(set(found_id_for_these)),
        m=len(species_names),
    ))

    print(found_id_for_these)

    print('Missing:')
    print(set(species_names) - set(found_id_for_these))
    print('You should remove these from input_species.csv')

    return 0

if __name__ == '__main__':
    sys.exit(main())