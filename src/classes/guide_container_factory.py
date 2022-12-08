from __future__ import annotations
import typing

if typing.TYPE_CHECKING:
    from classes.species import Species

import re
from Bio import SeqIO

from classes.gene import Gene
from classes.chromosome import Chromosome
from classes.guide_container import GuideContainer


class GuideContainerFactory:
    def __init__(self, guide_source: str) -> None:
        self.guide_source = guide_source


    def make_guide_containers(self, species_object: Species) -> list[GuideContainer]:
        guide_container_list: list[GuideContainer] = list()

        cds_path = species_object.get_cds_path()
        genome_path = species_object.get_genome_path()
        guide_scorer_obj = species_object.get_scorer()

        match self.guide_source:
            case 'from_orthogroups':
                records = list(SeqIO.parse(open(cds_path), 'fasta'))

                gene_regex = r'\[gene=(.*?)\]'
                tag_regex = r'\[locus_tag=(.*?)\]'
                protein_id_regex = r'\[protein_id=(.*?)\]'

                for id, cds_record in enumerate(records):
                    tag_match = re.search(tag_regex, cds_record.description)
                    locus_tag = tag_match.group(1) if tag_match is not None else 'Locus tag not available.'

                    gene_match = re.search(gene_regex, cds_record.description)
                    gene_name = gene_match.group(1) if gene_match is not None else 'Gene name not available.'

                    protein_id_match = re.search(protein_id_regex, cds_record.description)
                    protein_id = protein_id_match.group(1) if protein_id_match is not None else 'Protein ID not available.'

                    guide_container_list.append(Gene(
                        species=species_object,
                        guide_scorer=guide_scorer_obj,
                        sequence=str(cds_record.seq),
                        gene_name=gene_name,
                        locus_tag=locus_tag,
                        protein_id=protein_id, 
                        string_id=cds_record.id,
                        integer_id=id,
                        )
                    )

            case 'from_genome':
                records = list(SeqIO.parse(open(genome_path), 'fasta'))
                
                for id, chromosome_record in enumerate(records):
                    guide_container_list.append(Chromosome(
                        species=species_object,
                        guide_scorer=guide_scorer_obj,
                        sequence=str(chromosome_record.seq).upper(),
                        string_id=chromosome_record.id,
                        integer_id=id,
                        )
                    )
                    
            case _:
                print('No such source as {source}. Check config.yaml'.format(source=self.guide_source))
                raise NotImplementedError

        return guide_container_list