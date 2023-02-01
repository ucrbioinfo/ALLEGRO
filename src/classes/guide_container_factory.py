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
    def __init__(self) -> None:
        pass

    def make_guide_containers(self, species_object: Species) -> list[GuideContainer]:
        guide_container_list: list[GuideContainer] = list()

        cds_path = species_object.cds_path
        genome_path = species_object.genome_path
        guide_source = species_object.guide_source
        guide_scorer_obj = species_object.guide_scorer

        match guide_source:
            case 'from_orthogroups' | 'from_all_cds':
                records = list(SeqIO.parse(open(cds_path), 'fasta'))

                gene_regex = r'\[gene=(.*?)\]'
                tag_regex = r'\[locus_tag=(.*?)\]'
                protein_id_regex = r'\[protein_id=(.*?)\]'
                reference_species_regex = r'\[ref_species=(.*?)\]'
                orthologous_name_regex = r'\[orthologous_to_gene=(.*?)\]'
                orthologous_protein_regex = r'\[orthologous_to_ref_protein=(.*?)\]'

                for id, cds_record in enumerate(records):
                    # Locus tag example: [locus_tag=KLMA_50610], extracts KLMA_50610
                    tag_match = re.search(tag_regex, cds_record.description)
                    locus_tag = tag_match.group(1) if tag_match is not None else 'N/A'

                    # This gene's own name -- Usually N/A for unannotated CDS files
                    gene_match = re.search(gene_regex, cds_record.description)
                    gene_name = gene_match.group(1) if gene_match is not None else 'N/A'

                    # This gene's own protein id
                    protein_id_match = re.search(protein_id_regex, cds_record.description)
                    protein_id = protein_id_match.group(1) if protein_id_match is not None else 'N/A'

                    # For example, [orthologous_to_ref_protein=XP_022674739.1], extracts XP_022674739.1
                    ortho_prot_to_match = re.search(orthologous_protein_regex, cds_record.description)
                    ortho_prot_id = ortho_prot_to_match.group(1) if ortho_prot_to_match is not None else 'N/A'

                    # For example, [orthologous_to_gene=HIS7], extracts HIS7
                    ortho_gene_to_match = re.search(orthologous_name_regex, cds_record.description)
                    ortho_gene_name = ortho_gene_to_match.group(1) if ortho_gene_to_match is not None else 'N/A'

                    # For example, [ref_species=kluyveromyces_marxianus], extracts kluyveromyces_marxianus
                    reference_species_match = re.search(reference_species_regex, cds_record.description)
                    ref_species = reference_species_match.group(1) if reference_species_match is not None else 'N/A'

                    guide_container_list.append(Gene(
                        integer_id=id,
                        gene_name=gene_name,
                        locus_tag=locus_tag,
                        protein_id=protein_id,
                        species=species_object,
                        string_id=cds_record.id,
                        ref_species=ref_species,
                        sequence=str(cds_record.seq),
                        guide_scorer=guide_scorer_obj,
                        orthologous_to_prot=ortho_prot_id,
                        orthologous_to_gene=ortho_gene_name,
                        )
                    )

            case 'from_genome':
                records = list(SeqIO.parse(open(genome_path), 'fasta'))
                
                for id, chromosome_record in enumerate(records):
                    guide_container_list.append(Chromosome(
                        integer_id=id,
                        species=species_object,
                        guide_scorer=guide_scorer_obj,
                        string_id=chromosome_record.id,
                        sequence=str(chromosome_record.seq).upper(),
                        )
                    )
                    
            case _:
                print('No such source as {source}. Check config.yaml'.format(source=self.guide_source))
                raise ValueError

        return guide_container_list