import re
import sys
from Bio import SeqIO

from scorers.scorer_base import Scorer
from classes.guide_container import GuideContainer
from utils.shell_colors import bcolors


class GuideContainerFactory:
    def __init__(self) -> None:
        pass


    def make_guide_containers(
        self,
        species_name: str,
        records_path: str,
        guide_scorer_obj: Scorer
        ) -> list[GuideContainer]:
        
        guide_container_list: list[GuideContainer] = list()

        records = list()

        try:
            records = list(SeqIO.parse(open(records_path), 'fasta'))
        except FileNotFoundError:
            print(f'{bcolors.RED}> Error{bcolors.RESET}: Cannot find the fasta file {records_path} specified in the input csv for this species. Is the specified genome/CDS file for this species available? Exiting.')
            sys.exit(1)

        gene_regex = r'\[gene=(.*?)\]'
        tag_regex = r'\[locus_tag=(.*?)\]'
        protein_id_regex = r'\[protein_id=(.*?)\]'
        reference_species_regex = r'\[ref_species=(.*?)\]'
        orthologous_name_regex = r'\[orthologous_to_gene=(.*?)\]'
        orthologous_protein_regex = r'\[orthologous_to_ref_protein=(.*?)\]'

        for record in records:
            # Locus tag example: [locus_tag=KLMA_50610], extracts KLMA_50610
            tag_match = re.search(tag_regex, record.description)
            locus_tag = tag_match.group(1) if tag_match is not None else 'N/A'

            # This gene's own name -- Usually N/A for unannotated CDS files or chromosomes/scaffolds
            gene_match = re.search(gene_regex, record.description)
            gene_name = gene_match.group(1) if gene_match is not None else 'N/A'

            # This gene's own protein id e.g., XP_022674739.1
            protein_id_match = re.search(protein_id_regex, record.description)
            protein_id = protein_id_match.group(1) if protein_id_match is not None else 'N/A'

            # For example, [orthologous_to_ref_protein=XP_022674739.1], extracts XP_022674739.1
            ortho_prot_to_match = re.search(orthologous_protein_regex, record.description)
            ortho_prot_id = ortho_prot_to_match.group(1) if ortho_prot_to_match is not None else 'N/A'

            # For example, [orthologous_to_gene=HIS7], extracts HIS7
            ortho_gene_to_match = re.search(orthologous_name_regex, record.description)
            ortho_gene_name = ortho_gene_to_match.group(1) if ortho_gene_to_match is not None else 'N/A'

            # For example, [ref_species=kluyveromyces_marxianus], extracts kluyveromyces_marxianus
            reference_species_match = re.search(reference_species_regex, record.description)
            ref_species = reference_species_match.group(1) if reference_species_match is not None else 'N/A'

            guide_container_list.append(GuideContainer(
                gene_name=gene_name,
                locus_tag=locus_tag,
                protein_id=protein_id,
                species_name=species_name,
                string_id=record.id,
                ref_species=ref_species,
                sequence=str(record.seq).upper(),
                guide_scorer=guide_scorer_obj,
                orthologous_to_prot=ortho_prot_id,
                orthologous_to_gene=ortho_gene_name,
                )
            )

        return guide_container_list