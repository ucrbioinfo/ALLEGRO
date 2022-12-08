from __future__ import annotations
import typing

if typing.TYPE_CHECKING:
    from classes.guide_container import GuideContainer
    from classes.species import Species


class Guide:
    def __init__(
        self,
        pam: str,
        endonuclease: str,
        score: float,
        strand: str,
        sequence: str,
        container: GuideContainer,
        sequence_with_context: str = '',
        ) -> None:

        self.pam = pam
        self.endonuclease = endonuclease
        self.score = score
        self.strand = strand
        self.sequence = sequence
        self.container = container
        self.sequence_with_context = sequence_with_context

    
    def get_pam(self) -> str:
        return self.pam

    
    def get_endonuclease(self) -> str:
        return self.endonuclease


    def get_score(self) -> float:
        return self.score


    def get_strand(self) -> str:
        return self.strand


    def get_sequence(self) -> str:
        return self.sequence


    def get_container(self) -> GuideContainer:
        return self.container


    def get_sequence_with_context(self) -> str:
        return self.sequence_with_context

    
    def get_attributes_dict(self) -> dict:
        return dict({
            'pam': self.pam, 
            'endonuclease': self.endonuclease,
            'score': self.score,
            'strand': self.strand,
            'sequence': self.sequence,
            'container': self.container,
            'sequence_with_context': self.sequence_with_context,
        })

    
    def get_species(self) -> Species:
        return self.container.get_species()


    def get_species_name(self) -> str:
        return self.container.get_species().get_name()


    def get_location(self) -> str:
        return self.container.get_string_id()


    def get_info(self) -> list[str]:
        info_list: list[str] = list()
        attributes_dict = self.get_attributes_dict()

        info_list.append('Endonuclease and PAM: ' + attributes_dict['endonuclease'] + ' ' + attributes_dict['pam'])
        info_list.append('Sequence: ' + attributes_dict['sequence'])
        info_list.append('Score: ' + str(attributes_dict['score']))
        info_list.append('Species: ' + attributes_dict['container'].get_species().get_name())
        info_list.append('Location: ' + attributes_dict['container'].get_string_id())
        info_list.append('Gene: '  + attributes_dict['container'].get_gene_name())
        info_list.append('Strand: ' + attributes_dict['strand'])

        return info_list


    def print_info(self) -> None:
        info_list = self.get_info()

        for info in info_list:
            print(info)


    def write_info_to_open_file(self, file_handle: typing.TextIO) -> None:
        info_list = self.get_info()

        for info in info_list:
            file_handle.write(info)
            file_handle.write('\n')