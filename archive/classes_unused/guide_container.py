from classes.guide import Guide
from scorers.scorer_base import Scorer

from abc import ABC, abstractmethod


class GuideContainer(ABC):
    # the sequence of this guide container. This is a the CDS sequence
    # for a Gene object, and a chromosome sequence for a Chromosome object.
    sequence: str

    # return: the string id of this guide container.
    # E.g., NW_022983474.1 for some chromosome fasta entry.
    string_id: str
    species_name: str

    # The Guide scorer object assigned to this container.
    # See the 'src/scorers/' directory for a list.
    # Could be CHOPCHOP or DeepGuide or etc.
    # This option is affected by the setting in config.yaml
    guide_scorer: Scorer


    @abstractmethod
    def get_cas9_guides(self) -> list[Guide]:
        '''
        return: a list of Guide objects which have a score, sequence, and strand.
        Get these by calling self.guide_scorer.score_sequence(self) first. 
        '''


    @abstractmethod
    def get_attributes_dict(self) -> dict:
        '''
        return: the attributes dictionary for this container.
        For example,

        return dict({
            'genome_sequence': self.sequence,
            'genome_string_id': self.string_id,
            'genome_species_name': self.species.name,
        })
        '''