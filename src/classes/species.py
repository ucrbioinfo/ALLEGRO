import os

from classes.guide import Guide
from classes.guide_container import GuideContainer
from classes.guide_container_factory import GuideContainerFactory
from scorers.scorer_base import Scorer


class Species:
    def __init__(
        self,
        guide_scorer: Scorer,
        guide_container_factory: GuideContainerFactory,
        name: str, 
        id: int, 
        genome_path: str,
        cds_path: str,
        ) -> None:
        
        self.guide_scorer = guide_scorer
        self.guide_container_factory = guide_container_factory
        self.name = name
        self.id = id
        self.genome_path = genome_path
        self.cds_path = cds_path

        self.guide_containers: list[GuideContainer] = list()


    def get_cas9_guides(self) -> list[Guide]:
        self.guide_containers = self.guide_container_factory.make_guide_containers(self)

        cas9_guides_list: list[Guide] = list()
        for container in self.guide_containers:
            cas9_guides_list.extend(container.get_cas9_guides())

        return cas9_guides_list
    
    
    def get_scorer(self) -> Scorer:
        return self.guide_scorer

    
    def get_genome_path(self) -> str:
        return self.genome_path

    
    def get_cds_path(self) -> str:
        return self.cds_path


    def get_name(self) -> str:
        return self.name


    def get_id(self) -> int:
        return self.id

    
    def get_guide_containers(self) -> list[GuideContainer]:
        return self.guide_containers