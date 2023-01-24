from classes.guide import Guide
from scorers.scorer_base import Scorer
from classes.guide_container import GuideContainer
from classes.guide_container_factory import GuideContainerFactory


class Species:
    def __init__(
        self,
        id: int, 
        name: str,
        cds_path: str,
        genome_path: str,
        guide_scorer: Scorer,
        guide_container_factory: GuideContainerFactory,
        ) -> None:
        
        self.id = id
        self.name = name
        self.cds_path = cds_path
        self.genome_path = genome_path
        self.guide_scorer = guide_scorer
        self.guide_container_factory = guide_container_factory

        self.guide_containers: list[GuideContainer] = list()


    def get_cas9_guides(self) -> list[Guide]:
        self.guide_containers = self.guide_container_factory.make_guide_containers(self)

        cas9_guides_list: list[Guide] = list()
        for container in self.guide_containers:
            cas9_guides_list.extend(container.get_cas9_guides())

        return cas9_guides_list