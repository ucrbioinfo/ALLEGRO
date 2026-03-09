from allegro.classes.guide import Guide
from allegro.classes.guide_container import GuideContainer
from allegro.classes.guide_container_factory import GuideContainerFactory

class Species:
    __slots__ = [
        'name',
        'records_path',
        'bowtie_index_path',
        'guide_objects_list',
        'guide_containers_list',
        'guide_container_factory',
    ]

    name: str
    records_path: str
    bowtie_index_path: str
    guide_objects_list: list[Guide]
    guide_containers_list: list[GuideContainer]
    guide_container_factory: GuideContainerFactory

    def __init__(
        self,
        name: str,
        records_path: str,
        guide_container_factory: GuideContainerFactory,
        ) -> None:
        
        self.name = name
        self.records_path = records_path
        self.guide_container_factory = guide_container_factory
        self.guide_containers_list: list[GuideContainer]
        self.guide_objects_list: list[Guide]

    def make_guide_containers(self) -> None:
        self.guide_containers_list = self.guide_container_factory.make_guide_containers(
            species_name=self.name,
            records_path=self.records_path)

    def get_guides_from_containers(self) -> list[Guide]:
        self.guide_objects_list: list[Guide] = list()

        container_idx = 0
        for container in self.guide_containers_list:
            self.guide_objects_list.extend(container.get_guides())

        return self.guide_objects_list
