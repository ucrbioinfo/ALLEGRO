from classes.guide_container import GuideContainer
from scorers.scorer_base import Scorer


class DeepGuideWrapper(Scorer):
    def __init__(self) -> None:
        super().__init__()


    def score_sequence(self, guide_container: GuideContainer) -> list[tuple[str, str, float]]:
        return super().score_sequence(guide_container)