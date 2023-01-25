from classes.guide_container import GuideContainer


class Scorer:
    def __init__(self) -> None:
        pass
    

    def score_sequence(self, guide_container: GuideContainer) -> list[tuple[str, str, float, int]]:
        raise NotImplementedError