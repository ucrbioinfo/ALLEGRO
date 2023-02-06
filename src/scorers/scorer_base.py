from classes.guide_container import GuideContainer


class Scorer:
    def __init__(self, settings: dict[str, str]) -> None:
        pass
    

    def score_sequence(
        self,
        guide_container: GuideContainer,
        ) -> tuple[list[str], list[str], list[str], list[int], list[float]]:
        '''
        ## Args:
            * guide_container: Either a Gene or a Chromosome type guide container.
        
        ## Returns:
            A tuple of four lists:
            * The first list[str] is a list of the guides found in `sequence`.
            * The second list[str] is a list of guides in `sequence` with context around them.
            * The third list[str] is a list of 'F's and 'R's indicating on
              which strand, forward or reverse, each respective guide resides.
            * The fourth list[int] shows the locations of each guides in `sequence`.
            * The fifth list[float] indicates the efficiency scores of each guide.
        '''
        
        raise NotImplementedError