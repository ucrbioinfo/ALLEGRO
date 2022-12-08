from cover_set_parsers.coversets import Coversets


class SolverBase:
    def __init__(self, coverset_parser: Coversets, alpha:float, beta: float, solver_engine: str) -> None:
        raise NotImplementedError


    def solve(self) -> set[str]:
        raise NotImplementedError