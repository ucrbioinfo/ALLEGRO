from coverset_parsers.coversets_db import CoversetsDB
from coverset_parsers.coversets_ram import CoversetsRAM
from coverset_parsers.coversets_base import CoversetsBase


class CoversetsFactory:
    __slots__ = []

    def __init__(self) -> None:
        pass


    def make_coversets_parser(self, storage_destination: str) -> CoversetsBase:
        # match storage_destination:
        #     case 'ram':
        pass
