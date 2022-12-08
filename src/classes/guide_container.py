from __future__ import annotations
import typing

if typing.TYPE_CHECKING:
    from classes.species import Species

from classes.guide import Guide


class GuideContainer:
    def __init__(self) -> None:
        pass


    def get_species(self) -> Species:
        raise NotImplementedError

    
    def get_sequence(self) -> str:
        raise NotImplementedError

    
    def get_string_id(self) -> str:
        raise NotImplementedError

    
    def get_integer_id(self) -> int:
        raise NotImplementedError


    def get_cas9_guides(self) -> list[Guide]:
        raise NotImplementedError


    def get_guides(self) -> list[Guide]:
        raise NotImplementedError


    def get_attributes_dict(self) -> dict:
        raise NotImplementedError


    def get_gene_name(self) -> str:
        raise NotImplementedError