"""
formulae.py

A class containing formulae for each chapter in their respective objects, for use across all
problems in the book.

"""

from helpers.chapter_formulae.chapter_2 import AtomicAndNuclearPhysics
from helpers.chapter_formulae.chapter_3 import RadiationMatterInteraction


class Formulae:
    """Contains all equations in the book, separated by chapter"""

    ch2 = AtomicAndNuclearPhysics()
    ch3 = RadiationMatterInteraction()
