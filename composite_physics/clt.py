"""Compatibility wrapper for the original ABD entry point."""
from .laminate import assemble_laminate_stiffness


def compute_ABD_matrix(layup: list, materials: list):
    return assemble_laminate_stiffness(layup, materials).ABD
