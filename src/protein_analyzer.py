#!/usr/bin/env python3
"""
Protein analyzer module.
Handles biophysical property calculations and analysis.
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from typing import Dict, Any, Tuple
from src.config import ANALYSIS


class ProteinAnalyzer:
    """Analyzes protein sequences and calculates biophysical properties."""

    def __init__(self, sequence: str):
        """
        Initialize analyzer with protein sequence.

        Args:
            sequence: Protein amino acid sequence
        """
        # Clean sequence
        self.original_sequence = sequence
        self.sequence = self._clean_sequence(sequence)
        self.analysis = ProteinAnalysis(self.sequence)

    def _clean_sequence(self, sequence: str) -> str:
        """
        Clean sequence by replacing non-standard amino acids.

        Args:
            sequence: Raw protein sequence

        Returns:
            Cleaned sequence
        """
        for old, new in ANALYSIS['clean_nonstandard_aa'].items():
            sequence = sequence.replace(old, new)
        return sequence

    def get_molecular_weight(self) -> float:
        """Get molecular weight in kDa."""
        return self.analysis.molecular_weight() / 1000

    def get_isoelectric_point(self) -> float:
        """Get isoelectric point (pI)."""
        return self.analysis.isoelectric_point()

    def get_gravy(self) -> float:
        """Get GRAVY score (Grand average of hydropathy)."""
        return self.analysis.gravy()

    def get_aromaticity(self) -> float:
        """Get aromaticity percentage."""
        return self.analysis.aromaticity() * 100

    def get_instability_index(self) -> float:
        """Get instability index."""
        return self.analysis.instability_index()

    def get_cysteine_count(self) -> int:
        """Get number of cysteines."""
        return self.analysis.count_amino_acids().get('C', 0)

    def get_extinction_coefficients(self) -> Tuple[int, int]:
        """
        Get extinction coefficients.

        Returns:
            Tuple of (reduced, oxidized) coefficients
        """
        coeff_reduced, coeff_oxidized = self.analysis.molar_extinction_coefficient()
        return coeff_reduced, coeff_oxidized

    def get_secondary_structure_fractions(self) -> Tuple[float, float, float]:
        """
        Get predicted secondary structure fractions.

        Returns:
            Tuple of (helix%, turn%, sheet%)
        """
        helix, turn, sheet = self.analysis.secondary_structure_fraction()
        return helix * 100, turn * 100, sheet * 100

    def get_all_properties(self) -> Dict[str, Any]:
        """
        Get all calculated biophysical properties.

        Returns:
            Dictionary with all properties
        """
        helix, turn, sheet = self.get_secondary_structure_fractions()
        coeff_red, coeff_ox = self.get_extinction_coefficients()

        return {
            'sequence_length': len(self.sequence),
            'molecular_weight_kda': self.get_molecular_weight(),
            'isoelectric_point': self.get_isoelectric_point(),
            'gravy_score': self.get_gravy(),
            'aromaticity_percent': self.get_aromaticity(),
            'instability_index': self.get_instability_index(),
            'cysteine_count': self.get_cysteine_count(),
            'extinction_coeff_reduced': coeff_red,
            'extinction_coeff_oxidized': coeff_ox,
            'helix_predicted': helix,
            'turn_predicted': turn,
            'sheet_predicted': sheet,
        }

    def is_stable(self, threshold: float = 40.0) -> bool:
        """
        Check if protein is stable based on instability index.

        Args:
            threshold: Instability index threshold (default: 40)

        Returns:
            True if stable, False if unstable
        """
        return self.get_instability_index() < threshold
