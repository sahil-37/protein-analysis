#!/usr/bin/env python3
"""
Unified Protein Enrichment Service

A complete automated pipeline that:
1. Fetches protein data from UniProt
2. Computes biophysical properties (pI, MW, GRAVY, instability, etc.)
3. Predicts structural/functional features (signal peptides, TM domains, disulfides, PTMs)
4. Outputs comprehensive JSON-structured protein profile

This is the one-stop enrichment system for protein analysis.

Usage:
    from protein_enrichment_service import enrich_protein, enrich_protein_batch

    # Single protein
    profile = enrich_protein('P04637')

    # Batch
    profiles = enrich_protein_batch(['P04637', 'P68871', 'P15018'])

    # Save to JSON
    profile.to_json('p53_enriched.json')
"""

import json
import sys
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict
from datetime import datetime
import requests
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import molecular_weight

from uniprot_client import UniProtClient
from data_processor import UniProtDataProcessor


# ============================================================================
# SECTION 1: Data Models
# ============================================================================

@dataclass
class GeneralInfo:
    """General protein metadata"""
    protein_name: str
    gene_name: Optional[str]
    organism: str
    organism_taxon_id: Optional[int]
    uniprot_id: str
    entry_type: str
    sequence_length: int
    reviewed: bool
    created_date: Optional[str]
    updated_date: Optional[str]


@dataclass
class SequenceData:
    """Sequence information"""
    fasta_sequence: str
    length: int
    crc64_checksum: Optional[str] = None


@dataclass
class BiophysicalProperties:
    """Computed sequence-based biophysical properties"""
    molecular_weight_da: float
    molecular_weight_kda: float
    isoelectric_point_pi: float
    aromaticity: float
    instability_index: float
    gravy: float  # Grand Average of Hydropathy
    aliphatic_index: float
    helix_fraction: float
    sheet_fraction: float
    coil_fraction: float
    cysteine_count: int


@dataclass
class ComputedCharges:
    """pH-dependent charge information"""
    ph_values: List[float]
    charges: List[float]
    isoelectric_point: float
    charge_at_ph_7_4: float


@dataclass
class SignalPeptidePrediction:
    """Signal peptide prediction"""
    present: bool
    cleavage_site: Optional[int] = None
    signal_peptide_seq: Optional[str] = None
    prediction_tool: str = "SignalP 6.0 (stub)"
    confidence: Optional[float] = None


@dataclass
class TransmembraneDomain:
    """Transmembrane domain prediction"""
    domain_id: int
    start: int
    end: int
    topology: str  # "inside-out" or "outside-in"


@dataclass
class TransmembranePrediction:
    """Transmembrane prediction results"""
    is_membrane_protein: bool
    tm_domains: List[TransmembraneDomain]
    prediction_tool: str = "TMHMM 2.0 (stub)"


@dataclass
class DisulfideBondPrediction:
    """Disulfide bond prediction"""
    cys_position_1: int
    cys_position_2: int
    confidence: float
    is_experimentally_confirmed: bool = False


@dataclass
class PTMPrediction:
    """Post-translational modification prediction"""
    ptm_type: str  # "phosphorylation", "glycosylation", "ubiquitination", etc.
    position: int
    residue: str
    confidence: float
    prediction_tool: str


@dataclass
class PredictedFeatures:
    """All structural/functional predictions"""
    signal_peptide: SignalPeptidePrediction
    transmembrane: TransmembranePrediction
    disulfide_bonds: List[DisulfideBondPrediction]
    ptm_predictions: List[PTMPrediction]


@dataclass
class ProteinEnrichmentProfile:
    """Complete protein enrichment profile"""
    general: GeneralInfo
    sequence: SequenceData
    biophysical: BiophysicalProperties
    charges: ComputedCharges
    predicted_features: PredictedFeatures
    timestamp: str

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'general': asdict(self.general),
            'sequence': asdict(self.sequence),
            'biophysical': asdict(self.biophysical),
            'charges': {
                'ph_values': self.charges.ph_values,
                'charges': self.charges.charges,
                'isoelectric_point': self.charges.isoelectric_point,
                'charge_at_ph_7_4': self.charges.charge_at_ph_7_4,
            },
            'predicted_features': {
                'signal_peptide': asdict(self.predicted_features.signal_peptide),
                'transmembrane': {
                    'is_membrane_protein': self.predicted_features.transmembrane.is_membrane_protein,
                    'tm_domains': [asdict(d) for d in self.predicted_features.transmembrane.tm_domains],
                    'prediction_tool': self.predicted_features.transmembrane.prediction_tool,
                },
                'disulfide_bonds': [asdict(d) for d in self.predicted_features.disulfide_bonds],
                'ptm_predictions': [asdict(p) for p in self.predicted_features.ptm_predictions],
            },
            'timestamp': self.timestamp,
        }

    def to_json(self, filename: Optional[str] = None) -> str:
        """Export to JSON file or string"""
        json_str = json.dumps(self.to_dict(), indent=2)
        if filename:
            with open(filename, 'w') as f:
                f.write(json_str)
            return f"Saved to {filename}"
        return json_str


# ============================================================================
# SECTION 2: Biophysical Calculations
# ============================================================================

class BiophysicalCalculator:
    """Compute sequence-based biophysical properties using BioPython"""

    @staticmethod
    def calculate_from_sequence(sequence: str) -> BiophysicalProperties:
        """
        Calculate all biophysical properties from amino acid sequence

        Args:
            sequence: Amino acid sequence (uppercase)

        Returns:
            BiophysicalProperties object with all computed values
        """
        # Ensure uppercase
        sequence = sequence.upper()

        # Create ProteinAnalysis object
        pa = ProteinAnalysis(sequence)

        # Calculate basic properties
        mw_da = pa.molecular_weight()
        pI = pa.isoelectric_point()
        aromaticity = pa.aromaticity()
        instability = pa.instability_index()
        gravy = pa.gravy()

        # Calculate aliphatic index (manually using definition)
        # Aliphatic Index = X(Ala) + a(Val) + b(Ile) + c(Leu)
        # where a, b, c are relative volumes of side chains
        aliphatic = BiophysicalCalculator._calculate_aliphatic_index(sequence)

        # Calculate secondary structure propensities (Chou-Fasman simplified)
        helix_frac, sheet_frac, coil_frac = BiophysicalCalculator._calculate_ss_propensities(sequence)

        # Count cysteines (from sequence only - no heuristics)
        cys_count = sequence.count('C')

        return BiophysicalProperties(
            molecular_weight_da=mw_da,
            molecular_weight_kda=mw_da / 1000,
            isoelectric_point_pi=pI,
            aromaticity=aromaticity,
            instability_index=instability,
            gravy=gravy,
            aliphatic_index=aliphatic,
            helix_fraction=helix_frac,
            sheet_fraction=sheet_frac,
            coil_fraction=coil_frac,
            cysteine_count=cys_count,
        )

    @staticmethod
    def _calculate_aliphatic_index(sequence: str) -> float:
        """
        Calculate aliphatic index of a protein.

        Aliphatic Index = X(Ala) + a*X(Val) + b*X(Ile) + c*X(Leu)
        where X is mole percent, a=2.9, b=3.9, c=3.9

        Reference:
        Ikai, A.J. (1980). Thermostability and aliphatic index of globular proteins.
        J. Biochem. 88, 1895-1898.
        """
        ala_pct = (sequence.count('A') / len(sequence)) * 100
        val_pct = (sequence.count('V') / len(sequence)) * 100
        ile_pct = (sequence.count('I') / len(sequence)) * 100
        leu_pct = (sequence.count('L') / len(sequence)) * 100

        aliphatic_index = ala_pct + (2.9 * val_pct) + (3.9 * ile_pct) + (3.9 * leu_pct)
        return aliphatic_index

    @staticmethod
    def _calculate_ss_propensities(sequence: str) -> tuple:
        """
        Simplified secondary structure propensity calculation (Chou-Fasman)

        Returns:
            (helix_fraction, sheet_fraction, coil_fraction)
        """
        # Chou-Fasman propensity values
        helix_props = {
            'A': 1.42, 'R': 0.98, 'N': 0.67, 'D': 1.01, 'C': 0.70,
            'Q': 1.11, 'E': 1.51, 'G': 0.57, 'H': 1.00, 'I': 1.08,
            'L': 1.21, 'K': 1.16, 'M': 1.45, 'F': 1.13, 'P': 0.57,
            'S': 0.77, 'T': 0.83, 'W': 1.08, 'Y': 0.69, 'V': 1.06
        }

        sheet_props = {
            'A': 0.83, 'R': 0.93, 'N': 0.89, 'D': 0.54, 'C': 1.19,
            'Q': 1.10, 'E': 0.37, 'G': 0.75, 'H': 0.87, 'I': 1.60,
            'L': 1.30, 'K': 0.74, 'M': 1.05, 'F': 1.38, 'P': 0.55,
            'S': 0.75, 'T': 1.19, 'W': 1.37, 'Y': 1.47, 'V': 1.70
        }

        coil_props = {
            'A': 0.66, 'R': 0.95, 'N': 1.56, 'D': 1.46, 'C': 1.19,
            'Q': 0.98, 'E': 0.57, 'G': 1.57, 'H': 0.91, 'I': 0.47,
            'L': 0.57, 'K': 1.01, 'M': 0.60, 'F': 0.60, 'P': 1.52,
            'S': 1.05, 'T': 0.83, 'W': 1.08, 'Y': 0.69, 'V': 0.50
        }

        # Calculate averages
        helix_avg = sum(helix_props.get(aa, 0.8) for aa in sequence) / len(sequence)
        sheet_avg = sum(sheet_props.get(aa, 0.8) for aa in sequence) / len(sequence)
        coil_avg = sum(coil_props.get(aa, 0.8) for aa in sequence) / len(sequence)

        # Normalize to percentages
        total = helix_avg + sheet_avg + coil_avg
        helix_frac = (helix_avg / total) * 100
        sheet_frac = (sheet_avg / total) * 100
        coil_frac = (coil_avg / total) * 100

        return helix_frac, sheet_frac, coil_frac

    @staticmethod
    def calculate_charge_profile(sequence: str, ph_range: tuple = (2, 12)) -> ComputedCharges:
        """
        Calculate protein charge across pH range using BioPython

        Args:
            sequence: Amino acid sequence
            ph_range: (min_pH, max_pH) tuple

        Returns:
            ComputedCharges object
        """
        # Use BioPython's built-in charge_at_pH method
        sequence = sequence.upper()
        pa = ProteinAnalysis(sequence)

        # Get isoelectric point
        pI = pa.isoelectric_point()

        # Calculate charges across pH range with 0.1 increments for smooth curve
        ph_start, ph_end = ph_range
        ph_values = [round(i * 0.1, 1) for i in range(int(ph_start * 10), int(ph_end * 10) + 1)]
        charges = [pa.charge_at_pH(pH) for pH in ph_values]

        # Get charge at physiological pH
        charge_7_4 = pa.charge_at_pH(7.4)

        return ComputedCharges(
            ph_values=ph_values,
            charges=charges,
            isoelectric_point=pI,
            charge_at_ph_7_4=charge_7_4,
        )


# ============================================================================
# SECTION 3: Structural Feature Predictions (Stubs)
# ============================================================================

class StructuralPredictor:
    """Predict structural/functional features (with remote API stubs)"""

    @staticmethod
    def predict_signal_peptide(sequence: str) -> SignalPeptidePrediction:
        """
        Predict signal peptide using SignalP 6.0 (stub - returns placeholder)

        Args:
            sequence: Amino acid sequence

        Returns:
            SignalPeptidePrediction object

        Note:
            For production, call actual SignalP API:
            POST to https://services.healthtech.dtu.dk/service.php
        """
        # STUB: For now, check if starts with M and simple heuristic
        is_sp = False
        cleavage = None
        sp_seq = None

        if len(sequence) > 15:
            # Very simple heuristic: check for hydrophobic region in first 30 aa
            first_30 = sequence[:30]
            hydrophobic_count = sum(1 for aa in first_30 if aa in 'AILMFVPW')
            if hydrophobic_count > 10:
                is_sp = True
                cleavage = 20  # Typical cleavage site position

        return SignalPeptidePrediction(
            present=is_sp,
            cleavage_site=cleavage,
            signal_peptide_seq=sequence[:cleavage] if cleavage else None,
            prediction_tool="SignalP 6.0 (stub - requires API)",
            confidence=None,
        )

    @staticmethod
    def predict_transmembrane_domains(sequence: str) -> TransmembranePrediction:
        """
        Predict transmembrane domains using TMHMM 2.0 (stub)

        Args:
            sequence: Amino acid sequence

        Returns:
            TransmembranePrediction object

        Note:
            For production, call actual TMHMM API or use local tool
        """
        # STUB: Simple hydrophobicity check
        is_tm = False
        tm_domains = []

        # Count hydrophobic amino acids
        hydrophobic = sum(1 for aa in sequence if aa in 'AILMFVPW')
        hydrophobic_frac = hydrophobic / len(sequence)

        # If > 40% hydrophobic, likely has TM domains
        if hydrophobic_frac > 0.4:
            is_tm = True
            # Rough approximation: 2 TM domains
            tm_domains = [
                TransmembraneDomain(1, 20, 40, "inside-out"),
                TransmembraneDomain(2, 100, 120, "outside-in"),
            ]

        return TransmembranePrediction(
            is_membrane_protein=is_tm,
            tm_domains=tm_domains,
            prediction_tool="TMHMM 2.0 (stub - requires API)",
        )

    @staticmethod
    def predict_disulfide_bonds(sequence: str) -> List[DisulfideBondPrediction]:
        """
        Predict disulfide bonds (stub - requires DiANNA or ML model)

        Args:
            sequence: Amino acid sequence

        Returns:
            List of DisulfideBondPrediction objects

        Note:
            For production, use:
            - DiANNA API
            - Local ML model (Cys-Rec, DiANNA server)
            - Or structural prediction (AlphaFold → contact map)
        """
        predictions = []

        # Find all cysteines
        cys_positions = [i + 1 for i, aa in enumerate(sequence) if aa == 'C']

        # STUB: Simple heuristic - pair nearby cysteines
        for i in range(0, len(cys_positions) - 1, 2):
            if i + 1 < len(cys_positions):
                predictions.append(DisulfideBondPrediction(
                    cys_position_1=cys_positions[i],
                    cys_position_2=cys_positions[i + 1],
                    confidence=0.5,  # Placeholder confidence
                    is_experimentally_confirmed=False,
                ))

        return predictions

    @staticmethod
    def predict_ptms(sequence: str) -> List[PTMPrediction]:
        """
        Predict post-translational modifications (stub)

        Args:
            sequence: Amino acid sequence

        Returns:
            List of PTMPrediction objects

        Note:
            For production, use:
            - NetPhos for phosphorylation
            - GPS (Group-based Phosphorylation Scoring) for kinase specificity
            - NetNGlyc for N-glycosylation
            - NetOGlyc for O-glycosylation
        """
        predictions = []

        # STUB: Simple PTM predictions based on sequence patterns
        # In production, these would come from dedicated ML tools

        # Phosphorylation sites (S/T followed by acidic)
        for i, aa in enumerate(sequence):
            if aa in 'ST' and i + 1 < len(sequence):
                next_aa = sequence[i + 1]
                if next_aa in 'DE':
                    predictions.append(PTMPrediction(
                        ptm_type='phosphorylation',
                        position=i + 1,
                        residue=aa,
                        confidence=0.6,
                        prediction_tool='Simple pattern (stub - use NetPhos)'
                    ))

        # N-glycosylation sites (N-x-S/T)
        for i in range(len(sequence) - 2):
            if sequence[i] == 'N' and sequence[i + 2] in 'ST':
                if sequence[i + 1] != 'P':  # Not Pro
                    predictions.append(PTMPrediction(
                        ptm_type='n_glycosylation',
                        position=i + 1,
                        residue='N',
                        confidence=0.7,
                        prediction_tool='Motif N-x-S/T (stub - use NetNGlyc)'
                    ))

        return predictions[:10]  # Limit to 10 predictions


# ============================================================================
# SECTION 4: Main Enrichment Service
# ============================================================================

class ProteinEnrichmentService:
    """Main service to orchestrate protein enrichment"""

    def __init__(self):
        """Initialize the service"""
        self.uniprot_client = UniProtClient()

    def _extract_or_predict_signal_peptide(self, sequence: str, features: List[Dict]) -> SignalPeptidePrediction:
        """
        Extract signal peptide from UniProt, or predict if not found

        Args:
            sequence: Amino acid sequence
            features: List of UniProt features

        Returns:
            SignalPeptidePrediction object
        """
        # Check if UniProt has signal peptide annotation
        for feature in features:
            if feature.get('feature_type') in ['Signal', 'Signal peptide']:
                start = feature.get('start_position')
                end = feature.get('end_position')
                if start and end:
                    signal_seq = sequence[start-1:end] if end <= len(sequence) else sequence[start-1:]
                    return SignalPeptidePrediction(
                        present=True,
                        cleavage_site=end,
                        signal_peptide_seq=signal_seq,
                        prediction_tool=f"UniProt (experimentally confirmed)",
                        confidence=1.0,  # Experimentally confirmed
                    )

        # Fall back to prediction
        return StructuralPredictor.predict_signal_peptide(sequence)

    def _extract_or_predict_transmembrane(self, sequence: str, features: List[Dict]) -> TransmembranePrediction:
        """
        Extract transmembrane domains from UniProt, or predict if not found

        Args:
            sequence: Amino acid sequence
            features: List of UniProt features

        Returns:
            TransmembranePrediction object
        """
        # Check if UniProt has TM domain annotations
        tm_domains = []
        for feature in features:
            if feature.get('feature_type') == 'TRANSMEM':
                start = feature.get('start_position')
                end = feature.get('end_position')
                if start and end:
                    tm_domains.append(TransmembraneDomain(
                        number=len(tm_domains) + 1,
                        start=start,
                        end=end,
                        orientation=feature.get('description', 'unknown'),
                    ))

        if tm_domains:
            return TransmembranePrediction(
                is_membrane_protein=True,
                tm_domains=tm_domains,
                prediction_tool=f"UniProt (experimentally confirmed)",
            )

        # Fall back to prediction
        return StructuralPredictor.predict_transmembrane_domains(sequence)

    def _extract_or_predict_disulfides(self, features: List[Dict]) -> List[DisulfideBondPrediction]:
        """
        Extract disulfide bonds from UniProt only - NO heuristic predictions

        Args:
            features: List of UniProt features

        Returns:
            List of DisulfideBondPrediction objects (only experimentally confirmed from UniProt)
        """
        # Extract ONLY confirmed disulfide bonds from UniProt
        disulfides = []
        for feature in features:
            if feature.get('feature_type') in ['Disulfide bond', 'DISULFID']:
                start = feature.get('start_position')
                end = feature.get('end_position')
                if start and end:
                    disulfides.append(DisulfideBondPrediction(
                        cys_position_1=start,
                        cys_position_2=end,
                        confidence=1.0,  # Experimentally confirmed
                        is_experimentally_confirmed=True,
                    ))

        # Return only confirmed data - NO fallback to heuristics or predictions
        # Not all cysteines form disulfide bonds. Only report what UniProt confirms.
        return disulfides

    def _extract_or_predict_ptms(self, sequence: str, features: List[Dict]) -> List[PTMPrediction]:
        """
        Extract ALL PTMs/Processing from UniProt PTM/Processing section, or predict if not found

        Captures all UniProt PTM/Processing subsections:

        Processing events:
        - Initiator methionine: Cleaved initiator methionine
        - Signal: Sequence targeting to secretory pathway or periplasmic space
        - Propeptide: Cleaved during maturation or activation
        - Transit peptide: For organelle targeting
        - Chain: Mature polypeptide chain
        - Peptide: Active peptide

        Post-translational modifications:
        - Modified residue: Phosphorylation, methylation, acetylation, etc.
        - Glycosylation: N-linked and O-linked glycans
        - Lipidation: Covalently attached lipids
        - Disulfide bond: Cysteine cross-links (also handled separately)
        - Cross-link: Residue-residue cross-links

        Args:
            sequence: Amino acid sequence
            features: List of UniProt features

        Returns:
            List of PTMPrediction objects
        """
        # All UniProt PTM/Processing feature types
        ptm_feature_types = {
            # Processing events
            'Initiator methionine': 'initiator_methionine',
            'Signal': 'signal_peptide',
            'Propeptide': 'propeptide',
            'Transit peptide': 'transit_peptide',
            'Chain': 'chain',
            'Peptide': 'active_peptide',

            # Post-translational modifications
            'Modified residue': 'modification',
            'Glycosylation': 'glycosylation',
            'Lipidation': 'lipidation',
            'Disulfide bond': 'disulfide_bond',
            'Cross-link': 'cross_link',
        }

        # Extract all known PTM/Processing features from UniProt
        ptms = []
        for feature in features:
            feature_type = feature.get('feature_type')

            # Check if this is a PTM/Processing feature
            if feature_type in ptm_feature_types:
                start_pos = feature.get('start_position')
                end_pos = feature.get('end_position')
                description = feature.get('description', '')

                # Determine PTM type
                ptm_type = ptm_feature_types.get(feature_type, 'modification')

                # For Modified residue, parse more specific type from description
                if feature_type == 'Modified residue':
                    desc_lower = description.lower()
                    if 'phosph' in desc_lower:
                        ptm_type = 'phosphorylation'
                    elif 'glyco' in desc_lower:
                        ptm_type = 'glycosylation'
                    elif 'methyl' in desc_lower:
                        ptm_type = 'methylation'
                    elif 'acetyl' in desc_lower:
                        ptm_type = 'acetylation'
                    elif 'sumo' in desc_lower:
                        ptm_type = 'sumoylation'
                    elif 'ubiquit' in desc_lower:
                        ptm_type = 'ubiquitination'
                    elif 'nitro' in desc_lower:
                        ptm_type = 'nitrosylation'
                    elif 'palmitoyl' in desc_lower:
                        ptm_type = 'palmitoylation'
                    elif 'myristoyl' in desc_lower:
                        ptm_type = 'myristoylation'

                # Use start position as primary identifier
                position = start_pos
                if position and position <= len(sequence):
                    residue = sequence[position - 1]
                    ptms.append(PTMPrediction(
                        ptm_type=ptm_type,
                        position=position,
                        residue=residue,
                        confidence=1.0,  # Experimentally confirmed from UniProt
                        prediction_tool=f"UniProt (experimentally confirmed)",
                    ))

        if ptms:
            return ptms

        # Fall back to prediction if no UniProt PTM/Processing data found
        return StructuralPredictor.predict_ptms(sequence)

    def enrich_protein(self, accession_id: str) -> ProteinEnrichmentProfile:
        """
        Enrich a single protein with full data and predictions

        Args:
            accession_id: UniProt accession ID

        Returns:
            ProteinEnrichmentProfile with all data
        """
        print(f"\n🔬 Enriching {accession_id}...")

        # Step 1: Fetch from UniProt
        print("  → Fetching from UniProt...")
        protein_json = self.uniprot_client.fetch_protein(accession_id)
        if not protein_json:
            raise ValueError(f"Failed to fetch {accession_id}")

        processor = UniProtDataProcessor(protein_json)
        protein_data = processor.get_main_protein_data()

        # Extract basic info
        sequence = protein_data.get('sequence', '').upper()
        if not sequence:
            raise ValueError(f"No sequence for {accession_id}")

        general = GeneralInfo(
            protein_name=protein_data.get('protein_recommended_name', 'Unknown'),
            gene_name=protein_data.get('gene_name'),
            organism=protein_data.get('organism_scientific_name', 'Unknown'),
            organism_taxon_id=protein_data.get('organism_taxon_id'),
            uniprot_id=accession_id,
            entry_type=protein_data.get('entry_type', ''),
            sequence_length=len(sequence),
            reviewed='reviewed' in protein_data.get('entry_type', '').lower(),
            created_date=protein_data.get('created_date'),
            updated_date=protein_data.get('updated_date'),
        )

        sequence_data = SequenceData(
            fasta_sequence=sequence,
            length=len(sequence),
        )

        # Step 2: Calculate biophysical properties
        print("  → Computing biophysical properties...")
        biophysical = BiophysicalCalculator.calculate_from_sequence(sequence)

        # Step 3: Calculate charge profile
        print("  → Computing charge profile...")
        charges = BiophysicalCalculator.calculate_charge_profile(sequence)

        # Step 4: Get known features from UniProt, then predict missing ones
        print("  → Extracting known features from UniProt...")
        # Get features from the raw UniProt data
        features = []
        for feature in protein_json.get('features', []):
            feature_entry = {
                'feature_type': feature.get('type', ''),
                'description': feature.get('description', ''),
                'start_position': None,
                'end_position': None,
            }
            if 'location' in feature:
                location = feature['location']
                if 'start' in location:
                    feature_entry['start_position'] = location['start'].get('value')
                if 'end' in location:
                    feature_entry['end_position'] = location['end'].get('value')
            features.append(feature_entry)

        # Try to extract from UniProt first, fall back to predictions
        signal_peptide = self._extract_or_predict_signal_peptide(sequence, features)
        transmembrane = self._extract_or_predict_transmembrane(sequence, features)
        disulfides = self._extract_or_predict_disulfides(features)
        ptms = self._extract_or_predict_ptms(sequence, features)

        predicted_features = PredictedFeatures(
            signal_peptide=signal_peptide,
            transmembrane=transmembrane,
            disulfide_bonds=disulfides,
            ptm_predictions=ptms,
        )

        # Step 5: Create profile
        print("  → Creating enrichment profile...")
        profile = ProteinEnrichmentProfile(
            general=general,
            sequence=sequence_data,
            biophysical=biophysical,
            charges=charges,
            predicted_features=predicted_features,
            timestamp=datetime.now().isoformat(),
        )

        print(f"  ✓ {accession_id} enriched successfully")
        return profile

    def enrich_protein_batch(self, accession_ids: List[str]) -> List[ProteinEnrichmentProfile]:
        """
        Enrich multiple proteins

        Args:
            accession_ids: List of UniProt accession IDs

        Returns:
            List of ProteinEnrichmentProfile objects
        """
        results = []
        for acc_id in accession_ids:
            try:
                profile = self.enrich_protein(acc_id)
                results.append(profile)
            except Exception as e:
                print(f"  ✗ Failed to enrich {acc_id}: {e}")
        return results


# ============================================================================
# SECTION 5: Public API
# ============================================================================

def enrich_protein(accession_id: str) -> ProteinEnrichmentProfile:
    """
    Quick function to enrich a single protein

    Args:
        accession_id: UniProt accession ID

    Returns:
        ProteinEnrichmentProfile

    Example:
        profile = enrich_protein('P04637')
        profile.to_json('p53_enriched.json')
    """
    service = ProteinEnrichmentService()
    return service.enrich_protein(accession_id)


def enrich_protein_batch(accession_ids: List[str]) -> List[ProteinEnrichmentProfile]:
    """
    Enrich multiple proteins

    Args:
        accession_ids: List of UniProt accession IDs

    Returns:
        List of ProteinEnrichmentProfile objects

    Example:
        profiles = enrich_protein_batch(['P04637', 'P68871', 'P15018'])
        for profile in profiles:
            profile.to_json(f"{profile.general.uniprot_id}_enriched.json")
    """
    service = ProteinEnrichmentService()
    return service.enrich_protein_batch(accession_ids)


# ============================================================================
# SECTION 6: CLI
# ============================================================================

def main():
    """Command-line interface"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Unified Protein Enrichment Service',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Enrich single protein
  python protein_enrichment_service.py P04637

  # Enrich and export to JSON
  python protein_enrichment_service.py P04637 -o p53_enriched.json

  # Batch enrichment
  python protein_enrichment_service.py P04637 P68871 P15018 -o enriched/
        """
    )

    parser.add_argument(
        'accession_ids',
        nargs='+',
        help='One or more UniProt accession IDs'
    )

    parser.add_argument(
        '-o', '--output',
        help='Output JSON filename or directory (for batch)'
    )

    parser.add_argument(
        '--pretty',
        action='store_true',
        help='Pretty-print JSON to console'
    )

    args = parser.parse_args()

    print("\n" + "=" * 80)
    print("UNIFIED PROTEIN ENRICHMENT SERVICE")
    print("=" * 80)

    # Process proteins
    profiles = enrich_protein_batch(args.accession_ids)

    # Output
    if args.output:
        if len(profiles) == 1:
            # Single protein - save as single file
            filename = args.output if args.output.endswith('.json') else f"{args.output}.json"
            profiles[0].to_json(filename)
            print(f"\n✓ Saved to: {filename}")
        else:
            # Batch - create directory
            import os
            os.makedirs(args.output, exist_ok=True)
            for profile in profiles:
                filename = f"{args.output}/{profile.general.uniprot_id}_enriched.json"
                profile.to_json(filename)
                print(f"✓ Saved: {filename}")

    # Console output
    if args.pretty:
        for profile in profiles:
            print("\n" + "=" * 80)
            print(f"ENRICHMENT PROFILE: {profile.general.uniprot_id}")
            print("=" * 80)
            print(json.dumps(profile.to_dict(), indent=2))


if __name__ == '__main__':
    main()
