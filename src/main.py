#!/usr/bin/env python3
"""
Main entry point for protein analysis application.
Simple, modular interface for generating protein reports.
"""

import sys
import json
import os
from pathlib import Path
from typing import Dict, Any

from src.api_client import APIClient
from src.protein_analyzer import ProteinAnalyzer
from src.services.protein_enrichment_service import ProteinEnrichmentService
from src.data_manager import DataManager
from src.config import REPORT, ANALYSIS


class ProteinReportGenerator:
    """Main application for generating protein analysis reports."""

    def __init__(self):
        """Initialize report generator."""
        self.api_client = APIClient()
        self.enrichment_service = ProteinEnrichmentService()
        self.data_manager = DataManager()

    def generate_report(self, accession_id: str) -> bool:
        """
        Generate protein analysis report.

        Args:
            accession_id: UniProt accession ID

        Returns:
            True if successful, False otherwise
        """
        print(f"\n{'='*60}")
        print(f"  PROTEIN ANALYSIS REPORT GENERATOR")
        print(f"{'='*60}")
        print(f"\n🔬 Processing: {accession_id}")

        # 1. Fetch all protein data from APIs
        protein_data = self.api_client.fetch_all_protein_data(accession_id)

        if not protein_data['sequence']:
            print(f"❌ Error: Could not fetch sequence for {accession_id}")
            return False

        sequence = protein_data['sequence']
        features = protein_data['features']
        pdb_id = protein_data['pdb_id']

        # 2. Extract PTMs
        print(f"   🔧 Analyzing post-translational modifications...")
        ptm_summary = self._extract_ptm_summary(features) if features else {'total_ptms': 0, 'by_type': {}, 'positions': []}
        if ptm_summary:
            print(f"   ✅ Found {ptm_summary['total_ptms']} PTM(s)")

        # 3. Analyze biophysical properties
        print(f"   🧮 Calculating biophysical properties...")
        analyzer = ProteinAnalyzer(sequence)
        raw_properties = analyzer.get_all_properties()

        # Transform property keys to match template expectations
        # Note: ProtParam from BioPython provides helix, sheet, turn percentages
        # Coil is not explicitly calculated, so we don't include it
        helix = raw_properties['helix_predicted']
        turn = raw_properties['turn_predicted']
        sheet = raw_properties['sheet_predicted']

        # Calculate experimental secondary structure from UniProt annotations
        features_list = features if isinstance(features, list) else (features.get('features', []) if features else [])
        experimental_ss = self._calculate_experimental_secondary_structure(features_list, len(sequence))

        properties = {
            'molecular_weight': raw_properties['molecular_weight_kda'],
            'isoelectric_point': raw_properties['isoelectric_point'],
            'gravy_score': raw_properties['gravy_score'],
            'aromaticity': raw_properties['aromaticity_percent'],
            'instability_index': raw_properties['instability_index'],
            'cysteine_count': raw_properties['cysteine_count'],
            'extinction_coeff_reduced': raw_properties['extinction_coeff_reduced'],
            'extinction_coeff_oxidized': raw_properties['extinction_coeff_oxidized'],
            'helix_predicted': helix,
            'turn_predicted': turn,
            'sheet_predicted': sheet,
            'helix_experimental': experimental_ss['helix_experimental'],
            'turn_experimental': experimental_ss['turn_experimental'],
            'sheet_experimental': experimental_ss['sheet_experimental'],
        }
        print(f"   ✅ Calculated properties")
        if experimental_ss['helix_experimental'] > 0 or experimental_ss['sheet_experimental'] > 0:
            print(f"   ✅ Found experimental secondary structure annotations")

        # 4. Generate charge profile
        print(f"   🔋 Generating charge profile...")
        try:
            enriched_profile = self.enrichment_service.enrich_protein(accession_id)
            charge_profile_json = self._generate_charge_profile_json(enriched_profile)
            print(f"   ✅ Generated charge profile")
        except Exception as e:
            print(f"   ⚠️  Could not generate charge profile: {e}")
            charge_profile_json = json.dumps({
                "ph_values": [],
                "charges": [],
                "isoelectric_point": 0,
                "charge_at_ph_7_4": 0
            })

        # 5. Prepare data for database
        protein_name = self._extract_protein_name(protein_data['entry']) if protein_data['entry'] else "Unknown"
        gene_name = self._extract_gene_name(protein_data['entry']) if protein_data['entry'] else "N/A"
        organism = self._extract_organism(protein_data['entry']) if protein_data['entry'] else "Unknown"
        subcellular_location = self._extract_subcellular_location(protein_data['entry']) if protein_data['entry'] else "Unknown"

        charge_profile_dict = json.loads(charge_profile_json)

        # Handle features list (from UniProt API)
        features_list = features if isinstance(features, list) else (features.get('features', []) if features else [])

        database_data = {
            'protein_name': protein_name,
            'gene_name': gene_name,
            'organism': organism,
            'subcellular_location': subcellular_location,
            'sequence': sequence,
            'sequence_length': len(sequence),
            **properties,
            'pdb_id': pdb_id,
            'features': self._enhance_features(features_list),
            'ptm_summary': ptm_summary,
            'secondary_structure': {
                'helix_predicted': properties['helix_predicted'],
                'turn_predicted': properties['turn_predicted'],
                'sheet_predicted': properties['sheet_predicted'],
                'helix_experimental': properties['helix_experimental'],
                'turn_experimental': properties['turn_experimental'],
                'sheet_experimental': properties['sheet_experimental'],
            },
            'charge_profile': charge_profile_dict,
            'disulfide_bonds': [],
        }

        # 6. Save to database
        print(f"   💾 Saving data to database...")
        protein_id = self.data_manager.save_protein(accession_id, database_data)
        if not protein_id:
            return False

        # 7. Generate report
        print(f"   📝 Rendering template...")
        self._generate_html_report(accession_id, database_data)

        # 8. Save report metadata
        output_filename = f"{REPORT['prefix']}{accession_id}.html"
        output_path = Path(REPORT['output_dir']) / output_filename
        if output_path.exists():
            file_size_kb = output_path.stat().st_size / 1024
            self.data_manager.save_report(protein_id, output_filename, str(output_path), file_size_kb)

        print(f"   ✅ Report saved: {output_filename}")
        print(f"\n{'='*60}")

        return True

    def _extract_ptm_summary(self, features) -> Dict[str, Any]:
        """Extract PTM summary from features."""
        ptm_summary = {
            'total_ptms': 0,
            'by_type': {},
            'details': [],
            'peptides': []
        }

        # Handle features as a list (from UniProt API)
        features_list = features if isinstance(features, list) else (features.get('features', []) if isinstance(features, dict) else [])

        # PTM types to look for in UniProt
        ptm_types = {
            'Modified residue', 'Lipidation', 'Glycosylation', 'Phosphorylation',
            'N-linked glycosylation', 'O-linked glycosylation', 'Sulfation',
            'Acetylation', 'Amidation', 'GPI-anchor', 'Ubiquitination',
            'Sumoylation', 'Nitrosylation', 'S-nitrosylation', 'Neddylation',
            'Disulfide bond', 'Proteolytic cleavage', 'Cross-link', 'Hydroxylation',
            'ADP-ribosylation', 'Farnesylation', 'Palmitoylation', 'Methylation'
        }

        # Peptide types
        peptide_types = {
            'Signal peptide', 'Propeptide', 'Transit peptide', 'Peptide'
        }

        for feature in features_list:
            feature_type = feature.get('type', 'Unknown')

            # Check if this is a PTM feature
            if feature_type in ptm_types:
                location = feature.get('location', {})
                start_pos = 0
                end_pos = 0

                if isinstance(location, dict):
                    if 'start' in location:
                        start_pos = location['start'].get('value', 0)
                    if 'end' in location:
                        end_pos = location['end'].get('value', 0)

                description = feature.get('description', feature_type)

                ptm_summary['total_ptms'] += 1
                ptm_summary['by_type'][feature_type] = ptm_summary['by_type'].get(feature_type, 0) + 1
                ptm_summary['details'].append({
                    'position': start_pos,
                    'end_position': end_pos,
                    'type': feature_type,
                    'description': description
                })

            # Check if this is a peptide feature
            elif feature_type in peptide_types:
                location = feature.get('location', {})
                start_pos = 0
                end_pos = 0

                if isinstance(location, dict):
                    if 'start' in location:
                        start_pos = location['start'].get('value', 0)
                    if 'end' in location:
                        end_pos = location['end'].get('value', 0)

                description = feature.get('description', feature_type)

                ptm_summary['peptides'].append({
                    'start_position': start_pos,
                    'end_position': end_pos,
                    'type': feature_type,
                    'description': description
                })

        return ptm_summary

    def _generate_charge_profile_json(self, profile) -> str:
        """Generate charge profile JSON."""
        data = {
            "ph_values": profile.charges.ph_values,
            "charges": profile.charges.charges,
            "isoelectric_point": profile.charges.isoelectric_point,
            "charge_at_ph_7_4": profile.charges.charge_at_ph_7_4
        }
        return json.dumps(data)

    def _extract_protein_name(self, entry: Dict) -> str:
        """Extract protein name from UniProt entry."""
        try:
            prot_desc = entry.get('proteinDescription', {})
            rec_name = prot_desc.get('recommendedName', {})
            return rec_name.get('fullName', {}).get('value', 'Unknown')
        except:
            return 'Unknown'

    def _extract_gene_name(self, entry: Dict) -> str:
        """Extract gene name from UniProt entry."""
        try:
            genes = entry.get('genes', [])
            if genes and 'geneName' in genes[0]:
                return genes[0]['geneName'].get('value', 'N/A')
        except:
            pass
        return 'N/A'

    def _extract_organism(self, entry: Dict) -> str:
        """Extract organism name from UniProt entry."""
        try:
            organism = entry.get('organism', {})
            return organism.get('scientificName', 'Unknown')
        except:
            return 'Unknown'

    def _extract_subcellular_location(self, entry: Dict) -> str:
        """Extract subcellular location from UniProt entry."""
        try:
            comments = entry.get('comments', [])
            locations = []

            for comment in comments:
                if comment.get('commentType') == 'SUBCELLULAR LOCATION':
                    # Extract location values
                    subcell_locs = comment.get('subcellularLocations', [])
                    for loc in subcell_locs:
                        location_info = loc.get('location', {})
                        location_value = location_info.get('value', '')
                        if location_value:
                            locations.append(location_value)

            return '; '.join(locations) if locations else 'Unknown'
        except:
            return 'Unknown'

    def _calculate_experimental_secondary_structure(self, features: list, sequence_length: int) -> Dict[str, float]:
        """
        Calculate experimental secondary structure percentages from UniProt annotations.

        UniProt provides experimentally-determined or literature-based secondary structure
        assignments using features like 'Helix', 'Beta strand', and 'Turn'.

        Args:
            features: List of features from UniProt
            sequence_length: Length of protein sequence

        Returns:
            Dictionary with percentages for helix, sheet, turn, coil
        """
        # Track residues assigned to each structure type
        helix_residues = set()
        sheet_residues = set()
        turn_residues = set()

        for feature in features:
            feature_type = feature.get('type', '')
            location = feature.get('location', {})

            if not isinstance(location, dict):
                continue

            start = location.get('start', {})
            end = location.get('end', {})

            if not isinstance(start, dict) or not isinstance(end, dict):
                continue

            start_pos = start.get('value', 0)
            end_pos = end.get('value', 0)

            if start_pos == 0 or end_pos == 0:
                continue

            # Map feature types to secondary structure
            if feature_type == 'Helix':
                helix_residues.update(range(start_pos, end_pos + 1))
            elif feature_type == 'Beta strand':
                sheet_residues.update(range(start_pos, end_pos + 1))
            elif feature_type == 'Turn':
                turn_residues.update(range(start_pos, end_pos + 1))

        # Calculate percentages (only for explicitly annotated structures, not coil)
        total = sequence_length
        helix_pct = (len(helix_residues) / total * 100) if total > 0 else 0
        sheet_pct = (len(sheet_residues) / total * 100) if total > 0 else 0
        turn_pct = (len(turn_residues) / total * 100) if total > 0 else 0

        return {
            'helix_experimental': helix_pct,
            'sheet_experimental': sheet_pct,
            'turn_experimental': turn_pct,
        }

    def _get_observed_secondary_structure(self, pdb_id: str) -> Dict[str, float]:
        """
        Get observed secondary structure from PDB.

        Note: This would require PDB API calls and DSSP calculations.
        For now, return empty dict to indicate not available.
        """
        # TODO: Implement PDB structure analysis
        # This would involve:
        # 1. Fetching structure from PDB using pdb_id
        # 2. Running DSSP or similar algorithm
        # 3. Extracting secondary structure assignments
        # 4. Calculating percentages for each residue's assignment
        return {}

    def _enhance_features(self, features: list) -> list:
        """Add tooltips to features and normalize types for Nightingale viewer."""
        # Map UniProt feature types to standardized types for Nightingale viewer
        type_mapping = {
            'Chain': 'CHAIN',
            'Domain': 'DOMAIN',
            'Region': 'REGION',
            'Site': 'SITE',
            'Binding site': 'BINDING',
            'Disulfide bond': 'DISULFID',
            'Beta strand': 'STRAND',
            'Helix': 'HELIX',
            'Turn': 'TURN',
            'Modified residue': 'MOD_RES',
            'Glycosylation': 'CARBOHYD',
            'Initiator methionine': 'INIT_MET',
            'Natural variant': 'VAR_SEQ',
            'Peptide': 'PEPTIDE',
            'Sequence conflict': 'CONFLICT',
            'Mutagenesis': 'MUTAGEN',
            'Calcium binding': 'CA_BIND',
            'DNA binding': 'DNA_BIND',
            'Lipidation': 'LIPID',
            'Metal binding': 'METAL',
            'Zinc finger': 'ZN_FING',
            'Disulfide bond (establishment)': 'DISULFID',
            'Signal peptide': 'SIGNAL',
            'Propeptide': 'PROPEP',
            'Transmembrane': 'TRANSMEM',
            'Cross-link': 'CROSSLNK',
            'Alternative sequence': 'ALT_SEQ',
            'Motif': 'MOTIF',
            'Compositional bias': 'COMPBIAS',
            'Splice variant': 'SPLICEVAR',
            'Topological domain': 'TOPO_DOM',
            'Hydrophobic region': 'HYDROPHOBIC',
            'Coiled-coil': 'COIL_COIL',
            'Acetylation': 'ACETYL',
            'Phosphorylation': 'PHOSPHO',
            'Ubiquitination': 'UBIQUIT',
            'Sumoylation': 'SUMO',
            'S-nitrosylation': 'NITRO',
            'N-linked glycosylation': 'CARBOHYD',
            'O-linked glycosylation': 'CARBOHYD',
            'Hydroxylation': 'HYDROXY',
            'ADP-ribosylation': 'ADPRIB',
            'Farnesylation': 'FARNESL',
            'Palmitoylation': 'PALMI',
            'GPI-anchor': 'GPI',
            'Neddylation': 'NEDDYL',
            'PROTEOMICS_PTM': 'PROTEOMICS_PTM',
        }

        # Define which types are PTM (to be combined), peptide, or disulfide
        ptm_types = {
            'MOD_RES', 'CARBOHYD', 'ACETYL', 'PHOSPHO', 'UBIQUIT', 'SUMO',
            'NITRO', 'LIPID', 'CROSSLNK', 'HYDROXY', 'ADPRIB', 'FARNESL',
            'PALMI', 'GPI', 'NEDDYL', 'INIT_MET', 'PROTEOMICS_PTM'
        }
        peptide_types = {'PEPTIDE', 'SIGNAL', 'PROPEP'}

        enhanced = []
        for feature in features:
            # Handle both UniProt (location.start.value) and EBI (begin/end) formats
            location = feature.get('location', {})
            if isinstance(location, dict) and 'start' in location:
                # UniProt format
                start = int(location['start'].get('value', 0))
                end = int(location.get('end', {}).get('value', start))
            else:
                # EBI format
                start = int(feature.get('begin') or feature.get('start', 0))
                end = int(feature.get('end') or start)

            size = end - start + 1

            # Normalize the feature type for the Nightingale viewer
            original_type = feature.get('type', 'Unknown')
            standardized_type = type_mapping.get(original_type, original_type)

            # Assign category for proper grouping in viewer
            if standardized_type in ptm_types:
                category = 'PTM'
            elif standardized_type in peptide_types:
                category = 'PEPTIDE'
            elif standardized_type == 'DISULFID':
                category = 'DISULFID'
            else:
                category = 'STRUCTURE'

            # Build tooltip
            tooltip = f"{feature.get('description', original_type)}\n"
            tooltip += f"Position: {start}-{end} aa\nSize: {size} aa"

            feature_copy = feature.copy()
            feature_copy['start'] = start
            feature_copy['end'] = end
            feature_copy['title'] = tooltip
            feature_copy['type'] = standardized_type  # Override with standardized type
            feature_copy['category'] = category  # Add category for filtering
            enhanced.append(feature_copy)

        return enhanced

    def _generate_html_report(self, accession_id: str, data: Dict[str, Any]):
        """Generate HTML report from template."""
        try:
            # Read the HTML template
            template_path = Path(REPORT['template'])
            if not template_path.exists():
                print(f"   ⚠️  Template not found: {template_path}")
                return

            with open(template_path, 'r') as f:
                html_content = f.read()

            # Replace all placeholders with actual data
            replacements = {
                '{{ACCESSION_ID}}': accession_id,
                '{{PROTEIN_NAME}}': data.get('protein_name', 'Unknown'),
                '{{GENE_NAME}}': data.get('gene_name', 'N/A'),
                '{{ORGANISM}}': data.get('organism', 'Unknown'),
                '{{SEQUENCE_LENGTH}}': str(data.get('sequence_length', 0)),
                '{{MOLECULAR_WEIGHT}}': f"{data.get('molecular_weight', 0):.2f}",
                '{{ISOELECTRIC_POINT}}': f"{data.get('isoelectric_point', 0):.2f}",
                '{{GRAVY_SCORE}}': f"{data.get('gravy_score', 0):.4f}",
                '{{AROMATICITY}}': f"{data.get('aromaticity', 0):.4f}",
                '{{INSTABILITY_INDEX}}': f"{data.get('instability_index', 0):.2f}",
                '{{CYSTEINE_COUNT}}': str(data.get('cysteine_count', 0)),
                '{{PDB_ID}}': data.get('pdb_id', 'N/A'),
                '{{HELIX_FRACTION}}': f"{data.get('helix_predicted', 0):.2f}",
                '{{SHEET_FRACTION}}': f"{data.get('sheet_predicted', 0):.2f}",
                '{{COIL_FRACTION}}': f"{data.get('coil_predicted', 0):.2f}",
                '{{TOTAL_PTMS}}': str(data.get('ptm_summary', {}).get('total_ptms', 0)),
                '{{CHARGE_PROFILE_JSON}}': json.dumps(data.get('charge_profile', {})),
                '{{SEQUENCE}}': data.get('sequence', ''),
                '{{FEATURES_JSON}}': json.dumps(data.get('features', [])),
                '{{PROTPARAM_TABLE}}': self._generate_protparam_table(data),
                '{{PTM_SECTION}}': self._generate_ptm_section(data),
            }

            # Apply replacements
            for placeholder, value in replacements.items():
                html_content = html_content.replace(placeholder, str(value))

            # Write the output file
            output_filename = f"{REPORT['prefix']}{accession_id}.html"
            output_path = Path(REPORT['output_dir']) / output_filename

            # Create output directory if it doesn't exist
            output_path.parent.mkdir(parents=True, exist_ok=True)

            with open(output_path, 'w') as f:
                f.write(html_content)

            print(f"   ✅ Report generated: {output_path}")
        except Exception as e:
            print(f"   ⚠️  Error generating report: {e}")

    def _generate_protparam_table(self, data: Dict[str, Any]) -> str:
        """Generate ProtParam card-based layout HTML."""
        return f"""
        <div class="biophysical-container">
            <!-- Biophysical Properties Card -->
            <div class="property-card full-width">
                <div class="card-title">Biophysical Properties</div>
                <div class="property-grid-large">
                    <div class="property-item">
                        <div class="property-label">Molecular Weight</div>
                        <div class="property-value">{data.get('molecular_weight', 0):.2f} <span class="unit">kDa</span></div>
                    </div>
                    <div class="property-item">
                        <div class="property-label">Isoelectric Point</div>
                        <div class="property-value">{data.get('isoelectric_point', 0):.2f} <span class="unit">pI</span></div>
                    </div>
                    <div class="property-item">
                        <div class="property-label">GRAVY Score</div>
                        <div class="property-value">{data.get('gravy_score', 0):.4f}</div>
                        <div class="property-note">Hydropathy index</div>
                    </div>
                    <div class="property-item">
                        <div class="property-label">Instability Index</div>
                        <div class="property-value">{data.get('instability_index', 0):.2f}</div>
                        <div class="property-note">Stability prediction</div>
                    </div>
                    <div class="property-item">
                        <div class="property-label">Extinction Coeff. (Reduced)</div>
                        <div class="property-value">{data.get('extinction_coeff_reduced', 0):.0f} <span class="unit">M⁻¹cm⁻¹</span></div>
                    </div>
                    <div class="property-item">
                        <div class="property-label">Extinction Coeff. (Oxidized)</div>
                        <div class="property-value">{data.get('extinction_coeff_oxidized', 0):.0f} <span class="unit">M⁻¹cm⁻¹</span></div>
                    </div>
                    <div class="property-item">
                        <div class="property-label">Cysteine Count</div>
                        <div class="property-value">{data.get('cysteine_count', 0)}</div>
                    </div>
                    <div class="property-item">
                        <div class="property-label">Aromaticity</div>
                        <div class="property-value">{data.get('aromaticity', 0):.2f}%</div>
                    </div>
                </div>
            </div>

            <!-- Subcellular Location Card -->
            <div class="property-card full-width">
                <div class="card-title">Subcellular Location</div>
                <div class="location-text">{data.get('subcellular_location', 'Unknown')}</div>
            </div>

            <!-- Secondary Structure Card -->
            <div class="property-card full-width">
                <div class="card-title">Secondary Structure</div>
                <div class="structure-percentages">
                    <div class="structure-row">
                        <div class="structure-row-label">Theoretical</div>
                        <div class="structure-values">
                            <div class="structure-value-item">
                                <span class="structure-label">α-Helix</span>
                                <span class="structure-percent">{data.get('helix_predicted', 0):.1f}%</span>
                            </div>
                            <div class="structure-value-item">
                                <span class="structure-label">β-Sheet</span>
                                <span class="structure-percent">{data.get('sheet_predicted', 0):.1f}%</span>
                            </div>
                            <div class="structure-value-item">
                                <span class="structure-label">β-Turn</span>
                                <span class="structure-percent">{data.get('turn_predicted', 0):.1f}%</span>
                            </div>
                        </div>
                    </div>
                    <div class="structure-row">
                        <div class="structure-row-label">Experimental</div>
                        <div class="structure-values">
                            <div class="structure-value-item">
                                <span class="structure-label">α-Helix</span>
                                <span class="structure-percent">{data.get('helix_experimental', 0):.1f}%</span>
                            </div>
                            <div class="structure-value-item">
                                <span class="structure-label">β-Sheet</span>
                                <span class="structure-percent">{data.get('sheet_experimental', 0):.1f}%</span>
                            </div>
                            <div class="structure-value-item">
                                <span class="structure-label">β-Turn</span>
                                <span class="structure-percent">{data.get('turn_experimental', 0):.1f}%</span>
                            </div>
                        </div>
                    </div>
                </div>
                <div class="structure-note-text">Theoretical: DSSP/STRIDE algorithm | Experimental: UniProt annotations</div>
            </div>
        </div>
        """

    def _generate_ptm_section(self, data: Dict[str, Any]) -> str:
        """Generate detailed PTM section HTML."""
        ptm_summary = data.get('ptm_summary', {})
        ptm_details = ptm_summary.get('details', [])
        peptides = ptm_summary.get('peptides', [])
        ptm_by_type = ptm_summary.get('by_type', {})

        if not ptm_by_type and not peptides:
            return '<p>No PTMs or peptides found</p>'

        ptm_html = '<div class="ptm-container">'

        # PTM Overview
        if ptm_by_type:
            ptm_html += '<div class="ptm-overview">'
            ptm_html += '<h4 class="ptm-subtitle">Modified Residues</h4>'
            ptm_html += '<div class="ptm-type-grid">'
            for ptm_type, count in sorted(ptm_by_type.items()):
                ptm_html += f'<div class="ptm-type-badge"><span class="ptm-type-name">{ptm_type}</span><span class="ptm-count-badge">{count}</span></div>'
            ptm_html += '</div></div>'

        # PTM Details
        if ptm_details:
            ptm_html += '<div class="ptm-details-section">'
            ptm_html += '<h4 class="ptm-subtitle">Modification Sites</h4>'
            ptm_html += '<div class="ptm-list">'

            # Sort by position
            sorted_ptms = sorted(ptm_details, key=lambda x: x.get('position', 0))
            for ptm in sorted_ptms:
                position = ptm.get('position', 0)
                end_position = ptm.get('end_position', 0)
                ptm_type = ptm.get('type', 'Unknown')
                description = ptm.get('description', '')

                pos_text = f"{position}" if position == end_position else f"{position}-{end_position}"

                ptm_html += f'<div class="ptm-item">'
                ptm_html += f'  <div class="ptm-position">Pos {pos_text}</div>'
                ptm_html += f'  <div class="ptm-info">'
                ptm_html += f'    <div class="ptm-type">{ptm_type}</div>'
                if description:
                    ptm_html += f'    <div class="ptm-description">{description}</div>'
                ptm_html += f'  </div>'
                ptm_html += f'</div>'

            ptm_html += '</div></div>'

        ptm_html += '</div>'
        return ptm_html

    def close(self):
        """Close resources."""
        self.data_manager.close()

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python main.py <ACCESSION_ID> [<ACCESSION_ID2> ...]")
        print("Example: python main.py P04637 P15018")
        sys.exit(1)

    accession_ids = [arg.upper() for arg in sys.argv[1:]]

    with ProteinReportGenerator() as generator:
        for accession_id in accession_ids:
            try:
                success = generator.generate_report(accession_id)
                if not success:
                    print(f"⚠️  Failed to generate report for {accession_id}")
            except Exception as e:
                print(f"❌ Error processing {accession_id}: {e}")


if __name__ == '__main__':
    main()
