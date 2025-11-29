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
        properties = analyzer.get_all_properties()
        print(f"   ✅ Calculated properties")

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

        charge_profile_dict = json.loads(charge_profile_json)

        # Handle features list (from UniProt API)
        features_list = features if isinstance(features, list) else (features.get('features', []) if features else [])

        database_data = {
            'protein_name': protein_name,
            'gene_name': gene_name,
            'organism': organism,
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
                'helix_experimental': 0,
                'turn_experimental': 0,
                'sheet_experimental': 0,
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
            'positions': []
        }

        # Handle features as a list (from UniProt API)
        features_list = features if isinstance(features, list) else (features.get('features', []) if isinstance(features, dict) else [])

        # PTM types to look for in UniProt feature types
        ptm_types = {
            'Modified residue', 'Lipidation', 'Glycosylation', 'Phosphorylation',
            'N-linked glycosylation', 'O-linked glycosylation', 'Sulfation',
            'Acetylation', 'Amidation', 'GPI-anchor', 'Ubiquitination',
            'Sumoylation', 'Nitrosylation', 'S-nitrosylation', 'Neddylation',
            'Disulfide bond', 'Proteolytic cleavage', 'Signal peptide',
            'Propeptide', 'Transit peptide'
        }

        for feature in features_list:
            feature_type = feature.get('type', 'Unknown')

            # Check if this is a PTM feature (either by EBI category or UniProt type)
            is_ptm = (feature.get('category') == 'PTM' or feature_type in ptm_types)

            if is_ptm:
                # Get position - handle both UniProt (location.start.value) and EBI (begin) formats
                location = feature.get('location', {})
                if isinstance(location, dict) and 'start' in location:
                    position = location['start'].get('value', 0)
                else:
                    position = feature.get('begin', 0)

                description = feature.get('description', feature_type)

                ptm_summary['total_ptms'] += 1
                ptm_summary['by_type'][feature_type] = ptm_summary['by_type'].get(feature_type, 0) + 1
                ptm_summary['positions'].append({
                    'position': position,
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

    def _enhance_features(self, features: list) -> list:
        """Add tooltips to features."""
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

            tooltip = f"{feature.get('description', feature.get('type'))}\n"
            tooltip += f"Position: {start}-{end} aa\nSize: {size} aa"

            feature_copy = feature.copy()
            feature_copy['start'] = start
            feature_copy['end'] = end
            feature_copy['title'] = tooltip
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
        """Generate ProtParam table HTML."""
        return f"""
        <table class="properties-table">
            <tr><td>Molecular Weight</td><td>{data.get('molecular_weight', 0):.2f} Da</td></tr>
            <tr><td>Isoelectric Point</td><td>{data.get('isoelectric_point', 0):.2f}</td></tr>
            <tr><td>GRAVY Score</td><td>{data.get('gravy_score', 0):.4f}</td></tr>
            <tr><td>Aromaticity</td><td>{data.get('aromaticity', 0):.4f}</td></tr>
            <tr><td>Instability Index</td><td>{data.get('instability_index', 0):.2f}</td></tr>
            <tr><td>Cysteine Count</td><td>{data.get('cysteine_count', 0)}</td></tr>
        </table>
        """

    def _generate_ptm_section(self, data: Dict[str, Any]) -> str:
        """Generate PTM section HTML."""
        ptm_summary = data.get('ptm_summary', {})
        ptm_by_type = ptm_summary.get('by_type', {})

        if not ptm_by_type:
            return '<p>No PTMs found</p>'

        ptm_html = '<div class="ptm-types">'
        for ptm_type, count in ptm_by_type.items():
            ptm_html += f'<div class="ptm-type-item">{ptm_type}: {count}</div>'
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
