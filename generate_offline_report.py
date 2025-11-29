#!/usr/bin/env python3
"""
Generate offline protein analysis reports with Nightingale visualization.

Usage:
    python generate_offline_report.py P15018
    python generate_offline_report.py P04637
    python generate_offline_report.py P04637 P68871

Or edit ACCESSION_TO_PROCESS variable in the script.
"""

import sys
import requests
import json
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from protein_enrichment_service import ProteinEnrichmentService
from protein_database import ProteinDatabase
from uniprot_metadata_extractor import extract_uniprot_metadata
import os

# --- Configuration ---
ACCESSION_TO_PROCESS = "P15018"  # Change this to any UniProt ID

def get_uniprot_data(accession_id):
    """Fetches sequence and feature data from UniProt."""
    sequence_url = f"https://www.uniprot.org/uniprot/{accession_id}.fasta"
    features_url = f"https://www.ebi.ac.uk/proteins/api/features/{accession_id}"
    
    print(f"📥 Fetching sequence from UniProt...")
    seq_response = requests.get(sequence_url)
    if not seq_response.ok:
        print(f"❌ Error: Failed to fetch sequence for {accession_id}")
        print(f"   Status code: {seq_response.status_code}")
        return None, None, None
    
    # Extract sequence from FASTA - remove header and newlines
    fasta_lines = seq_response.text.split("\n")
    sequence = "".join(line.strip() for line in fasta_lines if not line.startswith(">"))
    
    print(f"📥 Fetching features from UniProt...")
    features_response = requests.get(features_url, headers={"Accept": "application/json"})
    if not features_response.ok:
        print(f"⚠️  Warning: Failed to fetch features for {accession_id}")
        print(f"   Status code: {features_response.status_code}")
        print(f"   Continuing without feature data...")
        return sequence, None, None
        
    features_json = features_response.json()
    
    # Get PDB structures
    print(f"📥 Fetching PDB structures...")
    pdb_id = get_pdb_id(accession_id)
    
    return sequence, features_json, pdb_id

def get_pdb_id(accession_id):
    """Fetches the first available PDB structure ID for this protein."""
    try:
        # Query UniProt entry for PDB cross-references
        entry_url = f"https://rest.uniprot.org/uniprotkb/{accession_id}.json"
        response = requests.get(entry_url)
        
        if not response.ok:
            print(f"   ⚠️  Could not fetch PDB structures")
            return None
        
        data = response.json()
        
        # Look for PDB cross-references
        if 'uniProtKBCrossReferences' in data:
            for xref in data['uniProtKBCrossReferences']:
                if xref.get('database') == 'PDB':
                    pdb_id = xref.get('id')
                    print(f"   ✅ Found PDB structure: {pdb_id}")
                    return pdb_id
        
        print(f"   ⚠️  No PDB structures available for this protein")
        return None
        
    except Exception as e:
        print(f"   ⚠️  Error fetching PDB: {e}")
        return None

def extract_ptm_summary(features):
    """
    Extract and summarize post-translational modifications (PTMs) from all PTM-related features.

    Extracts from: MOD_RES, CARBOHYD, DISULFID, LIPID, CROSSLNK, and other PTM types

    Args:
        features: Features JSON from EBI API

    Returns:
        Dictionary with PTM summary data
    """
    ptm_summary = {
        'total_ptms': 0,
        'by_type': {},
        'positions': []
    }

    if not features:
        return ptm_summary

    # Map EBI feature types to PTM categories
    ptm_feature_types = {
        'MOD_RES': 'Modified Residue',
        'CARBOHYD': 'Glycosylation',
        'DISULFID': 'Disulfide Bond',
        'LIPID': 'Lipidation',
        'CROSSLNK': 'Cross-link',
    }

    # PTM type classification mapping for MOD_RES descriptions
    mod_res_type_map = {
        'Phospho': 'Phosphorylation',
        'Acetyl': 'Acetylation',
        'Amidated': 'Amidation',
        'Hydroxyl': 'Hydroxylation',
        'Methylated': 'Methylation',
        'Sulfated': 'Sulfation',
        'Ubiquitinated': 'Ubiquitination',
        'Palmitoyl': 'Palmitoylation',
        'Myristoyl': 'Myristoylation',
        'Hydroxyproline': 'Hydroxylation'
    }

    for feature in features.get('features', []):
        # Use category field as primary indicator of PTM
        category = feature.get('category', '')

        # Only process features explicitly marked as PTM
        if category != 'PTM':
            continue

        feature_type = feature.get('type')
        begin = feature.get('begin')
        description = feature.get('description', '')

        # Determine PTM category based on feature type
        if feature_type == 'MOD_RES':
            # For MOD_RES, analyze the description for specific modification type
            ptm_type = 'Other'
            for key, value in mod_res_type_map.items():
                if key.lower() in description.lower():
                    ptm_type = value
                    break
        elif feature_type in ptm_feature_types:
            # For other known types, use the mapping
            ptm_type = ptm_feature_types[feature_type]
        else:
            # For any other PTM type not in our mapping, use feature type as category
            ptm_type = feature_type

        # Use full description if available, otherwise use feature type name
        display_description = description if description else ptm_feature_types.get(feature_type, feature_type)

        # Increment counts
        ptm_summary['total_ptms'] += 1
        if ptm_type not in ptm_summary['by_type']:
            ptm_summary['by_type'][ptm_type] = 0
        ptm_summary['by_type'][ptm_type] += 1

        # Store position and description
        ptm_summary['positions'].append({
            'position': begin,
            'type': ptm_type,
            'description': display_description,
            'feature_type': feature_type
        })

    # Sort positions by position number for better display
    ptm_summary['positions'].sort(key=lambda x: int(x['position']) if x['position'] else 0)

    return ptm_summary


def generate_ptm_section(ptm_summary):
    """
    Generate HTML section for PTM summary display.

    Args:
        ptm_summary: Dictionary with PTM summary data from extract_ptm_summary()

    Returns:
        HTML string for PTM display section
    """
    if ptm_summary['total_ptms'] == 0:
        return '<div class="ptm-section ptm-empty"><p style="color: var(--gray-600); text-align: center; padding: var(--space-4);">No post-translational modifications identified</p></div>'

    html = '<div class="ptm-section">'
    html += '<div class="ptm-header">'
    html += f'<div class="ptm-title">Post-Translational Modifications</div>'
    html += f'<div class="ptm-count">{ptm_summary["total_ptms"]} total PTMs</div>'
    html += '</div>'

    # PTM types overview
    if ptm_summary['by_type']:
        html += '<div class="ptm-overview">'
        for ptm_type, count in sorted(ptm_summary['by_type'].items(), key=lambda x: x[1], reverse=True):
            html += f'<div class="ptm-type-badge">{ptm_type} <span class="ptm-count-badge">{count}</span></div>'
        html += '</div>'

    # Detailed PTM list (all PTMs with scrollable container)
    html += '<div class="ptm-details">'
    html += '<div class="ptm-list">'

    for ptm in ptm_summary['positions']:
        html += f'''
        <div class="ptm-item">
            <div class="ptm-position">Position {ptm['position']}</div>
            <div class="ptm-info">
                <div class="ptm-type-label">{ptm['type']}</div>
                <div class="ptm-description">{ptm['description']}</div>
            </div>
        </div>
        '''

    html += '</div>'
    html += '</div>'
    html += '</div>'

    return html


def calculate_protparam_table(sequence, disulfide_count=0, features=None):
    """
    Calculates all ProtParam values and returns HTML grid with grouped properties.

    Args:
        sequence: Protein sequence string
        disulfide_count: Number of confirmed disulfide bonds from UniProt
        features: Features JSON from EBI API (for experimental secondary structure)
    """
    try:
        # Clean sequence: replace non-standard amino acids
        clean_sequence = sequence.replace("U", "C").replace("B", "X").replace("Z", "X")
        analysis = ProteinAnalysis(clean_sequence)

        # Calculate all properties
        mw = analysis.molecular_weight() / 1000  # Convert to kDa
        gravy = analysis.gravy()
        aromaticity = analysis.aromaticity() * 100
        instability = analysis.instability_index()
        cys_count = analysis.count_amino_acids().get("C", 0)

        # Predicted secondary structure fractions (from composition)
        sec_struc = analysis.secondary_structure_fraction()
        helix_pred = sec_struc[0] * 100
        turn_pred = sec_struc[1] * 100
        sheet_pred = sec_struc[2] * 100

        # Calculate experimental secondary structure from PDB features
        helix_exp = 0.0
        strand_exp = 0.0
        turn_exp = 0.0

        if features:
            total_aa = len(sequence)
            helix_aa = 0
            strand_aa = 0
            turn_aa = 0

            for feature in features.get('features', []):
                feature_type = feature.get('type')
                begin = int(feature.get('begin', 0))
                end = int(feature.get('end', begin))
                length = end - begin + 1

                if feature_type == 'HELIX':
                    helix_aa += length
                elif feature_type == 'STRAND':
                    strand_aa += length
                elif feature_type == 'TURN':
                    turn_aa += length

            helix_exp = (helix_aa / total_aa) * 100
            strand_exp = (strand_aa / total_aa) * 100
            turn_exp = (turn_aa / total_aa) * 100

        # Molar extinction coefficient (reduced and oxidized forms)
        extinction_coeff = analysis.molar_extinction_coefficient()
        extinction_reduced = extinction_coeff[0]  # Cysteines as -SH
        extinction_oxidized = extinction_coeff[1]  # Disulfide bonds

        # Display properties in three rows
        html = '<div class="properties-grid">'

        # First row: Molecular Weight, GRAVY Score, Aromaticity, Instability Index
        html += '<div class="property-group horizontal-row">'

        stable_status = "Stable" if instability < 40 else "Unstable"

        html += f'''
                <div class="property-item">
                    <div class="property-label">Molecular Weight</div>
                    <div class="property-value">{mw:.2f}<span class="property-unit">kDa</span></div>
                </div>
                <div class="property-item">
                    <div class="property-label">GRAVY Score</div>
                    <div class="property-value">{gravy:.3f}</div>
                </div>
                <div class="property-item">
                    <div class="property-label">Aromaticity</div>
                    <div class="property-value">{aromaticity:.2f}<span class="property-unit">%</span></div>
                </div>
                <div class="property-item">
                    <div class="property-label">Instability Index</div>
                    <div class="property-value">{instability:.2f}<span class="property-unit">{stable_status}</span></div>
                </div>
        '''

        html += '</div>'

        # Second row: Cysteine Count, Disulfide Bonds, Extinction Coefficients
        html += '<div class="property-group horizontal-row">'

        html += f'''
                <div class="property-item">
                    <div class="property-label">Cysteine Count</div>
                    <div class="property-value">{cys_count}<span class="property-unit">residues</span></div>
                </div>
                <div class="property-item">
                    <div class="property-label">Disulfide Bonds</div>
                    <div class="property-value">{disulfide_count}<span class="property-unit">confirmed</span></div>
                </div>
                <div class="property-item">
                    <div class="property-label">Extinction Coeff (Reduced)</div>
                    <div class="property-value">{extinction_reduced}<span class="property-unit">M⁻¹cm⁻¹</span></div>
                </div>
                <div class="property-item">
                    <div class="property-label">Extinction Coeff (Oxidized)</div>
                    <div class="property-value">{extinction_oxidized}<span class="property-unit">M⁻¹cm⁻¹</span></div>
                </div>
        '''

        html += '</div>'

        # Secondary Structure Comparison Table
        html += '<div class="sec-struct-table-wrapper">'
        html += '''
        <table class="sec-struct-comparison">
            <thead>
                <tr>
                    <th>Secondary Structure</th>
                    <th>Predicted (%)</th>
                    <th>Experimental (%)</th>
                </tr>
            </thead>
            <tbody>
        '''

        # Helix row
        exp_helix_display = f"{helix_exp:.2f}" if (features and helix_exp > 0) else "—"
        html += f'''
                <tr>
                    <td>Helix (α)</td>
                    <td class="pred-value">{helix_pred:.2f}</td>
                    <td class="exp-value">{exp_helix_display}</td>
                </tr>
        '''

        # Turn row
        turn_exp_display = f"{turn_exp:.2f}" if (features and turn_exp > 0) else "—"
        html += f'''
                <tr>
                    <td>Turn (β-turn)</td>
                    <td class="pred-value">{turn_pred:.2f}</td>
                    <td class="exp-value">{turn_exp_display}</td>
                </tr>
        '''

        # Sheet row
        exp_sheet_display = f"{strand_exp:.2f}" if (features and strand_exp > 0) else "—"
        html += f'''
                <tr>
                    <td>Sheet (β)</td>
                    <td class="pred-value">{sheet_pred:.2f}</td>
                    <td class="exp-value">{exp_sheet_display}</td>
                </tr>
        '''

        html += '''
            </tbody>
        </table>
        </div>
        '''

        return html

    except Exception as e:
        print(f"❌ Error in ProtParam calculation: {e}")
        return "<div style='color: red; padding: 1rem;'>Error calculating ProtParam data.</div>"

def get_uniprot_disulfide_bond_count(accession_id):
    """
    Get the count of confirmed disulfide bonds from UniProt REST API.

    Args:
        accession_id: UniProt accession ID

    Returns:
        Number of confirmed disulfide bonds (int), 0 if not available
    """
    try:
        entry_url = f"https://rest.uniprot.org/uniprotkb/{accession_id}.json"
        response = requests.get(entry_url)

        if not response.ok:
            return 0

        data = response.json()

        # Extract disulfide bond count from extraAttributes
        if 'extraAttributes' in data:
            count_by_feature = data['extraAttributes'].get('countByFeatureType', {})
            return count_by_feature.get('Disulfide bond', 0)

        return 0

    except Exception as e:
        print(f"   ⚠️  Could not fetch disulfide bond count: {e}")
        return 0


def generate_charge_profile_json(profile):
    """
    Generate charge profile JSON data for template injection.

    Args:
        profile: ProteinEnrichmentProfile object

    Returns:
        JSON string with charge profile data
    """
    data = {
        "ph_values": profile.charges.ph_values,
        "charges": profile.charges.charges,
        "isoelectric_point": profile.charges.isoelectric_point,
        "charge_at_ph_7_4": profile.charges.charge_at_ph_7_4
    }
    return json.dumps(data)


def enhance_features_with_tooltips(features_list):
    """
    Enhance features with tooltips containing basic information.

    This function:
    1. Converts 'begin' and 'end' to 'start' and 'end'
    2. Creates tooltip text with feature information

    Args:
        features_list: List of feature dictionaries from EBI API

    Returns:
        Enhanced features list with 'title' field added
    """
    # Convert begin/end to start/end and build tooltips
    enhanced_features = []
    for feature in features_list:
        start = int(feature.get('begin') or feature.get('start', 0))
        end = int(feature.get('end') or start)
        size = end - start + 1

        ft_type = feature.get('type', 'Unknown')
        description = feature.get('description', ft_type)
        ft_id = feature.get('ftId', 'N/A')

        # Build basic tooltip
        tooltip_text = f"{description}\nPosition: {start}-{end} aa\nSize: {size} aa"
        tooltip_text += f"\nFeature ID: {ft_id}"

        enhanced_feature = feature.copy()
        enhanced_feature['start'] = start
        enhanced_feature['end'] = end
        enhanced_feature['title'] = tooltip_text
        enhanced_features.append(enhanced_feature)

    # Return all enhanced features with tooltips specific to each feature
    return enhanced_features


def process_protein(accession_id):
    """
    Process a single protein and generate offline HTML report.

    Args:
        accession_id: UniProt accession ID
    """
    print(f"\n🔬 Processing: {accession_id}")

    # 1. Fetch UniProt data
    sequence, features, pdb_id = get_uniprot_data(accession_id)
    if not sequence:
        print(f"❌ Failed to fetch sequence for {accession_id}. Skipping.")
        return

    print(f"   ✅ Fetched sequence ({len(sequence)} amino acids)")

    if not features:
        print("   ⚠️  Continuing without feature data...")
        features = {"features": []}
    else:
        print(f"   ✅ Fetched {len(features.get('features', []))} features")

    # 2. Get disulfide bond count from UniProt
    print(f"   🔗 Fetching confirmed disulfide bonds from UniProt...")
    disulfide_count = get_uniprot_disulfide_bond_count(accession_id)
    if disulfide_count > 0:
        print(f"   ✅ Found {disulfide_count} confirmed disulfide bond(s)")

    # 3. Extract PTM summary
    print(f"   🔧 Analyzing post-translational modifications...")
    ptm_summary = extract_ptm_summary(features)
    ptm_html = generate_ptm_section(ptm_summary)
    if ptm_summary['total_ptms'] > 0:
        print(f"   ✅ Found {ptm_summary['total_ptms']} post-translational modification(s)")

    # 4. Calculate ProtParam properties
    print(f"   🧮 Calculating biophysical properties...")
    protparam_table_html = calculate_protparam_table(sequence, disulfide_count, features)

    # 4. Generate enriched profile with charge data using service
    print(f"   🔋 Generating charge profile...")
    try:
        service = ProteinEnrichmentService()
        profile = service.enrich_protein(accession_id)
        charge_profile_json = generate_charge_profile_json(profile)
        print(f"   ✅ Generated charge profile (pH 2.0 - 12.0)")
    except Exception as e:
        print(f"   ⚠️  Could not generate charge profile: {e}")
        charge_profile_json = json.dumps({
            "ph_values": [],
            "charges": [],
            "isoelectric_point": 0,
            "charge_at_ph_7_4": 0
        })

    # 4. Read the HTML template
    template_file = "offline_template.html"
    try:
        with open(template_file, "r", encoding="utf-8") as f:
            template_content = f.read()
    except FileNotFoundError:
        print(f"\n❌ Error: {template_file} not found in current directory.")
        print(f"   Make sure the template file exists!")
        return

    # 5. Inject data into template
    print(f"   📝 Rendering template...")

    # Prepare data for injection
    # Enhance features with tooltips and overlap detection before JSON serialization
    raw_features = features.get('features', [])
    enhanced_features = enhance_features_with_tooltips(raw_features)
    features_json_string = json.dumps(enhanced_features)
    sequence_length = str(len(sequence))
    pdb_id_value = pdb_id if pdb_id else ""
    has_structure = "true" if pdb_id else "false"

    # Get protein metadata from general info
    protein_name = profile.general.protein_name if profile.general.protein_name else "Unknown"
    gene_name = profile.general.gene_name if profile.general.gene_name else "N/A"
    organism = profile.general.organism if profile.general.organism else "Unknown"

    report_content = template_content
    report_content = report_content.replace("{{ACCESSION_ID}}", accession_id)
    report_content = report_content.replace("{{PROTEIN_NAME}}", protein_name)
    report_content = report_content.replace("{{GENE_NAME}}", gene_name)
    report_content = report_content.replace("{{ORGANISM}}", organism)
    report_content = report_content.replace("{{PTM_SECTION}}", ptm_html)
    report_content = report_content.replace("{{PROTPARAM_TABLE}}", protparam_table_html)
    report_content = report_content.replace("{{SEQUENCE}}", sequence)
    report_content = report_content.replace("{{SEQUENCE_LENGTH}}", sequence_length)
    report_content = report_content.replace("{{FEATURES_JSON}}", features_json_string)
    report_content = report_content.replace("{{PDB_ID}}", pdb_id_value)
    report_content = report_content.replace("{{HAS_STRUCTURE}}", has_structure)
    report_content = report_content.replace("{{CHARGE_PROFILE_JSON}}", charge_profile_json)

    # 6. Save data to database
    print(f"   💾 Saving data to database...")
    protein_id = None
    try:
        db = ProteinDatabase()

        # Recalculate biophysical properties for database storage
        clean_sequence = sequence.replace("U", "C").replace("B", "X").replace("Z", "X")
        analysis = ProteinAnalysis(clean_sequence)

        mw = analysis.molecular_weight() / 1000  # Convert to kDa
        gravy = analysis.gravy()
        aromaticity = analysis.aromaticity() * 100
        instability = analysis.instability_index()
        cys_count = analysis.count_amino_acids().get("C", 0)
        extinction_coeff = analysis.molar_extinction_coefficient()
        extinction_reduced = extinction_coeff[0]
        extinction_oxidized = extinction_coeff[1]

        # Predicted secondary structure fractions
        sec_struc_frac = analysis.secondary_structure_fraction()
        helix_pred = sec_struc_frac[0] * 100
        turn_pred = sec_struc_frac[1] * 100
        sheet_pred = sec_struc_frac[2] * 100

        # Experimental secondary structure from features
        helix_exp = 0.0
        strand_exp = 0.0
        turn_exp = 0.0

        if features and 'features' in features:
            total_aa = len(sequence)
            helix_aa = 0
            strand_aa = 0
            turn_aa = 0

            for feature in features.get('features', []):
                feature_type = feature.get('type')
                begin = int(feature.get('begin', 0))
                end = int(feature.get('end', begin))
                length = end - begin + 1

                if feature_type == 'HELIX':
                    helix_aa += length
                elif feature_type == 'STRAND':
                    strand_aa += length
                elif feature_type == 'TURN':
                    turn_aa += length

            if total_aa > 0:
                helix_exp = (helix_aa / total_aa) * 100
                strand_exp = (strand_aa / total_aa) * 100
                turn_exp = (turn_aa / total_aa) * 100

        # Extract charge profile data
        charge_profile_dict = json.loads(charge_profile_json)

        # Extract comprehensive UniProt metadata
        print(f"   📚 Extracting UniProt metadata...")
        uniprot_metadata = extract_uniprot_metadata(accession_id)
        if uniprot_metadata:
            print(f"   ✅ Extracted UniProt metadata")
        else:
            print(f"   ⚠️  Could not extract UniProt metadata")
            uniprot_metadata = None

        # Prepare protein data for storage
        protein_data = {
            'protein_name': protein_name,
            'gene_name': gene_name,
            'organism': organism,
            'sequence': sequence,
            'sequence_length': len(sequence),
            'molecular_weight': mw,
            'isoelectric_point': charge_profile_dict.get('isoelectric_point'),
            'gravy_score': gravy,
            'aromaticity': aromaticity,
            'instability_index': instability,
            'cysteine_count': cys_count,
            'extinction_coeff_reduced': extinction_reduced,
            'extinction_coeff_oxidized': extinction_oxidized,
            'pdb_id': pdb_id,
            'features': enhanced_features,
            'ptm_summary': ptm_summary,
            'secondary_structure': {
                'helix_predicted': helix_pred,
                'turn_predicted': turn_pred,
                'sheet_predicted': sheet_pred,
                'helix_experimental': helix_exp,
                'turn_experimental': turn_exp,
                'sheet_experimental': strand_exp
            },
            'charge_profile': charge_profile_dict,
            'disulfide_bonds': [],
            'uniprot_metadata': uniprot_metadata
        }

        # Save to database
        protein_id = db.save_protein(accession_id, protein_data)
        print(f"   ✅ Data saved to database (ID: {protein_id})")
        db.close()
    except Exception as e:
        print(f"   ⚠️  Warning: Could not save to database: {e}")
        protein_id = None

    # 7. Write the final report
    output_filename = f"report_{accession_id}.html"
    with open(output_filename, "w", encoding="utf-8") as f:
        f.write(report_content)

    print(f"   ✅ Report saved: {output_filename}")

    # Save report metadata to database
    try:
        db = ProteinDatabase()
        file_size_kb = os.path.getsize(output_filename) / 1024
        db.save_report(protein_id, output_filename, os.path.abspath(output_filename), file_size_kb)
        db.close()
    except Exception as e:
        print(f"   ⚠️  Warning: Could not save report metadata: {e}")

    return output_filename


def main():
    # Check if accessions provided as command-line arguments
    if len(sys.argv) > 1:
        accession_ids = [arg.upper() for arg in sys.argv[1:]]  # Convert to uppercase
    else:
        print("Usage: python generate_offline_report.py <UNIPROT_ID> [<UNIPROT_ID2> ...]")
        print("\nExample:")
        print("  python generate_offline_report.py P15018")
        print("  python generate_offline_report.py P04637 P68871 P15018")
        print("\nOr edit ACCESSION_TO_PROCESS in the script to set a default.")
        print(f"\nUsing default: {ACCESSION_TO_PROCESS}")
        accession_ids = [ACCESSION_TO_PROCESS]

    print("=" * 60)
    print(f"  OFFLINE PROTEIN REPORT GENERATOR")
    print("=" * 60)

    output_files = []
    for accession_id in accession_ids:
        output_file = process_protein(accession_id)
        if output_file:
            output_files.append(output_file)

    if output_files:
        print("\n" + "=" * 60)
        print(f"📄 Generated {len(output_files)} report(s)")
        print("=" * 60)
        print("\n📖 USAGE INSTRUCTIONS:")
        print(f"\n  1. Start a local web server:")
        print(f"     python -m http.server 8000")
        print(f"\n  2. Open in your browser (examples):")
        for output_file in output_files:
            print(f"     http://localhost:8000/{output_file}")
        print(f"\n  3. View the interactive protein report with charge profile!")
        print()
    else:
        print("\n❌ No reports were generated successfully.")

if __name__ == "__main__":
    main()