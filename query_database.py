#!/usr/bin/env python3
"""
Database Query Utility
View and query stored protein data
"""

import sys
from protein_database import ProteinDatabase
import json


def print_separator(title=""):
    """Print a formatted separator"""
    if title:
        print(f"\n{'='*70}")
        print(f"  {title}")
        print(f"{'='*70}")
    else:
        print(f"{'='*70}")


def list_proteins():
    """List all stored proteins"""
    db = ProteinDatabase()
    proteins = db.get_all_proteins(limit=50)

    if not proteins:
        print("No proteins stored in database yet.")
        db.close()
        return

    print_separator("STORED PROTEINS")
    print(f"{'Accession':<15} {'Name':<30} {'Organism':<20} {'Length':<10}")
    print("-" * 75)

    for protein in proteins:
        accession, name, organism, length = protein
        name_short = (name[:27] + "...") if name and len(name) > 30 else (name or "Unknown")
        org_short = (organism[:17] + "...") if organism and len(organism) > 20 else (organism or "Unknown")
        print(f"{accession:<15} {name_short:<30} {org_short:<20} {length:<10}")

    print(f"\nTotal proteins: {len(proteins)}")
    db.close()


def view_protein(accession_id):
    """View detailed information about a protein"""
    db = ProteinDatabase()
    protein = db.get_protein(accession_id)

    if not protein:
        print(f"❌ Protein {accession_id} not found in database.")
        db.close()
        return

    print_separator(f"PROTEIN: {accession_id}")

    # Basic info
    print("\n📋 BASIC INFORMATION:")
    print(f"  Name: {protein['protein_name']}")
    print(f"  Gene: {protein['gene_name']}")
    print(f"  Organism: {protein['organism']}")
    print(f"  Sequence Length: {protein['sequence_length']} aa")
    print(f"  PDB ID: {protein['pdb_id'] or 'N/A'}")

    # Biophysical properties
    print("\n⚗️  BIOPHYSICAL PROPERTIES:")
    print(f"  Molecular Weight: {protein['molecular_weight']:.2f} kDa" if protein['molecular_weight'] else "  Molecular Weight: N/A")
    print(f"  Isoelectric Point: {protein['isoelectric_point']:.2f}" if protein['isoelectric_point'] else "  Isoelectric Point: N/A")
    print(f"  GRAVY Score: {protein['gravy_score']:.3f}" if protein['gravy_score'] else "  GRAVY Score: N/A")
    print(f"  Aromaticity: {protein['aromaticity']:.2f}%" if protein['aromaticity'] else "  Aromaticity: N/A")
    print(f"  Instability Index: {protein['instability_index']:.2f}" if protein['instability_index'] else "  Instability Index: N/A")
    print(f"  Cysteine Count: {protein['cysteine_count']}")
    print(f"  Extinction Coeff (Reduced): {protein['extinction_coeff_reduced']} M⁻¹cm⁻¹")
    print(f"  Extinction Coeff (Oxidized): {protein['extinction_coeff_oxidized']} M⁻¹cm⁻¹")

    # Features
    protein_id = protein['protein_id']
    features = db.get_protein_features(protein_id)
    print(f"\n🧬 FEATURES: {len(features)} total")
    feature_types = {}
    for feature in features:
        ftype = feature['feature_type']
        feature_types[ftype] = feature_types.get(ftype, 0) + 1
    for ftype, count in sorted(feature_types.items()):
        print(f"  {ftype}: {count}")

    # PTMs
    ptm_details = db.get_ptm_details(protein_id)
    print(f"\n⚡ PTMs: {len(ptm_details)} total")
    if ptm_details:
        ptm_types = {}
        for ptm in ptm_details:
            ptype = ptm['ptm_type']
            ptm_types[ptype] = ptm_types.get(ptype, 0) + 1
        for ptype, count in sorted(ptm_types.items()):
            print(f"  {ptype}: {count}")

    # Charge profile
    charge_data = db.get_charge_profile(protein_id)
    if charge_data['ph_values']:
        print(f"\n📈 CHARGE PROFILE: {len(charge_data['ph_values'])} pH points")
        print(f"  pH Range: {charge_data['ph_values'][0]:.1f} - {charge_data['ph_values'][-1]:.1f}")
        print(f"  Charge Range: {min(charge_data['charges']):.2f} to {max(charge_data['charges']):.2f}")

    # Reports
    reports = db.get_reports(protein_id)
    print(f"\n📄 GENERATED REPORTS: {len(reports)}")
    for report in reports:
        filename, path, size, generated = report
        print(f"  {filename} ({size:.1f} KB) - {generated}")

    db.close()


def search_proteins(search_term):
    """Search proteins by name, gene, or accession"""
    db = ProteinDatabase()
    results = db.search_proteins(search_term)

    if not results:
        print(f"No proteins found matching '{search_term}'")
        db.close()
        return

    print_separator(f"SEARCH RESULTS FOR '{search_term}'")
    print(f"{'Accession':<15} {'Name':<30} {'Gene':<15} {'Length':<10}")
    print("-" * 70)

    for protein in results:
        accession, name, gene, organism, length = protein
        name_short = (name[:27] + "...") if name and len(name) > 30 else (name or "Unknown")
        gene_short = (gene[:12] + "...") if gene and len(gene) > 15 else (gene or "N/A")
        print(f"{accession:<15} {name_short:<30} {gene_short:<15} {length:<10}")

    print(f"\nFound {len(results)} protein(s)")
    db.close()


def get_statistics():
    """Show database statistics"""
    db = ProteinDatabase()
    stats = db.get_statistics()

    print_separator("DATABASE STATISTICS")
    print(f"  Total Proteins: {stats['total_proteins']}")
    print(f"  Total Features: {stats['total_features']}")
    print(f"  Total PTMs: {stats['total_ptms']}")
    print(f"  Total Reports Generated: {stats['total_reports']}")
    print(f"  Database Size: {stats['db_size_mb']} MB")

    db.close()


def print_usage():
    """Print usage information"""
    print("""
Database Query Tool Usage:

  python query_database.py list              # List all stored proteins
  python query_database.py view <ACCESSION>  # View protein details
  python query_database.py search <TERM>     # Search proteins
  python query_database.py stats              # Show database statistics

Examples:
  python query_database.py list
  python query_database.py view P04637
  python query_database.py search hemoglobin
  python query_database.py stats
    """)


def main():
    """Main entry point"""
    if len(sys.argv) < 2:
        print_usage()
        return

    command = sys.argv[1].lower()

    if command == "list":
        list_proteins()
    elif command == "view":
        if len(sys.argv) < 3:
            print("❌ Please provide accession ID: python query_database.py view <ACCESSION>")
            return
        view_protein(sys.argv[2].upper())
    elif command == "search":
        if len(sys.argv) < 3:
            print("❌ Please provide search term: python query_database.py search <TERM>")
            return
        search_proteins(" ".join(sys.argv[2:]))
    elif command == "stats":
        get_statistics()
    else:
        print(f"❌ Unknown command: {command}")
        print_usage()


if __name__ == "__main__":
    main()
