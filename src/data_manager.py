#!/usr/bin/env python3
"""
Data manager module.
Handles database operations and data persistence.
"""

from typing import Dict, Any, Optional, List
from src.services.protein_database import ProteinDatabase
from src.services.uniprot_metadata_extractor import extract_uniprot_metadata


class DataManager:
    """Manages data persistence and retrieval."""

    def __init__(self, db_path: str = 'protein_reports.db'):
        """
        Initialize data manager.

        Args:
            db_path: Path to SQLite database
        """
        self.db = ProteinDatabase(db_path)

    def save_protein(self, accession_id: str, protein_data: Dict[str, Any]) -> Optional[int]:
        """
        Save complete protein data to database.

        Args:
            accession_id: UniProt accession ID
            protein_data: Dictionary with all protein information

        Returns:
            Protein ID from database or None if failed
        """
        try:
            # Extract metadata from UniProt
            print(f"   📚 Extracting UniProt metadata...")
            metadata = extract_uniprot_metadata(accession_id)
            if metadata:
                print(f"   ✅ Extracted UniProt metadata")
                protein_data['uniprot_metadata'] = metadata
            else:
                print(f"   ⚠️  Could not extract UniProt metadata")

            # Save to database
            protein_id = self.db.save_protein(accession_id, protein_data)
            if protein_id:
                print(f"   ✅ Data saved to database (ID: {protein_id})")
                return protein_id
            else:
                print(f"   ⚠️  Failed to save data to database")
                return None

        except Exception as e:
            print(f"   ⚠️  Error saving to database: {e}")
            return None

    def save_report(self, protein_id: int, report_filename: str,
                   report_path: str, file_size_kb: float) -> bool:
        """
        Save report metadata to database.

        Args:
            protein_id: Database protein ID
            report_filename: Name of report file
            report_path: Full path to report file
            file_size_kb: Report file size in KB

        Returns:
            True if successful, False otherwise
        """
        try:
            self.db.save_report(protein_id, report_filename, report_path, file_size_kb)
            return True
        except Exception as e:
            print(f"   ⚠️  Error saving report metadata: {e}")
            return False

    def get_protein(self, accession_id: str) -> Optional[Dict[str, Any]]:
        """Get protein data by accession ID."""
        return self.db.get_protein(accession_id)

    def list_proteins(self, limit: int = None) -> List:
        """Get list of all stored proteins."""
        return self.db.get_all_proteins(limit)

    def search_proteins(self, search_term: str) -> List:
        """Search proteins by name, gene, or accession."""
        return self.db.search_proteins(search_term)

    def get_statistics(self) -> Dict[str, Any]:
        """Get database statistics."""
        return self.db.get_statistics()

    def close(self):
        """Close database connection."""
        self.db.close()

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()
