#!/usr/bin/env python3
"""
Protein Database Management Module
Handles storing and retrieving protein data from SQLite database
"""

import sqlite3
import json
from datetime import datetime
from pathlib import Path


class ProteinDatabase:
    """Database manager for storing and retrieving protein data"""

    def __init__(self, db_path="protein_reports.db"):
        """Initialize database connection"""
        self.db_path = db_path
        self.connection = None
        self.initialize_db()

    def initialize_db(self):
        """Initialize database with schema"""
        self.connection = sqlite3.connect(self.db_path)
        self.connection.row_factory = sqlite3.Row

        # Read and execute schema
        schema_file = Path(__file__).parent / "database_schema.sql"
        if schema_file.exists():
            with open(schema_file, "r") as f:
                schema = f.read()
            self.connection.executescript(schema)
            self.connection.commit()
            print(f"✓ Database initialized: {self.db_path}")
        else:
            print(f"⚠️  Schema file not found: {schema_file}")

    def close(self):
        """Close database connection"""
        if self.connection:
            self.connection.close()

    def protein_exists(self, accession_id):
        """Check if protein already exists in database"""
        cursor = self.connection.cursor()
        cursor.execute("SELECT protein_id FROM proteins WHERE accession_id = ?", (accession_id,))
        return cursor.fetchone() is not None

    def save_protein(self, accession_id, protein_data):
        """
        Save complete protein data to database

        Args:
            accession_id: UniProt accession ID
            protein_data: Dictionary containing all protein information
        """
        cursor = self.connection.cursor()

        try:
            # Insert or update main protein record
            cursor.execute("""
                INSERT OR REPLACE INTO proteins (
                    accession_id, protein_name, gene_name, organism, sequence,
                    sequence_length, molecular_weight, isoelectric_point,
                    gravy_score, aromaticity, instability_index, cysteine_count,
                    extinction_coeff_reduced, extinction_coeff_oxidized, pdb_id, has_structure
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                accession_id,
                protein_data.get('protein_name'),
                protein_data.get('gene_name'),
                protein_data.get('organism'),
                protein_data.get('sequence', ''),
                protein_data.get('sequence_length', 0),
                protein_data.get('molecular_weight'),
                protein_data.get('isoelectric_point'),
                protein_data.get('gravy_score'),
                protein_data.get('aromaticity'),
                protein_data.get('instability_index'),
                protein_data.get('cysteine_count'),
                protein_data.get('extinction_coeff_reduced'),
                protein_data.get('extinction_coeff_oxidized'),
                protein_data.get('pdb_id'),
                1 if protein_data.get('pdb_id') else 0
            ))

            protein_id = cursor.lastrowid

            # Save secondary structure data
            if 'secondary_structure' in protein_data:
                ss = protein_data['secondary_structure']
                cursor.execute("""
                    INSERT OR REPLACE INTO secondary_structure
                    (protein_id, helix_predicted, turn_predicted, sheet_predicted,
                     helix_experimental, turn_experimental, sheet_experimental)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                """, (
                    protein_id,
                    ss.get('helix_predicted'),
                    ss.get('turn_predicted'),
                    ss.get('sheet_predicted'),
                    ss.get('helix_experimental'),
                    ss.get('turn_experimental'),
                    ss.get('sheet_experimental')
                ))

            # Save charge profile data
            if 'charge_profile' in protein_data:
                cursor.execute("DELETE FROM charge_profile WHERE protein_id = ?", (protein_id,))
                for ph, charge in zip(protein_data['charge_profile']['ph_values'],
                                     protein_data['charge_profile']['charges']):
                    cursor.execute("""
                        INSERT INTO charge_profile (protein_id, ph_value, charge)
                        VALUES (?, ?, ?)
                    """, (protein_id, ph, charge))

            # Save PTM summary
            if 'ptm_summary' in protein_data:
                ptm = protein_data['ptm_summary']
                cursor.execute("DELETE FROM ptm_summary WHERE protein_id = ?", (protein_id,))
                cursor.execute("""
                    INSERT INTO ptm_summary (protein_id, total_ptms)
                    VALUES (?, ?)
                """, (protein_id, ptm.get('total_ptms', 0)))

                ptm_id = cursor.lastrowid

                # Save PTM type counts
                if 'by_type' in ptm:
                    for ptm_type, count in ptm['by_type'].items():
                        cursor.execute("""
                            INSERT INTO ptm_type_counts (ptm_id, ptm_type, count)
                            VALUES (?, ?, ?)
                        """, (ptm_id, ptm_type, count))

                # Save individual PTM details
                if 'positions' in ptm:
                    for ptm_detail in ptm['positions']:
                        cursor.execute("""
                            INSERT INTO ptm_details (ptm_id, position, ptm_type, description, feature_type)
                            VALUES (?, ?, ?, ?, ?)
                        """, (
                            ptm_id,
                            ptm_detail.get('position'),
                            ptm_detail.get('type'),
                            ptm_detail.get('description'),
                            ptm_detail.get('feature_type')
                        ))

            # Save protein features
            if 'features' in protein_data:
                cursor.execute("DELETE FROM protein_features WHERE protein_id = ?", (protein_id,))
                for feature in protein_data['features']:
                    start = int(feature.get('start') or feature.get('begin', 0))
                    end = int(feature.get('end') or start)
                    length = end - start + 1

                    cursor.execute("""
                        INSERT INTO protein_features
                        (protein_id, feature_type, category, description, start_position,
                         end_position, feature_length, ft_id)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """, (
                        protein_id,
                        feature.get('type'),
                        feature.get('category'),
                        feature.get('description'),
                        start,
                        end,
                        length,
                        feature.get('ftId')
                    ))

            # Save disulfide bonds
            if 'disulfide_bonds' in protein_data:
                cursor.execute("DELETE FROM disulfide_bonds WHERE protein_id = ?", (protein_id,))
                for bond in protein_data['disulfide_bonds']:
                    cursor.execute("""
                        INSERT INTO disulfide_bonds (protein_id, position_1, position_2)
                        VALUES (?, ?, ?)
                    """, (protein_id, bond.get('pos1'), bond.get('pos2')))

            # Save UniProt metadata if provided
            if 'uniprot_metadata' in protein_data:
                self.save_uniprot_metadata(protein_id, protein_data['uniprot_metadata'])

            self.connection.commit()
            return protein_id

        except Exception as e:
            self.connection.rollback()
            print(f"❌ Error saving protein data: {e}")
            raise

    def save_report(self, protein_id, report_filename, report_path, file_size_kb):
        """Save generated report metadata"""
        cursor = self.connection.cursor()
        cursor.execute("""
            INSERT INTO generated_reports (protein_id, report_filename, report_path, file_size_kb)
            VALUES (?, ?, ?, ?)
        """, (protein_id, report_filename, report_path, file_size_kb))
        self.connection.commit()

    def get_protein(self, accession_id):
        """Retrieve protein data by accession ID"""
        cursor = self.connection.cursor()
        cursor.execute("SELECT * FROM proteins WHERE accession_id = ?", (accession_id,))
        return cursor.fetchone()

    def get_all_proteins(self, limit=None):
        """Get list of all stored proteins"""
        cursor = self.connection.cursor()
        query = "SELECT accession_id, protein_name, organism, sequence_length FROM proteins ORDER BY created_at DESC"
        if limit:
            query += f" LIMIT {limit}"
        cursor.execute(query)
        return cursor.fetchall()

    def get_protein_features(self, protein_id):
        """Get all features for a protein"""
        cursor = self.connection.cursor()
        cursor.execute("""
            SELECT feature_type, category, description, start_position, end_position
            FROM protein_features
            WHERE protein_id = ?
            ORDER BY start_position
        """, (protein_id,))
        return cursor.fetchall()

    def get_ptm_details(self, protein_id):
        """Get PTM details for a protein"""
        cursor = self.connection.cursor()
        cursor.execute("""
            SELECT pd.position, pd.ptm_type, pd.description, COUNT(*) as count
            FROM ptm_details pd
            JOIN ptm_summary ps ON pd.ptm_id = ps.ptm_id
            WHERE ps.protein_id = ?
            GROUP BY pd.position, pd.ptm_type
            ORDER BY pd.position
        """, (protein_id,))
        return cursor.fetchall()

    def get_charge_profile(self, protein_id):
        """Get charge profile data for a protein"""
        cursor = self.connection.cursor()
        cursor.execute("""
            SELECT ph_value, charge FROM charge_profile
            WHERE protein_id = ?
            ORDER BY ph_value
        """, (protein_id,))
        rows = cursor.fetchall()
        return {
            'ph_values': [row[0] for row in rows],
            'charges': [row[1] for row in rows]
        }

    def get_reports(self, protein_id):
        """Get all reports generated for a protein"""
        cursor = self.connection.cursor()
        cursor.execute("""
            SELECT report_filename, report_path, file_size_kb, generated_at
            FROM generated_reports
            WHERE protein_id = ?
            ORDER BY generated_at DESC
        """, (protein_id,))
        return cursor.fetchall()

    def search_proteins(self, search_term):
        """Search proteins by name, gene, or accession"""
        cursor = self.connection.cursor()
        search_pattern = f"%{search_term}%"
        cursor.execute("""
            SELECT accession_id, protein_name, gene_name, organism, sequence_length
            FROM proteins
            WHERE accession_id LIKE ? OR protein_name LIKE ? OR gene_name LIKE ?
            ORDER BY created_at DESC
        """, (search_pattern, search_pattern, search_pattern))
        return cursor.fetchall()

    def get_statistics(self):
        """Get database statistics"""
        cursor = self.connection.cursor()

        stats = {}

        # Total proteins
        cursor.execute("SELECT COUNT(*) FROM proteins")
        stats['total_proteins'] = cursor.fetchone()[0]

        # Total features
        cursor.execute("SELECT COUNT(*) FROM protein_features")
        stats['total_features'] = cursor.fetchone()[0]

        # Total PTMs
        cursor.execute("SELECT SUM(total_ptms) FROM ptm_summary")
        result = cursor.fetchone()[0]
        stats['total_ptms'] = result if result else 0

        # Total reports generated
        cursor.execute("SELECT COUNT(*) FROM generated_reports")
        stats['total_reports'] = cursor.fetchone()[0]

        # Database size
        db_size = Path(self.db_path).stat().st_size / (1024 * 1024)  # MB
        stats['db_size_mb'] = round(db_size, 2)

        return stats

    def save_uniprot_metadata(self, protein_id, metadata):
        """
        Save comprehensive UniProt metadata for a protein

        Args:
            protein_id: Database protein ID
            metadata: Dictionary containing UniProt metadata (secondary accessions, genes, names, comments, etc.)
        """
        cursor = self.connection.cursor()

        try:
            # Save secondary accessions
            if 'secondary_accessions' in metadata:
                cursor.execute("DELETE FROM protein_secondary_accessions WHERE protein_id = ?", (protein_id,))
                for sec_acc in metadata['secondary_accessions']:
                    cursor.execute("""
                        INSERT INTO protein_secondary_accessions (protein_id, secondary_accession)
                        VALUES (?, ?)
                    """, (protein_id, sec_acc))

            # Save gene names
            if 'genes' in metadata:
                cursor.execute("DELETE FROM protein_genes WHERE protein_id = ?", (protein_id,))
                for gene in metadata['genes']:
                    cursor.execute("""
                        INSERT INTO protein_genes (protein_id, gene_name, gene_type)
                        VALUES (?, ?, ?)
                    """, (protein_id, gene.get('gene_name'), gene.get('gene_type', 'primary')))

            # Save protein names
            if 'protein_names' in metadata:
                cursor.execute("DELETE FROM protein_names WHERE protein_id = ?", (protein_id,))
                for name in metadata['protein_names']:
                    cursor.execute("""
                        INSERT INTO protein_names (protein_id, name_type, name_value)
                        VALUES (?, ?, ?)
                    """, (protein_id, name.get('name_type'), name.get('name_value')))

            # Save comments
            if 'comments' in metadata:
                cursor.execute("DELETE FROM protein_comments WHERE protein_id = ?", (protein_id,))
                for comment in metadata['comments']:
                    cursor.execute("""
                        INSERT INTO protein_comments (protein_id, comment_type, comment_text)
                        VALUES (?, ?, ?)
                    """, (protein_id, comment.get('comment_type'), comment.get('comment_text')))

            # Save cross-references
            if 'cross_references' in metadata:
                cursor.execute("DELETE FROM protein_cross_references WHERE protein_id = ?", (protein_id,))
                for xref in metadata['cross_references']:
                    cursor.execute("""
                        INSERT INTO protein_cross_references (protein_id, database_name, database_id, properties)
                        VALUES (?, ?, ?, ?)
                    """, (protein_id, xref.get('database_name'), xref.get('database_id'),
                          json.dumps(xref.get('properties', {}))))

            # Save keywords
            if 'keywords' in metadata:
                cursor.execute("DELETE FROM protein_keywords WHERE protein_id = ?", (protein_id,))
                for keyword in metadata['keywords']:
                    cursor.execute("""
                        INSERT INTO protein_keywords (protein_id, keyword_id_code, keyword_name)
                        VALUES (?, ?, ?)
                    """, (protein_id, keyword.get('keyword_id'), keyword.get('keyword_name')))

            # Save GO terms
            if 'go_terms' in metadata:
                cursor.execute("DELETE FROM protein_go_terms WHERE protein_id = ?", (protein_id,))
                for go in metadata['go_terms']:
                    cursor.execute("""
                        INSERT INTO protein_go_terms (protein_id, go_id, go_term, go_aspect, evidence_code)
                        VALUES (?, ?, ?, ?, ?)
                    """, (protein_id, go.get('go_id'), go.get('go_term'), go.get('go_aspect'), go.get('evidence_code')))

            # Save interactions
            if 'interactions' in metadata:
                cursor.execute("DELETE FROM protein_interactions WHERE protein_id = ?", (protein_id,))
                for interaction in metadata['interactions']:
                    cursor.execute("""
                        INSERT INTO protein_interactions (protein_id, interactant_accession, interactant_name,
                                                        interaction_type, experiments)
                        VALUES (?, ?, ?, ?, ?)
                    """, (protein_id, interaction.get('interactant_accession'), interaction.get('interactant_name'),
                          interaction.get('interaction_type'), interaction.get('experiments')))

            # Save references
            if 'references' in metadata:
                cursor.execute("DELETE FROM protein_references WHERE protein_id = ?", (protein_id,))
                for ref in metadata['references']:
                    cursor.execute("""
                        INSERT INTO protein_references (protein_id, pubmed_id, doi, title, authors, journal, publication_date)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                    """, (protein_id, ref.get('pubmed_id'), ref.get('doi'), ref.get('title'),
                          ref.get('authors'), ref.get('journal'), ref.get('publication_date')))

            self.connection.commit()
            return True

        except Exception as e:
            self.connection.rollback()
            print(f"❌ Error saving UniProt metadata: {e}")
            return False


def init_database(db_path="protein_reports.db"):
    """Initialize database and return connection"""
    db = ProteinDatabase(db_path)
    return db


if __name__ == "__main__":
    # Test database initialization
    db = ProteinDatabase()
    print("✓ Database initialized successfully")

    stats = db.get_statistics()
    print("\nDatabase Statistics:")
    for key, value in stats.items():
        print(f"  {key}: {value}")

    db.close()
