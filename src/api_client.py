#!/usr/bin/env python3
"""
API Client module for UniProt, EBI, and PDB interactions.
Handles all external API calls.
"""

import requests
from typing import Dict, Any, Optional, Tuple
from src.config import UNIPROT_API, EBI_API, TIMEOUT, RETRIES


class APIClient:
    """Unified API client for biological databases."""

    def __init__(self, timeout: int = TIMEOUT):
        """
        Initialize API client.

        Args:
            timeout: Request timeout in seconds
        """
        self.timeout = timeout
        self.headers = {'Accept': 'application/json'}

    def fetch_protein_sequence(self, accession_id: str) -> Optional[str]:
        """
        Fetch protein sequence from UniProt.

        Args:
            accession_id: UniProt accession ID (e.g., 'P04637')

        Returns:
            Protein sequence string or None if failed
        """
        try:
            url = UNIPROT_API['sequence_url'].format(accession=accession_id)
            response = requests.get(url, timeout=self.timeout)
            response.raise_for_status()

            # Extract sequence from FASTA
            fasta_lines = response.text.split('\n')
            sequence = ''.join(
                line.strip() for line in fasta_lines
                if not line.startswith('>')
            )

            return sequence if sequence else None

        except requests.RequestException as e:
            print(f"❌ Error fetching sequence for {accession_id}: {e}")
            return None

    def fetch_protein_features(self, accession_id: str) -> Optional[Dict[str, Any]]:
        """
        Fetch protein features from EBI API.

        Args:
            accession_id: UniProt accession ID

        Returns:
            Features JSON dictionary or None if failed
        """
        try:
            url = EBI_API['features_url'].format(accession=accession_id)
            response = requests.get(url, headers=self.headers, timeout=self.timeout)
            response.raise_for_status()
            return response.json()

        except requests.RequestException as e:
            print(f"⚠️  Warning: Failed to fetch features for {accession_id}: {e}")
            return None

    def fetch_uniprot_entry(self, accession_id: str) -> Optional[Dict[str, Any]]:
        """
        Fetch complete UniProt entry (REST API).

        Args:
            accession_id: UniProt accession ID

        Returns:
            UniProt entry JSON or None if failed
        """
        try:
            url = UNIPROT_API['entry_url'].format(accession=accession_id)
            response = requests.get(url, headers=self.headers, timeout=self.timeout)
            response.raise_for_status()
            return response.json()

        except requests.RequestException as e:
            print(f"⚠️  Warning: Failed to fetch entry for {accession_id}: {e}")
            return None

    def fetch_pdb_structure(self, accession_id: str) -> Optional[str]:
        """
        Fetch PDB structure ID from UniProt entry.

        Args:
            accession_id: UniProt accession ID

        Returns:
            PDB ID or None if not found
        """
        try:
            data = self.fetch_uniprot_entry(accession_id)
            if not data:
                return None

            # Look for PDB cross-references
            for xref in data.get('uniProtKBCrossReferences', []):
                if xref.get('database') == 'PDB':
                    return xref.get('id')

            return None

        except Exception as e:
            print(f"⚠️  Could not fetch PDB structure: {e}")
            return None

    def fetch_all_protein_data(self, accession_id: str) -> Dict[str, Any]:
        """
        Fetch all protein data (sequence, features, metadata, PDB).

        Args:
            accession_id: UniProt accession ID

        Returns:
            Dictionary containing all fetched data
        """
        data = {
            'accession_id': accession_id,
            'sequence': None,
            'features': None,
            'entry': None,
            'pdb_id': None,
        }

        print(f"📥 Fetching protein data for {accession_id}...")

        data['sequence'] = self.fetch_protein_sequence(accession_id)
        if data['sequence']:
            print(f"   ✅ Fetched sequence ({len(data['sequence'])} aa)")

        # Fetch UniProt entry (primary source for features)
        data['entry'] = self.fetch_uniprot_entry(accession_id)
        if data['entry']:
            print(f"   ✅ Fetched UniProt entry")
            # Extract features from UniProt entry (authoritative source)
            data['features'] = data['entry'].get('features', [])
            if data['features']:
                print(f"   ✅ Extracted {len(data['features'])} features from UniProt")

        # Try EBI as fallback if UniProt features are incomplete
        if not data['features']:
            ebi_features = self.fetch_protein_features(accession_id)
            if ebi_features:
                # Restructure EBI response to match UniProt format
                data['features'] = ebi_features.get('features', [])
                print(f"   ✅ Fetched {len(data['features'])} features from EBI (fallback)")

        # Get PDB ID from UniProt if available
        if data['entry']:
            data['pdb_id'] = self.fetch_pdb_structure(accession_id)
            if data['pdb_id']:
                print(f"   ✅ Found PDB structure: {data['pdb_id']}")

        return data
