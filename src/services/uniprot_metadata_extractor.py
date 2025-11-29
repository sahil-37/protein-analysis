#!/usr/bin/env python3
"""
UniProt Metadata Extractor
Extracts comprehensive UniProt metadata from REST API for database storage
"""

import requests
from typing import Dict, Any, List, Optional


class UniProtMetadataExtractor:
    """Extract comprehensive UniProt metadata from UniProt REST API"""

    BASE_URL = "https://rest.uniprot.org/uniprotkb"

    def __init__(self):
        """Initialize extractor"""
        self.headers = {"Accept": "application/json"}

    def extract_metadata(self, accession_id: str) -> Optional[Dict[str, Any]]:
        """
        Extract comprehensive metadata for a protein from UniProt

        Args:
            accession_id: UniProt accession ID

        Returns:
            Dictionary containing all extracted metadata or None if failed
        """
        try:
            url = f"{self.BASE_URL}/{accession_id}.json"
            response = requests.get(url, headers=self.headers, timeout=10)
            response.raise_for_status()
            data = response.json()

            metadata = {
                'secondary_accessions': self._extract_secondary_accessions(data),
                'genes': self._extract_genes(data),
                'protein_names': self._extract_protein_names(data),
                'comments': self._extract_comments(data),
                'cross_references': self._extract_cross_references(data),
                'keywords': self._extract_keywords(data),
                'go_terms': self._extract_go_terms(data),
                'interactions': self._extract_interactions(data),
                'references': self._extract_references(data),
            }

            return metadata

        except requests.RequestException as e:
            print(f"⚠️  Could not fetch UniProt metadata for {accession_id}: {e}")
            return None

    def _extract_secondary_accessions(self, data: Dict) -> List[str]:
        """Extract secondary accessions"""
        return data.get('secondaryAccessions', [])

    def _extract_genes(self, data: Dict) -> List[Dict[str, str]]:
        """Extract gene names including synonyms"""
        genes = []

        for gene_obj in data.get('genes', []):
            # Primary gene name
            if 'geneName' in gene_obj:
                genes.append({
                    'gene_name': gene_obj['geneName'].get('value', ''),
                    'gene_type': 'primary'
                })

            # Gene synonyms
            for synonym in gene_obj.get('synonyms', []):
                genes.append({
                    'gene_name': synonym.get('value', ''),
                    'gene_type': 'synonym'
                })

            # ORF (Open Reading Frame) names
            for orf in gene_obj.get('orfNames', []):
                genes.append({
                    'gene_name': orf.get('value', ''),
                    'gene_type': 'orf'
                })

        return genes

    def _extract_protein_names(self, data: Dict) -> List[Dict[str, str]]:
        """Extract protein names (recommended, alternative, short names, etc.)"""
        names = []
        protein_desc = data.get('proteinDescription', {})

        # Recommended name
        if 'recommendedName' in protein_desc:
            rec_name = protein_desc['recommendedName']
            if 'fullName' in rec_name:
                names.append({
                    'name_type': 'recommended_full',
                    'name_value': rec_name['fullName'].get('value', '')
                })
            if 'shortNames' in rec_name:
                for short_name in rec_name['shortNames']:
                    names.append({
                        'name_type': 'recommended_short',
                        'name_value': short_name.get('value', '')
                    })

        # Alternative names
        for alt_name in protein_desc.get('alternativeNames', []):
            if 'fullName' in alt_name:
                names.append({
                    'name_type': 'alternative_full',
                    'name_value': alt_name['fullName'].get('value', '')
                })
            if 'shortNames' in alt_name:
                for short_name in alt_name['shortNames']:
                    names.append({
                        'name_type': 'alternative_short',
                        'name_value': short_name.get('value', '')
                    })

        # Submitted name
        if 'submittedName' in protein_desc:
            sub_name = protein_desc['submittedName']
            if 'fullName' in sub_name:
                names.append({
                    'name_type': 'submitted_full',
                    'name_value': sub_name['fullName'].get('value', '')
                })

        return names

    def _extract_comments(self, data: Dict) -> List[Dict[str, str]]:
        """Extract protein comments (function, activity, location, etc.)"""
        comments = []

        for comment in data.get('comments', []):
            comment_type = comment.get('commentType', 'unknown')
            comment_text = ''

            # Extract text depending on comment type
            if comment_type == 'FUNCTION':
                comment_text = comment.get('texts', [{}])[0].get('value', '')
            elif comment_type == 'CATALYTIC_ACTIVITY':
                reaction = comment.get('reaction', {})
                comment_text = reaction.get('name', '')
            elif comment_type == 'SUBCELLULAR_LOCATION':
                locations = []
                for location in comment.get('locations', []):
                    locations.append(location.get('value', ''))
                comment_text = '; '.join(locations)
            elif comment_type == 'TISSUE_SPECIFICITY':
                comment_text = comment.get('texts', [{}])[0].get('value', '')
            elif comment_type == 'DISEASE':
                disease = comment.get('disease', {})
                comment_text = f"{disease.get('name', '')} - {disease.get('description', '')}"
            elif comment_type == 'SIMILARITY':
                comment_text = comment.get('texts', [{}])[0].get('value', '')
            else:
                comment_text = comment.get('text', '') or str(comment.get('texts', [{}])[0].get('value', ''))

            if comment_text:
                comments.append({
                    'comment_type': comment_type,
                    'comment_text': comment_text
                })

        return comments

    def _extract_cross_references(self, data: Dict) -> List[Dict[str, Any]]:
        """Extract cross-references to other databases"""
        xrefs = []

        for xref in data.get('uniProtKBCrossReferences', []):
            database = xref.get('database', '')
            xref_id = xref.get('id', '')

            # Extract properties for certain databases
            properties = {}
            if 'properties' in xref:
                for prop in xref['properties']:
                    key = prop.get('key', '')
                    value = prop.get('value', '')
                    if key:
                        properties[key] = value

            xrefs.append({
                'database_name': database,
                'database_id': xref_id,
                'properties': properties
            })

        return xrefs

    def _extract_keywords(self, data: Dict) -> List[Dict[str, str]]:
        """Extract keywords"""
        keywords = []

        for keyword in data.get('keywords', []):
            keywords.append({
                'keyword_id': keyword.get('id', ''),
                'keyword_name': keyword.get('name', '')
            })

        return keywords

    def _extract_go_terms(self, data: Dict) -> List[Dict[str, str]]:
        """Extract Gene Ontology (GO) terms"""
        go_terms = []

        for xref in data.get('uniProtKBCrossReferences', []):
            if xref.get('database') == 'GO':
                go_id = xref.get('id', '')
                properties = {p['key']: p['value'] for p in xref.get('properties', [])}

                go_terms.append({
                    'go_id': go_id,
                    'go_term': properties.get('term', ''),
                    'go_aspect': properties.get('aspect', ''),
                    'evidence_code': properties.get('evidence', '')
                })

        return go_terms

    def _extract_interactions(self, data: Dict) -> List[Dict[str, Any]]:
        """Extract protein-protein interactions"""
        interactions = []

        for comment in data.get('comments', []):
            if comment.get('commentType') == 'INTERACTION':
                for interaction in comment.get('interactions', []):
                    interactant = interaction.get('interactantOne', {})
                    interaction_type = interaction.get('interactionType', '')

                    # Try to get accession ID
                    interactant_acc = ''
                    if 'uniProtKBAccession' in interactant:
                        interactant_acc = interactant['uniProtKBAccession']
                    elif 'chainId' in interactant:
                        interactant_acc = interactant['chainId']

                    interactions.append({
                        'interactant_accession': interactant_acc,
                        'interactant_name': interactant.get('name', ''),
                        'interaction_type': interaction_type,
                        'experiments': len(interaction.get('experiments', []))
                    })

        return interactions

    def _extract_references(self, data: Dict) -> List[Dict[str, str]]:
        """Extract publication references"""
        references = []

        for ref in data.get('references', []):
            citation = ref.get('citation', {})
            citation_type = citation.get('citationType', '')

            # Extract common fields
            title = citation.get('title', '')

            # Extract authors - can be authoringGroup (list) or authors (list of dicts/strings)
            authors_list = ''
            authoring_group = citation.get('authoringGroup', [])
            if isinstance(authoring_group, list) and authoring_group:
                if isinstance(authoring_group[0], dict):
                    authors_list = authoring_group[0].get('name', '')
            elif isinstance(authoring_group, dict):
                authors_list = authoring_group.get('name', '')

            if not authors_list:
                authors = citation.get('authors', [])
                # Authors can be a list of strings or a list of dicts
                author_names = []
                for a in authors:
                    if isinstance(a, dict):
                        author_names.append(a.get('name', ''))
                    else:
                        author_names.append(str(a))
                authors_list = '; '.join(author_names)

            journal = citation.get('journal', '')
            publication_date = citation.get('publicationDate', '')

            # Extract identifiers
            pubmed_id = ''
            doi = ''

            for db_ref in citation.get('dbReferences', []):
                if db_ref.get('type') == 'PubMed':
                    pubmed_id = db_ref.get('id', '')
                elif db_ref.get('type') == 'DOI':
                    doi = db_ref.get('id', '')

            if title or pubmed_id or doi:
                references.append({
                    'pubmed_id': pubmed_id,
                    'doi': doi,
                    'title': title,
                    'authors': authors_list,
                    'journal': journal,
                    'publication_date': publication_date
                })

        return references


def extract_uniprot_metadata(accession_id: str) -> Optional[Dict[str, Any]]:
    """
    Convenience function to extract UniProt metadata

    Args:
        accession_id: UniProt accession ID

    Returns:
        Dictionary containing metadata or None if failed
    """
    extractor = UniProtMetadataExtractor()
    return extractor.extract_metadata(accession_id)


if __name__ == "__main__":
    import json

    # Test extraction
    test_accession = "P04637"
    print(f"Extracting metadata for {test_accession}...")

    metadata = extract_uniprot_metadata(test_accession)
    if metadata:
        print(json.dumps(metadata, indent=2))
    else:
        print("Failed to extract metadata")
