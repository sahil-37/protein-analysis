#!/usr/bin/env python3
"""
Configuration module for protein analysis application.
Centralized settings and constants.
"""

import os
from pathlib import Path

# ============================================================================
# API ENDPOINTS
# ============================================================================

UNIPROT_API = {
    'sequence_url': 'https://www.uniprot.org/uniprot/{accession}.fasta',
    'entry_url': 'https://rest.uniprot.org/uniprotkb/{accession}.json',
}

EBI_API = {
    'features_url': 'https://www.ebi.ac.uk/proteins/api/features/{accession}',
}

# ============================================================================
# DATABASE SETTINGS
# ============================================================================

DATABASE = {
    'path': 'protein_reports.db',
    'schema': 'data/database_schema.sql',
}

# ============================================================================
# REPORT SETTINGS
# ============================================================================

REPORT = {
    'template': 'templates/offline_template.html',
    'output_dir': '.',
    'prefix': 'report_',
}

# ============================================================================
# ANALYSIS SETTINGS
# ============================================================================

ANALYSIS = {
    'pH_min': 2.0,
    'pH_max': 12.0,
    'pH_step': 0.1,
    'clean_nonstandard_aa': {
        'U': 'C',  # Selenocysteine -> Cysteine
        'B': 'X',  # Aspartic acid or Asparagine -> Unknown
        'Z': 'X',  # Glutamic acid or Glutamine -> Unknown
    }
}

# ============================================================================
# LOGGING SETTINGS
# ============================================================================

LOGGING = {
    'verbose': True,
    'show_progress': True,
}

# ============================================================================
# RUNTIME DIRECTORIES
# ============================================================================

PROJECT_DIR = Path(__file__).parent
DATA_DIR = PROJECT_DIR / 'data'
TEMP_DIR = PROJECT_DIR / '.temp'

# Create directories if they don't exist
DATA_DIR.mkdir(exist_ok=True)
TEMP_DIR.mkdir(exist_ok=True)

# ============================================================================
# PROTEIN FEATURE SETTINGS
# ============================================================================

PTM_CATEGORIES = [
    'MOD_RES',      # Modified residue
    'CARBOHYD',     # Carbohydrate attachment
    'DISULFID',     # Disulfide bond
    'LIPID',        # Lipid attachment
    'CROSSLNK',     # Cross-link
]

# ============================================================================
# DEFAULT VALUES
# ============================================================================

DEFAULT_ACCESSION = 'P15018'
TIMEOUT = 10  # seconds
RETRIES = 3

# ============================================================================
# FUNCTION: Get Full Path
# ============================================================================

def get_config_path(key_path):
    """
    Get a configuration value by dot-notation path.

    Example:
        get_config_path('DATABASE.path') -> 'protein_reports.db'
    """
    parts = key_path.split('.')
    config = globals()

    for part in parts:
        if part in config:
            config = config[part]
        else:
            raise KeyError(f"Configuration key not found: {key_path}")

    return config
