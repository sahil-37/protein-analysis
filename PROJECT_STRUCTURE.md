# UniProt Lookup - Project Structure Guide

## Overview

This is a modular Python application for comprehensive protein analysis using UniProt data. The project is organized for clarity, maintainability, and scalability, with a clear separation of concerns across different modules.

**Total Size**: ~186 KB
**Core Files**: 14 essential modules
**Architecture**: Modular with 5 core modules + supporting services

---

## Directory Structure

```
uniprot_lookup/
├── src/                              # Main application code
│   ├── __init__.py                   # Package initialization
│   ├── config.py                     # Configuration & constants (ROOT MODULE)
│   ├── api_client.py                 # External API interactions
│   ├── protein_analyzer.py           # Biophysical property calculations
│   ├── data_manager.py               # Data persistence interface
│   ├── main.py                       # Application entry point & orchestration
│   └── services/                     # Supporting service modules
│       ├── __init__.py               # Package initialization
│       ├── protein_database.py       # SQLite database operations
│       ├── protein_enrichment_service.py  # Protein enrichment pipeline
│       └── uniprot_metadata_extractor.py  # UniProt metadata extraction
│
├── data/                             # Data files & schema
│   └── database_schema.sql           # SQLite schema (18 tables)
│
├── templates/                        # HTML report templates
│   └── offline_template.html         # Interactive protein report (69 KB)
│
├── docs/                             # Documentation
│   ├── README.md                     # User guide
│   ├── ARCHITECTURE.md               # Detailed architecture docs
│   └── SUMMARY.md                    # Project summary
│
├── tests/                            # Test files (placeholder)
│
├── requirements.txt                  # Python dependencies
├── .env.example                      # Environment template
├── .gitignore                        # Git ignore rules
└── PROJECT_STRUCTURE.md              # This file

```

---

## Module Descriptions

### Core Modules (src/)

#### 1. **config.py** (3.4 KB) - Configuration Root
- **Purpose**: Centralized configuration and constants
- **Dependencies**: None (depends on no other modules)
- **Key Variables**:
  - `UNIPROT_API`: UniProt REST API endpoints
  - `EBI_API`: EBI Proteins API endpoints
  - `DATABASE`: Database configuration
  - `REPORT`: Report generation settings
  - `ANALYSIS`: Analysis parameters (pH range, cleaning rules)
  - `PTM_CATEGORIES`: Post-translational modification types
  - `PROJECT_DIR`, `DATA_DIR`, `TEMP_DIR`: Directory paths
- **Usage**: Imported by all other modules for configuration values

**Key Configuration:**
```python
DATABASE = {
    'path': 'protein_reports.db',
    'schema': 'data/database_schema.sql',  # Moved to data/ directory
}

REPORT = {
    'template': 'templates/offline_template.html',  # Updated path
    'output_dir': '.',
    'prefix': 'report_',
}
```

---

#### 2. **api_client.py** (4.7 KB) - API Client
- **Purpose**: Unified interface for all external API calls
- **Dependencies**: config.py
- **Class**: `APIClient`
- **Methods**:
  - `fetch_protein_sequence(accession_id)`: Get FASTA sequence from UniProt
  - `fetch_protein_features(accession_id)`: Get protein features from EBI
  - `fetch_uniprot_entry(accession_id)`: Get complete UniProt JSON entry
  - `fetch_pdb_structure(accession_id)`: Get PDB ID from UniProt
  - `fetch_all_protein_data(accession_id)`: Fetch all data in one call
- **Imports**: `requests`, `config.py`
- **Usage**: Called by main.py and protein_enrichment_service.py

**Example Usage:**
```python
from src.api_client import APIClient

client = APIClient()
sequence = client.fetch_protein_sequence('P04637')
entry = client.fetch_uniprot_entry('P04637')
all_data = client.fetch_all_protein_data('P04637')
```

---

#### 3. **protein_analyzer.py** (3.8 KB) - Biophysical Analysis
- **Purpose**: Calculate biophysical properties from protein sequences
- **Dependencies**: config.py, BioPython
- **Class**: `ProteinAnalyzer`
- **Methods**:
  - `get_molecular_weight()`: Molecular weight in Da/kDa
  - `get_isoelectric_point()`: Isoelectric point (pI)
  - `get_gravy()`: Grand Average of Hydropathy
  - `get_aromaticity()`: Aromatic amino acid percentage
  - `get_instability_index()`: Protein stability indicator
  - `get_cysteine_count()`: Count of cysteine residues
  - `get_extinction_coefficients()`: UV absorption
  - `get_secondary_structure_fractions()`: Helix/Sheet/Coil predictions
  - `get_all_properties()`: Complete analysis in one call
  - `is_stable()`: Boolean stability prediction
- **Imports**: `Bio.SeqUtils.ProtParam`, `config.py`
- **Usage**: Called by main.py for protein analysis

**Example Usage:**
```python
from src.protein_analyzer import ProteinAnalyzer

analyzer = ProteinAnalyzer('MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYTNGK...')
mw = analyzer.get_molecular_weight()
pi = analyzer.get_isoelectric_point()
properties = analyzer.get_all_properties()
```

---

#### 4. **data_manager.py** (3.3 KB) - Data Persistence
- **Purpose**: High-level interface for saving and retrieving protein data
- **Dependencies**: config.py, services/protein_database.py
- **Class**: `DataManager`
- **Methods**:
  - `save_protein(protein_data)`: Store protein analysis results
  - `save_report(report_path, metadata)`: Store generated report
  - `get_protein(accession_id)`: Retrieve protein data
  - `list_proteins()`: List all stored proteins
  - `search_proteins(query)`: Search database
  - `get_statistics()`: Get database statistics
- **Imports**: `protein_database.py`, `uniprot_metadata_extractor.py`
- **Usage**: Called by main.py for data persistence

**Example Usage:**
```python
from src.data_manager import DataManager

dm = DataManager()
dm.save_protein(protein_data)
protein = dm.get_protein('P04637')
stats = dm.get_statistics()
```

---

#### 5. **main.py** (9.0 KB) - Application Orchestrator
- **Purpose**: Main entry point and workflow orchestration
- **Dependencies**: All other core modules
- **Class**: `ProteinReportGenerator`
- **Key Methods**:
  - `generate_report(accession_id)`: Complete workflow orchestration
  - `generate_batch_reports(accession_ids)`: Process multiple proteins
- **Imports**: All other modules (api_client, protein_analyzer, data_manager, etc.)
- **Usage**: Run from command line or import for programmatic use

**Example Usage:**
```bash
python src/main.py P04637
```

**Programmatic Usage:**
```python
from src.main import ProteinReportGenerator

generator = ProteinReportGenerator()
report_path = generator.generate_report('P04637')
```

---

### Service Modules (src/services/)

#### 6. **protein_database.py** (18 KB) - Database Operations
- **Purpose**: SQLite database management and CRUD operations
- **Dependencies**: config.py (for schema path)
- **Class**: `ProteinDatabase`
- **Features**:
  - Automatic schema initialization from `data/database_schema.sql`
  - 18 normalized tables (9 calculated, 9 metadata)
  - Connection pooling and row factory setup
  - Comprehensive CRUD methods
- **Key Path**: Schema file at `data/database_schema.sql` (updated from parent directory)
- **Schema Tables**:
  - **Calculated Data**: molecular_weight, isoelectric_point, charge, stability, etc.
  - **UniProt Metadata**: genes, proteins, organisms, comments, interactions, etc.
- **Usage**: Used by data_manager.py

**Database Connection Example:**
```python
from src.services.protein_database import ProteinDatabase

db = ProteinDatabase('protein_reports.db')
db.initialize_db()  # Creates schema if needed
db.close()
```

---

#### 7. **protein_enrichment_service.py** (32 KB) - Enrichment Pipeline
- **Purpose**: Complete protein enrichment workflow
- **Dependencies**: api_client.py, BioPython
- **Class**: `ProteinEnrichmentService`
- **Features**:
  - Biophysical property calculation (BiophysicalCalculator)
  - Structural feature prediction (StructuralPredictor)
  - PTM and feature extraction from UniProt
  - Charge profile generation
- **Key Data Classes**:
  - `ProteinEnrichmentProfile`: Complete enrichment output
  - `BiophysicalProperties`: Calculated properties
  - `ComputedCharges`: pH-dependent charge data
  - `PredictedFeatures`: Structural predictions
- **Public Functions**:
  - `enrich_protein(accession_id)`: Single protein enrichment
  - `enrich_protein_batch(accession_ids)`: Batch enrichment
- **Usage**: Standalone enrichment service or via main.py

**Example Usage:**
```python
from src.services.protein_enrichment_service import enrich_protein

profile = enrich_protein('P04637')
profile.to_json('p53_enriched.json')
```

---

#### 8. **uniprot_metadata_extractor.py** (12 KB) - Metadata Extraction
- **Purpose**: Extract comprehensive UniProt metadata from REST API
- **Dependencies**: requests library
- **Class**: `UniProtMetadataExtractor`
- **Extraction Methods**:
  - Secondary accessions
  - Gene names and synonyms
  - Protein names
  - Organism information
  - Comments and cross-references
  - Keywords and GO terms
  - Protein interactions
  - References and authors
- **Usage**: Used by protein_enrichment_service.py and data_manager.py

**Example Usage:**
```python
from src.services.uniprot_metadata_extractor import UniProtMetadataExtractor

extractor = UniProtMetadataExtractor()
metadata = extractor.extract_metadata('P04637')
```

---

## Data & Templates

### data/ Directory
- **database_schema.sql** (7.4 KB): SQLite schema with 18 tables
  - Stores calculated biophysical properties
  - Stores extracted UniProt metadata
  - Indexed for fast queries

### templates/ Directory
- **offline_template.html** (69 KB): Interactive protein report
  - Nightingale protein viewer integration
  - PTM visualization with vertical alignment
  - Biophysical property display
  - Responsive design

---

## Dependency Graph

```
config.py (ROOT - no dependencies)
    ↓
    ├─ api_client.py
    │   ├─ main.py
    │   └─ protein_enrichment_service.py
    │
    ├─ protein_analyzer.py
    │   └─ main.py
    │
    ├─ data_manager.py
    │   ├─ main.py
    │   ├─ protein_database.py
    │   └─ uniprot_metadata_extractor.py
    │
    └─ main.py
        └─ (orchestrates all modules)
```

---

## How Files Are Used

### Typical Workflow

1. **User runs main.py** with protein accession ID
   ```bash
   python src/main.py P04637
   ```

2. **main.py (orchestrator)** delegates to:
   - `api_client.py`: Fetch sequence and metadata from UniProt
   - `protein_analyzer.py`: Calculate biophysical properties
   - `protein_enrichment_service.py`: Enrich with predictions
   - `data_manager.py`: Save results to database
   - Generates HTML report using template

3. **api_client.py** makes HTTP requests to:
   - UniProt REST API (sequence, JSON entry)
   - EBI Proteins API (features)

4. **protein_database.py** handles:
   - Schema initialization from `data/database_schema.sql`
   - Storing calculated properties
   - Storing extracted metadata

5. **protein_enrichment_service.py** provides:
   - BiophysicalCalculator for property calculations
   - StructuralPredictor for feature predictions
   - PTM extraction from UniProt features

6. **Report generation** uses:
   - `templates/offline_template.html` for HTML structure
   - Nightingale viewer for protein visualization
   - Data from all calculations above

---

## Configuration & Paths

### Path Resolution

All paths are configured in `config.py` and resolve from project root:

```python
PROJECT_DIR = Path(__file__).parent  # src/
# Then go up to project root:
# src/../../ (parent directory of src/) = project root

DATABASE['schema'] = 'data/database_schema.sql'
REPORT['template'] = 'templates/offline_template.html'
```

### Database Path

In `protein_database.py:28`:
```python
schema_file = Path(__file__).parent.parent.parent / "data" / "database_schema.sql"
# Resolves to: project_root/data/database_schema.sql
```

---

## Python Dependencies

See `requirements.txt`:
- `requests`: HTTP API calls
- `biopython`: Biophysical calculations
- `jinja2`: Template rendering (if used)
- Python 3.8+

---

## Import Patterns

### From Modules
```python
# In src/main.py
from src.config import REPORT, ANALYSIS
from src.api_client import APIClient
from src.protein_analyzer import ProteinAnalyzer
from src.data_manager import DataManager
```

### From Services
```python
# In src/data_manager.py
from src.services.protein_database import ProteinDatabase
from src.services.uniprot_metadata_extractor import UniProtMetadataExtractor
```

---

## Key Design Principles

1. **Separation of Concerns**: Each module has a single responsibility
2. **No Circular Dependencies**: All dependencies form a directed acyclic graph
3. **Configuration Root**: All settings centralized in config.py
4. **Single Entry Point**: main.py orchestrates the workflow
5. **Service Pattern**: Supporting functionality in services/ directory
6. **Data Isolation**: Database operations isolated in database module
7. **API Abstraction**: All external calls go through api_client.py

---

## Development Notes

### Adding New Features

1. **New calculation**: Add to `protein_analyzer.py`
2. **New API source**: Add to `api_client.py`
3. **New database table**: Update `data/database_schema.sql` and `protein_database.py`
4. **New report feature**: Update `templates/offline_template.html`

### Testing

Place test files in `tests/` directory following naming convention `test_*.py`

### Running

```bash
# Single protein
python src/main.py P04637

# With options
python src/main.py P04637 --output report.html

# Batch processing
python src/main.py P04637 P68871 P15018
```

---

## File Sizes Summary

| File | Size | Purpose |
|------|------|---------|
| offline_template.html | 69 KB | HTML report template |
| protein_enrichment_service.py | 32 KB | Enrichment pipeline |
| protein_database.py | 18 KB | Database operations |
| uniprot_metadata_extractor.py | 12 KB | Metadata extraction |
| main.py | 9.0 KB | Orchestrator |
| api_client.py | 4.7 KB | API client |
| database_schema.sql | 7.4 KB | Database schema |
| config.py | 3.4 KB | Configuration |
| data_manager.py | 3.3 KB | Data persistence |
| protein_analyzer.py | 3.8 KB | Biophysical analysis |
| README.md | 11 KB | User guide |
| ARCHITECTURE.md | 9.3 KB | Architecture docs |
| SUMMARY.md | 2.5 KB | Project summary |
| **Total** | **~186 KB** | |

---

## Quick Reference

### Import Everything
```python
from src.config import *
from src.api_client import APIClient
from src.protein_analyzer import ProteinAnalyzer
from src.data_manager import DataManager
from src.main import ProteinReportGenerator
```

### Generate Report Programmatically
```python
from src.main import ProteinReportGenerator

gen = ProteinReportGenerator()
report_path = gen.generate_report('P04637')
print(f"Report saved to: {report_path}")
```

### Access Database
```python
from src.data_manager import DataManager

dm = DataManager()
protein = dm.get_protein('P04637')
all_proteins = dm.list_proteins()
```

---

**Last Updated**: November 29, 2024
**Version**: 2.0 (Reorganized with modular structure)
