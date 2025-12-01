# UniProt Lookup - Comprehensive Protein Analysis & Reporting Tool

A modular Python application for comprehensive protein analysis using UniProt data, biophysical property calculations, and interactive HTML report generation.

**Latest Version**: 2.0 (Reorganized with modular architecture)
**Status**: Production Ready ✓
**Project Size**: 186 KB (down from 1.6 MB)

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Features](#features)
3. [Project Structure](#project-structure)
4. [Installation](#installation)
5. [Usage Guide](#usage-guide)
6. [API Reference](#api-reference)
7. [Database](#database)
8. [Report Generation](#report-generation)
9. [Examples](#examples)
10. [Troubleshooting](#troubleshooting)

---

## Quick Start

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Generate Your First Protein Report

```bash
# Generate report for a single protein (Hemoglobin)
python -m src.main P69905

# Or use a different accession ID (p53)
python -m src.main P04637
```

### 3. Check Generated Files

The report is saved as `report_P69905.html` and protein data is stored in `protein_reports.db`

### Done! 🎉

That's it! The protein analysis report has been generated and data stored in the database.

---

## Features

### Core Features

✓ **UniProt Data Fetching**
- Fetch complete protein information from UniProt REST API
- Extract sequence, metadata, features, and annotations
- Extract subcellular location information
- Support for batch processing multiple proteins

✓ **Biophysical Property Calculations**
- Molecular weight (Da & kDa)
- Isoelectric point (pI) using Henderson-Hasselbalch
- GRAVY score (hydropathy)
- Aromaticity
- Instability index
- Secondary structure (theoretical from DSSP/STRIDE, experimental from UniProt)
- Extinction coefficients (reduced & oxidized forms)
- pH-dependent charge profile (pH 2-12)

✓ **Post-Translational Modification (PTM) Analysis**
- Extract all PTM types from UniProt
- Visualize PTM locations on protein sequence
- Support for:
  - Modified residues (MOD_RES)
  - Carbohydrate attachments (CARBOHYD)
  - Disulfide bonds (DISULFID)
  - Lipid attachments (LIPID)
  - Cross-links (CROSSLNK)

✓ **Protein Enrichment**
- Signal peptide prediction
- Transmembrane domain prediction
- Disulfide bond prediction
- PTM/processing feature extraction
- Charge profile generation

✓ **Database Storage**
- SQLite database with 18 normalized tables
- Store all calculated biophysical properties
- Store extracted UniProt metadata (location, comments, GO terms, references, etc.)
- Full synchronization: all report data is backed by database
- Indexed for fast queries
- Support for batch storage and retrieval

✓ **Interactive Report Generation**
- Professional HTML reports with Nightingale viewer
- Interactive protein sequence visualization
- PTM visualization with vertical alignment
- Biophysical properties display
- Charge profile graphs
- Responsive design (desktop & mobile)

---

## Project Structure

```
uniprot_lookup/
├── src/                              # Main application code
│   ├── __init__.py
│   ├── config.py                     # Configuration & constants (ROOT)
│   ├── api_client.py                 # External API interactions
│   ├── protein_analyzer.py           # Biophysical calculations
│   ├── data_manager.py               # Data persistence interface
│   ├── main.py                       # Entry point & orchestration
│   └── services/
│       ├── __init__.py
│       ├── protein_database.py       # SQLite database operations
│       ├── protein_enrichment_service.py  # Protein enrichment pipeline
│       └── uniprot_metadata_extractor.py  # UniProt metadata extraction
│
├── data/                             # Data files & schema
│   └── database_schema.sql           # SQLite schema (18 tables)
│
├── templates/                        # HTML report templates
│   └── offline_template.html         # Interactive protein report
│
├── docs/                             # Documentation
│   ├── README.md                     # Original user guide
│   ├── ARCHITECTURE.md               # Architecture documentation
│   └── SUMMARY.md                    # Project summary
│
├── tests/                            # Test files (placeholder)
│
├── PROJECT_STRUCTURE.md              # Detailed structure guide
├── REORGANIZATION_SUMMARY.md         # Changes from v1.0 to v2.0
├── requirements.txt                  # Python dependencies
├── .env.example                      # Environment template
├── .gitignore                        # Git ignore rules
└── README_COMPREHENSIVE.md           # This file

```

---

## Installation

### Prerequisites

- Python 3.8+
- pip (Python package manager)

### Step 1: Clone/Download Project

```bash
cd uniprot_lookup
```

### Step 2: Install Dependencies

```bash
pip install -r requirements.txt
```

**Required packages:**
- `requests` - HTTP API calls
- `biopython` - Biophysical calculations
- `jinja2` - Template rendering (optional)

### Step 3: Verify Installation

```bash
python3 -c "from src.config import *; from src.api_client import APIClient; print('✓ All imports working!')"
```

Expected output:
```
✓ All imports working!
```

---

## Usage Guide

### Command Line Interface

#### Generate a Single Protein Report

```bash
python -m src.main P69905
```

**Output:**
```
============================================================
  PROTEIN ANALYSIS REPORT GENERATOR
============================================================

🔬 Processing: P69905
📥 Fetching protein data for P69905...
   ✅ Fetched sequence (142 aa)
   ✅ Fetched features
   ✅ Fetched UniProt entry
   ✅ Found PDB structure: 1A00
   🧮 Calculating biophysical properties...
   ✅ Calculated properties
   🔋 Generating charge profile...
   💾 Saving data to database...
   ✅ Data saved to database (ID: 1)
   📝 Rendering template...
   ✅ Report saved: report_P69905.html

============================================================
```

#### Generate Multiple Protein Reports (Batch)

```bash
python -m src.main P69905 P04637 P12345
```

This will generate reports for three proteins sequentially.

### Python API

#### Basic Usage

```python
from src.main import ProteinReportGenerator

# Initialize generator
generator = ProteinReportGenerator()

# Generate report for single protein
report_path = generator.generate_report('P69905')
print(f"Report saved to: {report_path}")

# Generate reports for multiple proteins
proteins = ['P69905', 'P04637', 'P12345']
generator.generate_batch_reports(proteins)
```

#### Using Individual Modules

```python
# Fetch data from APIs
from src.api_client import APIClient

client = APIClient()
sequence = client.fetch_protein_sequence('P69905')
entry = client.fetch_uniprot_entry('P69905')
features = client.fetch_protein_features('P69905')

# Calculate biophysical properties
from src.protein_analyzer import ProteinAnalyzer

analyzer = ProteinAnalyzer(sequence)
mw = analyzer.get_molecular_weight()
pi = analyzer.get_isoelectric_point()
properties = analyzer.get_all_properties()

# Store in database
from src.data_manager import DataManager

dm = DataManager()
dm.save_protein({
    'accession_id': 'P69905',
    'sequence': sequence,
    'properties': properties
})

# Retrieve from database
protein = dm.get_protein('P69905')
```

#### Using Enrichment Service

```python
from src.services.protein_enrichment_service import enrich_protein

# Get complete enrichment profile
profile = enrich_protein('P69905')

# Export to JSON
profile.to_json('enriched_protein.json')

# Access profile data
print(f"Protein: {profile.general.protein_name}")
print(f"Molecular Weight: {profile.biophysical.molecular_weight_da} Da")
print(f"pI: {profile.biophysical.isoelectric_point_pi}")
print(f"PTMs: {len(profile.predicted_features.ptm_predictions)}")
```

---

## API Reference

### src.config (Configuration)

**Root module** - All other modules depend on this. No external dependencies.

```python
from src.config import *

# API Endpoints
UNIPROT_API['sequence_url']  # UniProt FASTA endpoint
UNIPROT_API['entry_url']     # UniProt JSON endpoint
EBI_API['features_url']      # EBI Features endpoint

# Database
DATABASE['path']             # Database file path
DATABASE['schema']           # Schema file path

# Report
REPORT['template']           # HTML template path
REPORT['output_dir']         # Report output directory

# Analysis Parameters
ANALYSIS['pH_min']           # Minimum pH for charge calculation
ANALYSIS['pH_max']           # Maximum pH
ANALYSIS['pH_step']          # pH increment

# Directories
PROJECT_DIR                  # Project root
DATA_DIR                     # data/ directory
TEMP_DIR                     # Temporary files directory
```

### src.api_client (API Client)

Unified interface for all external API calls.

```python
from src.api_client import APIClient

client = APIClient(timeout=10)

# Fetch sequence
sequence = client.fetch_protein_sequence('P69905')
# Returns: 'MSVLSPADKTNVKAAWGKVG...'

# Fetch features from EBI
features = client.fetch_protein_features('P69905')
# Returns: {'features': [...], ...}

# Fetch UniProt JSON entry
entry = client.fetch_uniprot_entry('P69905')
# Returns: Full UniProt JSON object

# Fetch PDB structure ID
pdb_id = client.fetch_pdb_structure('P69905')
# Returns: '1A00'

# Fetch all data at once
all_data = client.fetch_all_protein_data('P69905')
# Returns: {'sequence': '...', 'features': {...}, 'entry': {...}, 'pdb_id': '...'}
```

### src.protein_analyzer (Biophysical Analysis)

Calculates biophysical properties from protein sequences.

```python
from src.protein_analyzer import ProteinAnalyzer

analyzer = ProteinAnalyzer('MSVLSPADKTNVKAAWGKVG...')

# Individual property calculations
mw = analyzer.get_molecular_weight()           # Returns: 16500.5 (Da)
pi = analyzer.get_isoelectric_point()          # Returns: 8.72
gravy = analyzer.get_gravy()                   # Returns: 0.0479
aroma = analyzer.get_aromaticity()             # Returns: 0.0704
inst = analyzer.get_instability_index()        # Returns: 6.97
cys = analyzer.get_cysteine_count()            # Returns: 1
ext_coeff = analyzer.get_extinction_coefficients()  # Returns: (9970, 9970)
ss = analyzer.get_secondary_structure_fractions()  # Returns: (helix%, sheet%, coil%)

# Get all properties at once
props = analyzer.get_all_properties()
# Returns: {
#     'molecular_weight': 16500.5,
#     'isoelectric_point': 8.72,
#     'gravy': 0.0479,
#     'aromaticity': 0.0704,
#     'instability_index': 6.97,
#     'cysteine_count': 1,
#     'extinction_coeff_reduced': 9970,
#     'extinction_coeff_oxidized': 9970,
#     'helix_fraction': 40.14,
#     'sheet_fraction': 35.92,
#     'coil_fraction': 26.06
# }

# Check stability
is_stable = analyzer.is_stable()  # Returns: True/False
```

### src.data_manager (Data Persistence)

High-level interface for database operations.

```python
from src.data_manager import DataManager

dm = DataManager()

# Save protein
dm.save_protein(protein_data)

# Get protein by accession ID
protein = dm.get_protein('P69905')

# List all stored proteins
all_proteins = dm.list_proteins()

# Search proteins
results = dm.search_proteins('hemoglobin')

# Get statistics
stats = dm.get_statistics()
```

### src.services.protein_database (Database Operations)

SQLite database management.

```python
from src.services.protein_database import ProteinDatabase

db = ProteinDatabase('protein_reports.db')

# Initialize (creates schema if needed)
db.initialize_db()

# CRUD operations
db.save_protein(protein_data)
protein = db.get_protein('P69905')
all_proteins = db.list_proteins()
db.update_protein('P69905', new_data)
db.delete_protein('P69905')

# Queries
exists = db.protein_exists('P69905')
count = db.count_proteins()

# Close connection
db.close()
```

### src.services.protein_enrichment_service (Enrichment)

Complete protein enrichment pipeline.

```python
from src.services.protein_enrichment_service import enrich_protein, enrich_protein_batch

# Single protein
profile = enrich_protein('P69905')
# Returns: ProteinEnrichmentProfile object

# Multiple proteins
profiles = enrich_protein_batch(['P69905', 'P04637'])
# Returns: List of ProteinEnrichmentProfile objects

# Access profile data
print(profile.general.protein_name)
print(profile.sequence.fasta_sequence)
print(profile.biophysical.molecular_weight_da)
print(profile.charges.charge_at_ph_7_4)
print(profile.predicted_features.ptm_predictions)

# Export to JSON
json_string = profile.to_json()
profile.to_json('output.json')  # Save to file
```

### src.services.uniprot_metadata_extractor (Metadata Extraction)

Extracts comprehensive UniProt metadata.

```python
from src.services.uniprot_metadata_extractor import UniProtMetadataExtractor

extractor = UniProtMetadataExtractor()

# Extract all metadata
metadata = extractor.extract_metadata('P69905')
# Returns: {
#     'secondary_accessions': [...],
#     'genes': [...],
#     'protein_names': [...],
#     'comments': [...],
#     'cross_references': [...],
#     'keywords': [...],
#     'go_terms': [...],
#     'interactions': [...],
#     'references': [...]
# }
```

---

## Database

### Schema Overview

The SQLite database contains **18 tables** organized into two groups. **All data displayed in reports is synchronized with the database**, ensuring consistency and enabling persistent storage and querying.

#### Calculated Data Tables (9)
1. **proteins** - Main protein information with calculated properties
   - Accession ID, protein name, gene name, organism
   - **Subcellular location** (extracted from UniProt)
   - Sequence and sequence length
   - Biophysical properties: molecular weight, isoelectric point, GRAVY score, aromaticity, instability index, cysteine count
   - Extinction coefficients (reduced and oxidized forms)
   - PDB ID and structure availability flag
   - Created/updated timestamps
2. **secondary_structure** - Predicted and experimental secondary structure percentages
   - Helix, sheet, and turn percentages (theoretical from DSSP/STRIDE)
   - Helix, sheet, and turn percentages (experimental from UniProt annotations)
3. **charge_profile** - pH-dependent charge values (pH 2.0 to 12.0)
4. **ptm_summary** - Total PTM count per protein
5. **ptm_details** - Individual PTM site information (position, type, description)
6. **ptm_type_counts** - Count of each PTM type
7. **disulfide_bonds** - Disulfide bond information (if present)
8. **protein_features** - General sequence features (domains, regions, sites, variants, peptides, etc.)
9. **generated_reports** - Report metadata and file paths

#### UniProt Metadata Tables (9)
10. **protein_secondary_accessions** - Alternative accession IDs
11. **protein_genes** - Gene names and synonyms
12. **protein_names** - Protein name variants
13. **protein_comments** - Functional annotations (FUNCTION, TISSUE SPECIFICITY, DISEASE, CAUTION, etc.)
14. **protein_cross_references** - References to external databases (PDB, InterPro, SMART, etc.)
15. **protein_keywords** - Keywords/topics
16. **protein_go_terms** - Gene Ontology annotations
17. **protein_interactions** - Protein-protein interactions
18. **protein_references** - Publication information

### Data Synchronization

**Every data point displayed in the HTML report is backed by the database.** This includes:
- **40+ synchronized data points** across all major sections
- **Biophysical properties** (8 metrics)
- **Secondary structure** (6 values: 3 theoretical + 3 experimental)
- **Charge profile** (101+ pH-charge data points)
- **PTM information** (counts and detailed positions)
- **Protein features** (325+ entries including domains, sites, variants, peptides)
- **Metadata** (name, organism, location, sequence)

### Example Queries

```sql
-- Get stored proteins with core properties
SELECT accession_id, protein_name, organism, subcellular_location,
       molecular_weight, isoelectric_point
FROM proteins;

-- Find proteins by subcellular location
SELECT accession_id, protein_name
FROM proteins
WHERE subcellular_location LIKE '%Nucleus%';

-- Find proteins by molecular weight
SELECT * FROM proteins WHERE molecular_weight > 50000;

-- Get secondary structure for a protein (theoretical vs experimental)
SELECT p.accession_id, p.protein_name,
       ss.helix_predicted, ss.sheet_predicted, ss.turn_predicted,
       ss.helix_experimental, ss.sheet_experimental, ss.turn_experimental
FROM proteins p
JOIN secondary_structure ss ON p.protein_id = ss.protein_id;

-- Get charge at specific pH (physiological pH 7.4)
SELECT p.accession_id, p.protein_name, cp.charge
FROM proteins p
JOIN charge_profile cp ON p.protein_id = cp.protein_id
WHERE cp.ph_value = 7.4;

-- List all PTMs for a protein
SELECT position, ptm_type, description FROM ptm_details WHERE ptm_id = 1;

-- Count PTM types for all proteins
SELECT p.accession_id, ptc.ptm_type, ptc.count
FROM proteins p
JOIN ptm_summary ps ON p.protein_id = ps.protein_id
JOIN ptm_type_counts ptc ON ps.ptm_id = ptc.ptm_id
ORDER BY p.accession_id, ptc.count DESC;

-- Get proteins with specific GO terms
SELECT DISTINCT p.accession_id, p.protein_name
FROM proteins p
JOIN protein_go_terms g ON p.protein_id = g.protein_id
WHERE g.go_term LIKE '%kinase%';

-- Find proteins with disulfide bonds
SELECT p.accession_id, p.protein_name, COUNT(*) as disulfide_count
FROM proteins p
JOIN disulfide_bonds d ON p.protein_id = d.protein_id
GROUP BY p.protein_id;

-- Get all features for a protein
SELECT feature_type, category, description, start_position, end_position
FROM protein_features
WHERE protein_id = 1
ORDER BY start_position;
```

---

## Report Generation

### Understanding the Report

The generated HTML report contains:

1. **Protein Information**
   - Accession ID, Name, Gene Name
   - Organism, Sequence Length
   - **Subcellular Location** (extracted from UniProt)
   - PDB Structure ID (if available)
   - Interactive amino acid sequence viewer

2. **Biophysical Properties**
   - Molecular Weight (Da & kDa)
   - Isoelectric Point (pI)
   - GRAVY Score
   - Aromaticity
   - Instability Index
   - Cysteine Count
   - Extinction Coefficients (Reduced & Oxidized forms in M⁻¹cm⁻¹)

3. **Secondary Structure**
   - Theoretical (DSSP/STRIDE):
     - α-Helix percentage
     - β-Sheet percentage
     - Turn percentage
   - Experimental (UniProt annotations):
     - α-Helix percentage
     - β-Sheet percentage
     - Turn percentage

4. **Charge Profile**
   - Interactive graph of charge vs pH
   - Charge at physiological pH (7.4)
   - Isoelectric point marked

5. **Post-Translational Modifications**
   - Interactive protein sequence viewer
   - PTM annotations
   - Vertical alignment for overlapping PTMs
   - Disulfide bonds highlighted

6. **UniProt Metadata**
   - Gene names and synonyms
   - Protein names (recommended, alternative, short)
   - Functional annotations
   - Cross-references
   - Gene Ontology terms
   - Protein interactions

### Accessing the Report

After running:
```bash
python -m src.main P69905
```

Open `report_P69905.html` in your web browser to view the interactive report.

---

## Examples

### Example 1: Analyze Hemoglobin

```bash
# Generate report
python -m src.main P69905

# Check database
sqlite3 protein_reports.db "SELECT accession_id, protein_name FROM proteins WHERE accession_id = 'P69905';"
```

### Example 2: Analyze p53 (Tumor Suppressor)

```python
from src.main import ProteinReportGenerator

gen = ProteinReportGenerator()
report_path = gen.generate_report('P04637')
print(f"Report: {report_path}")
```

### Example 3: Batch Analysis

```python
from src.services.protein_enrichment_service import enrich_protein_batch

proteins = ['P69905', 'P04637', 'P12345']
profiles = enrich_protein_batch(proteins)

for profile in profiles:
    print(f"{profile.general.uniprot_id}: {profile.general.protein_name}")
    print(f"  MW: {profile.biophysical.molecular_weight_kda:.1f} kDa")
    print(f"  pI: {profile.biophysical.isoelectric_point_pi:.2f}")
    print()
```

### Example 4: Direct Database Query

```python
import sqlite3

conn = sqlite3.connect('protein_reports.db')
cursor = conn.cursor()

# Get all stored proteins
cursor.execute("SELECT accession_id, protein_name FROM proteins")
for acc_id, name in cursor.fetchall():
    print(f"{acc_id}: {name}")

conn.close()
```

### Example 5: Custom Analysis

```python
from src.api_client import APIClient
from src.protein_analyzer import ProteinAnalyzer

client = APIClient()
sequence = client.fetch_protein_sequence('P69905')

analyzer = ProteinAnalyzer(sequence)
properties = analyzer.get_all_properties()

for key, value in properties.items():
    print(f"{key}: {value}")
```

---

## Troubleshooting

### Issue: "ModuleNotFoundError: No module named 'src'"

**Solution**: Make sure you're running from the project root and using the correct Python path:

```bash
# Correct way
python -m src.main P69905

# Also correct
python3 src/main.py P69905

# Incorrect (will fail)
cd src && python main.py P69905
```

### Issue: "Database is locked"

**Solution**: Close any other connections to the database and try again:

```bash
# Check if database file is being used
lsof protein_reports.db

# Delete and recreate if needed
rm protein_reports.db
python -m src.main P69905
```

### Issue: "Connection timeout from UniProt API"

**Solution**: Try again later or increase timeout:

```python
from src.api_client import APIClient

client = APIClient(timeout=30)  # Increase timeout to 30 seconds
```

### Issue: "Template file not found"

**Solution**: Ensure you're running from project root and the file structure is correct:

```bash
# Verify template exists
ls -la templates/offline_template.html

# Verify from correct directory
pwd  # Should show: .../uniprot_lookup
```

### Issue: "No biophysical properties calculated"

**Solution**: Ensure BioPython is installed:

```bash
pip install biopython --upgrade
```

---

## Configuration

### Environment Variables

Create `.env` file (optional):

```bash
cp .env.example .env
```

Edit `.env` if needed. Default configuration is suitable for most use cases.

### Customizing Settings

Edit `src/config.py` to change:

- API endpoints
- Database path
- Report output directory
- Analysis pH range
- PTM categories
- Timeout and retry settings

---

## Performance Notes

- **Single protein analysis**: ~5-10 seconds (includes API calls)
- **Database insert**: <100ms per protein
- **Report generation**: ~2-5 seconds per protein
- **Typical HTML report size**: 150-300 KB

---

## Architecture

See [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) for detailed architecture documentation.

**Key Design Principles:**
- **Modular**: Each module has a single responsibility
- **Scalable**: Easy to add new features
- **Testable**: Clear dependencies and interfaces
- **Maintainable**: Well-documented and organized

**Dependency Graph:**
```
config.py (ROOT - no dependencies)
    ↓
    ├─ api_client.py
    ├─ protein_analyzer.py
    ├─ data_manager.py → protein_database.py, uniprot_metadata_extractor.py
    └─ main.py (orchestrates all modules)
```

---

## Version History

### v2.0 (November 2024)
- Complete repository reorganization
- Modular architecture with src/ structure
- 18 database tables (9 calculated + 9 metadata)
- Comprehensive enrichment service
- Interactive HTML reports with Nightingale viewer
- Documentation: PROJECT_STRUCTURE.md, REORGANIZATION_SUMMARY.md

### v1.0 (Earlier)
- Initial flat project structure
- Basic UniProt fetching and database storage

---

## Contributing

Contributions welcome! Areas for enhancement:
- Additional PTM prediction models
- AlphaFold integration for structure prediction
- Additional visualization options
- Performance optimization for batch processing
- Unit tests

---

## Support

For issues or questions:
1. Check [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) for detailed module documentation
2. Review [REORGANIZATION_SUMMARY.md](REORGANIZATION_SUMMARY.md) for migration guide
3. Check [ARCHITECTURE.md](docs/ARCHITECTURE.md) for architecture details
4. Run troubleshooting commands in [Troubleshooting](#troubleshooting) section

---

## License

MIT License - Feel free to use and modify as needed.

---

## Citation

If you use this tool in research, please cite:

```bibtex
@software{uniprot_lookup_2024,
  title={UniProt Lookup - Comprehensive Protein Analysis Tool},
  year={2024},
  version={2.0}
}
```

---

## Quick Reference

### Common Commands

```bash
# Generate single report
python -m src.main P69905

# Generate multiple reports
python -m src.main P69905 P04637 P12345

# Check Python imports
python3 -c "from src.config import *; from src.api_client import APIClient; print('✓')"

# Check database
sqlite3 protein_reports.db "SELECT COUNT(*) FROM proteins;"

# View stored protein
sqlite3 protein_reports.db "SELECT accession_id, protein_name FROM proteins;"
```

### Common Python Patterns

```python
# One-liner to generate report
from src.main import ProteinReportGenerator; ProteinReportGenerator().generate_report('P69905')

# Get biophysical properties
from src.api_client import APIClient; from src.protein_analyzer import ProteinAnalyzer
seq = APIClient().fetch_protein_sequence('P69905')
props = ProteinAnalyzer(seq).get_all_properties()

# Store and retrieve from database
from src.data_manager import DataManager
dm = DataManager()
dm.save_protein({'accession_id': 'P69905', 'properties': props})
protein = dm.get_protein('P69905')
```

---

**Happy analyzing! 🧬**

For detailed documentation, see:
- [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) - Complete file guide
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) - Architecture details
- [REORGANIZATION_SUMMARY.md](REORGANIZATION_SUMMARY.md) - What changed
