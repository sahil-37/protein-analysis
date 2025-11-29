# UniProt Lookup & Database Storage

This module fetches protein data from the UniProt REST API and stores it in a relational database with a normalized schema.

**Database Options:**
- **SQLite** (default) - Zero configuration, perfect for development and small-scale use
- **PostgreSQL** (optional) - For production and large-scale deployments

## Features

### UniProt Data Fetching & Storage
- **Fetch protein data** from UniProt API by accession ID
- **Zero-setup option** - Uses SQLite by default (no database installation needed!)
- **Flatten complex JSON** into relational database tables
- **Comprehensive schema** covering:
  - Main protein information
  - Gene names and synonyms
  - Protein features (domains, sites, regions)
  - Comments (function, catalytic activity, etc.)
  - Cross-references to other databases
  - Gene Ontology (GO) terms
  - Keywords
  - Publications/References
- **Batch processing** support for multiple proteins
- **Error handling** with retries and rate limiting
- **Idempotent operations** - safe to re-run

### Protein Design Analysis (NEW!)
- **Isoelectric Point (pI) calculation** using Henderson-Hasselbalch equation
- **Secondary structure prediction** using Chou-Fasman algorithm
- **Net charge calculation** at any specified pH
- **Expression & purification recommendations**:
  - Ion-exchange chromatography strategy (cation vs anion)
  - Optimal buffer pH ranges
  - Solubility and aggregation risk assessment
  - Fusion tag recommendations (His, GST, MBP, NusA)
  - Temperature and growth media suggestions
  - E. coli strain recommendations
  - Protein refolding strategies if needed

## Project Structure

```
uniprot_lookup/
├── main.py                         # Main orchestration module
├── uniprot_client.py               # UniProt API client
├── data_processor.py               # JSON flattening logic
├── database.py                     # PostgreSQL database handler
├── database_sqlite.py              # SQLite database handler (default)
├── protein_design_analyzer.py      # Protein design & analysis tool (NEW!)
├── schema.sql                      # PostgreSQL schema definition
├── requirements.txt                # Python dependencies
├── .env.example                    # Example environment configuration
├── QUICKSTART.md                   # Quick start guide for UniProt lookup
├── PROTEIN_DESIGN_QUICKSTART.md    # Quick start for protein design analyzer
├── PROTEIN_DESIGN_GUIDE.md         # Detailed guide for protein design
├── README.md                       # This file
└── VISUALIZATION_GUIDE.md          # Guide for data visualization
```

## Quick Start (SQLite - No Setup Required!)

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Initialize database
python main.py --init

# 3. Fetch a protein
python main.py P12345
```

That's it! The data is now stored in `uniprot.db`.

See [QUICKSTART.md](QUICKSTART.md) for more examples and queries.

## Full Installation Guide

### Option 1: SQLite (Default - Recommended for Getting Started)

1. **Install Dependencies**
   ```bash
   pip install -r requirements.txt
   ```

2. **Initialize Database**
   ```bash
   python main.py --init
   ```

3. **Start Using**
   ```bash
   python main.py P04637
   ```

### Option 2: PostgreSQL (For Production Use)

1. **Install Dependencies**
   ```bash
   pip install -r requirements.txt
   ```

2. **Set Up PostgreSQL Database**
   ```bash
   createdb uniprot_db
   ```

3. **Configure Environment**
   ```bash
   cp .env.example .env
   # Edit .env with your PostgreSQL credentials
   ```

4. **Initialize Database Schema**
   ```bash
   python main.py --init
   ```

## Usage

### Command Line

#### Fetch a single protein:

```bash
python main.py P12345
```

#### Fetch multiple proteins:

```bash
python main.py P12345 Q9Y6K9 P04637
```

### As a Python Library

```python
from main import fetch_and_save_protein, fetch_and_save_multiple_proteins

# Fetch and save a single protein
success = fetch_and_save_protein('P12345')

# Fetch and save multiple proteins
accessions = ['P12345', 'Q9Y6K9', 'P04637']
results = fetch_and_save_multiple_proteins(accessions)

print(f"Successful: {results['successful']}")
print(f"Failed: {results['failed']}")
```

### Custom Database Configuration

```python
from main import fetch_and_save_protein

db_config = {
    'host': 'localhost',
    'port': '5432',
    'dbname': 'my_database',
    'user': 'my_user',
    'password': 'my_password'
}

fetch_and_save_protein('P12345', db_config=db_config)
```

### Using Individual Components

#### Fetch data from UniProt API:

```python
from uniprot_client import UniProtClient

client = UniProtClient()
protein_data = client.fetch_protein('P12345')
```

#### Process and flatten JSON data:

```python
from data_processor import UniProtDataProcessor

processor = UniProtDataProcessor(protein_data)
flattened = processor.get_all_flattened_data()

# Access specific data
print(flattened['protein'])          # Main protein info
print(flattened['features'])         # Protein features
print(flattened['go_terms'])         # GO terms
```

#### Save to database:

```python
from database import UniProtDatabase

with UniProtDatabase() as db:
    db.save_protein(flattened)
```

## Database Schema

The database consists of 10 main tables:

1. **proteins** - Main protein information
2. **protein_secondary_accessions** - Alternative accession IDs
3. **protein_genes** - Gene names and synonyms
4. **protein_names** - Protein names (recommended, alternative, short)
5. **protein_features** - Sequence features (domains, sites, etc.)
6. **protein_comments** - Functional annotations
7. **protein_cross_references** - Links to other databases
8. **protein_keywords** - Keywords
9. **protein_go_terms** - Gene Ontology annotations
10. **protein_references** - Publications

### Example Queries

```sql
-- Get all proteins for a specific organism
SELECT accession, protein_recommended_name, gene_name
FROM proteins
WHERE organism_taxon_id = 9606;  -- Human

-- Get all GO terms for a protein
SELECT go_id, go_term, go_aspect
FROM protein_go_terms
WHERE accession = 'P12345';

-- Find proteins with specific features
SELECT DISTINCT p.accession, p.protein_recommended_name
FROM proteins p
JOIN protein_features f ON p.accession = f.accession
WHERE f.feature_type = 'Domain';

-- Get all cross-references to PDB
SELECT accession, database_id
FROM protein_cross_references
WHERE database_name = 'PDB';
```

## API Reference

### UniProtClient

```python
client = UniProtClient(max_retries=3, retry_delay=2)

# Fetch single protein
data = client.fetch_protein('P12345')

# Fetch multiple proteins
results = client.fetch_multiple_proteins(['P12345', 'Q9Y6K9'])
```

### UniProtDataProcessor

```python
processor = UniProtDataProcessor(uniprot_json_data)

# Get all flattened data
all_data = processor.get_all_flattened_data()

# Or get specific sections
protein_info = processor.get_main_protein_data()
features = processor.get_features()
go_terms = processor.get_go_terms()
```

### UniProtDatabase

```python
# Using context manager (recommended)
with UniProtDatabase() as db:
    db.save_protein(flattened_data)
    exists = db.protein_exists('P12345')
    protein = db.get_protein('P12345')

# Manual connection management
db = UniProtDatabase()
db.connect()
db.save_protein(flattened_data)
db.disconnect()
```

## Error Handling

The module includes comprehensive error handling:

- **Network errors**: Automatic retries with exponential backoff
- **Rate limiting**: Respects API rate limits with delays
- **Database errors**: Transaction rollback on failures
- **Missing data**: Graceful handling of optional fields

## Examples

### Example 1: Fetch Human TP53 (tumor suppressor p53)

```bash
python main.py P04637
```

### Example 2: Batch Processing

```python
from main import fetch_and_save_multiple_proteins

# Fetch multiple cancer-related proteins
cancer_proteins = [
    'P04637',  # TP53
    'P01112',  # KRAS
    'P42336',  # PIK3CA
    'Q9Y6K9',  # NF2
]

results = fetch_and_save_multiple_proteins(cancer_proteins)
```

### Example 3: Check if Protein Exists Before Fetching

```python
from database import UniProtDatabase
from main import fetch_and_save_protein

with UniProtDatabase() as db:
    if not db.protein_exists('P12345'):
        fetch_and_save_protein('P12345', skip_if_exists=False)
```

## Troubleshooting

### Database Connection Errors

Ensure PostgreSQL is running and credentials in `.env` are correct:

```bash
psql -U postgres -d uniprot_db -c "SELECT 1;"
```

### API Rate Limiting

If you encounter rate limiting, increase the retry delay:

```python
client = UniProtClient(max_retries=5, retry_delay=5)
```

### Missing Dependencies

Reinstall dependencies:

```bash
pip install -r requirements.txt --upgrade
```

## License

MIT License

## Contributing

Contributions welcome! Please submit issues and pull requests.

## Protein Design Analysis

The `protein_design_analyzer.py` tool provides biophysical analysis for protein expression optimization.

### Quick Start

```bash
# Analyze a protein
python protein_design_analyzer.py P04637

# With options
python protein_design_analyzer.py P04637 --ph 7.4 --pka-table lehninger --show-ss-string
```

### Features

- **Isoelectric Point (pI)**: Calculates protein charge at zero using Henderson-Hasselbalch equation
- **Secondary Structure**: Predicts % helix, % sheet, % coil using Chou-Fasman algorithm
- **Net Charge**: Calculates protein charge at any specified pH
- **Design Suggestions**: Automatic recommendations for expression and purification

### Documentation

- [PROTEIN_DESIGN_QUICKSTART.md](PROTEIN_DESIGN_QUICKSTART.md) - Quick reference
- [PROTEIN_DESIGN_GUIDE.md](PROTEIN_DESIGN_GUIDE.md) - Comprehensive guide with theory and examples

### Example

```bash
$ python protein_design_analyzer.py P04637 --ph 7.4

Accession ID:           P04637
Protein Name:           Cellular tumor antigen p53
Protein Length:         393 amino acids
Isoelectric Point (pI):  6.37
Net Charge at pH 7.4:    -4.69

SECONDARY STRUCTURE:
  α-helix:    32.3%
  β-sheet:    23.9%
  Coil:       43.8%

DESIGN SUGGESTIONS:
  → Use CATION EXCHANGE chromatography
  → Bind at pH 6.5-7.5, elute with 0.1-0.3 M NaCl
  → Good solubility expected (high coil content)
  → Consider His-tag or MBP-tag for purification
```

## References

### UniProt Lookup
- [UniProt REST API Documentation](https://www.uniprot.org/help/api)
- [UniProt Data Structure](https://www.uniprot.org/help/uniprotkb)

### Protein Design Analysis
- **Henderson-Hasselbalch**: Hasselbalch, K. A. (1917)
- **Chou-Fasman**: Chou, P. Y., & Fasman, G. D. (1978). Empirical predictions of protein conformation
- **pKa Values**: EMBOSS and Lehninger Principles of Biochemistry
