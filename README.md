# UniProt Lookup - Protein Analysis & Report Generation

A modular Python application for comprehensive protein analysis using UniProt data, biophysical calculations, and interactive HTML report generation.

**Version**: 2.0 | **Status**: Production Ready ✓ | **Size**: ~186 KB

---

## 🚀 Quick Start (30 seconds)

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Generate a report
python -m src.main P69905

# 3. View the report
open report_P69905.html
```

Done! You now have a comprehensive protein analysis report.

---

## 📚 Documentation

Choose your path based on what you need:

### For New Users
**Start here**: [README_COMPREHENSIVE.md](README_COMPREHENSIVE.md)
- Overview of features
- Installation guide
- Complete API reference
- Examples and best practices

### For Generating Reports
**Go here**: [REPORT_GENERATION_GUIDE.md](REPORT_GENERATION_GUIDE.md)
- Command-line usage
- Python API examples
- Report customization
- Troubleshooting

### For Technical Details
**Check here**: [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)
- File organization
- Module descriptions
- Dependency graph
- Configuration guide

---

## ✨ Features

✓ **UniProt Data Fetching** - Sequence, metadata, features, annotations
✓ **Biophysical Analysis** - MW, pI, GRAVY, instability index, secondary structure
✓ **PTM Analysis** - Post-translational modifications with visualization
✓ **Database Storage** - SQLite with 18 normalized tables
✓ **Interactive Reports** - Professional HTML with Nightingale viewer
✓ **Batch Processing** - Generate multiple reports at once
✓ **CLI & Python API** - Use from command line or import as library

---

## 📋 Project Structure

```
uniprot_lookup/
├── src/                              # Application code
│   ├── config.py                     # Configuration
│   ├── api_client.py                 # API interactions
│   ├── protein_analyzer.py           # Biophysical analysis
│   ├── data_manager.py               # Data management
│   ├── main.py                       # Entry point
│   └── services/
│       ├── protein_database.py       # Database operations
│       ├── protein_enrichment_service.py
│       └── uniprot_metadata_extractor.py
│
├── data/                             # Database schema
│   └── database_schema.sql
│
├── templates/                        # HTML templates
│   └── offline_template.html
│
├── README_COMPREHENSIVE.md           # Complete guide
├── REPORT_GENERATION_GUIDE.md        # Report generation
├── PROJECT_STRUCTURE.md              # Technical details
├── requirements.txt                  # Dependencies
└── .env.example                      # Configuration template
```

---

## 🎯 Common Tasks

### Generate Report for Single Protein

```bash
python -m src.main P69905
```

### Generate Reports for Multiple Proteins

```bash
python -m src.main P69905 P04637 P01308
```

### Use as Python Library

```python
from src.main import ProteinReportGenerator

gen = ProteinReportGenerator()
gen.generate_report('P69905')
```

### Query Database

```bash
sqlite3 protein_reports.db "SELECT accession_id, protein_name FROM proteins;"
```

---

## 🔧 Installation

### Prerequisites
- Python 3.8+
- pip

### Steps

```bash
# Clone/navigate to project
cd uniprot_lookup

# Install dependencies
pip install -r requirements.txt

# Verify installation
python3 -c "from src.config import *; from src.api_client import APIClient; print('✓ Ready!')"
```

---

## 📖 API Quick Reference

### Generate Reports

```python
from src.main import ProteinReportGenerator

gen = ProteinReportGenerator()
gen.generate_report('P69905')              # Single protein
gen.generate_batch_reports(['P04637', 'P15018'])  # Batch
```

### Fetch Data

```python
from src.api_client import APIClient

client = APIClient()
sequence = client.fetch_protein_sequence('P69905')
entry = client.fetch_uniprot_entry('P69905')
all_data = client.fetch_all_protein_data('P69905')
```

### Analyze Proteins

```python
from src.protein_analyzer import ProteinAnalyzer

analyzer = ProteinAnalyzer(sequence)
properties = analyzer.get_all_properties()
print(f"MW: {properties['molecular_weight']} Da")
print(f"pI: {properties['isoelectric_point']}")
```

### Store & Retrieve

```python
from src.data_manager import DataManager

dm = DataManager()
dm.save_protein(data)
protein = dm.get_protein('P69905')
all_proteins = dm.list_proteins()
```

---

## 📊 Database

- **Type**: SQLite
- **File**: `protein_reports.db`
- **Tables**: 18 normalized tables
- **Size**: Typical ~100 KB per 100 proteins

**Includes:**
- Biophysical properties (MW, pI, GRAVY, etc.)
- Secondary structure predictions
- Charge profiles (pH-dependent)
- PTM information
- UniProt metadata

---

## 🐛 Troubleshooting

### "ModuleNotFoundError: No module named 'src'"
```bash
# Make sure you're running from project root
cd uniprot_lookup
python -m src.main P69905
```

### "Failed to fetch protein data"
```bash
# Verify accession ID (e.g., P69905, P04637)
# Check internet connection
# Try a known ID: P04637 (p53)
```

### "BioPython not installed"
```bash
pip install biopython --upgrade
```

See [REPORT_GENERATION_GUIDE.md](REPORT_GENERATION_GUIDE.md#troubleshooting) for more solutions.

---

## 📚 Documentation Guide

| Document | Purpose | When to Use |
|----------|---------|------------|
| **README_COMPREHENSIVE.md** | Complete guide with everything | First time setup, API reference |
| **REPORT_GENERATION_GUIDE.md** | Report generation focused | Generating reports, examples |
| **PROJECT_STRUCTURE.md** | Technical architecture | Understanding code structure |

---

## 💡 Examples

### Example 1: Hemoglobin Analysis

```bash
python -m src.main P69905
open report_P69905.html
```

### Example 2: Batch Analysis

```python
from src.main import ProteinReportGenerator

proteins = ['P69905', 'P04637', 'P01308']  # Hemoglobin, p53, Insulin
gen = ProteinReportGenerator()

for protein_id in proteins:
    gen.generate_report(protein_id)
    print(f"✓ {protein_id}")
```

### Example 3: Custom Analysis

```python
from src.api_client import APIClient
from src.protein_analyzer import ProteinAnalyzer

client = APIClient()
sequence = client.fetch_protein_sequence('P69905')

analyzer = ProteinAnalyzer(sequence)
print(f"Sequence: {sequence}")
print(f"Length: {len(sequence)} aa")
print(f"MW: {analyzer.get_molecular_weight():.2f} Da")
print(f"pI: {analyzer.get_isoelectric_point():.2f}")
print(f"Stability: {analyzer.is_stable()}")
```

---

## 📦 Dependencies

- `requests` - HTTP API calls
- `biopython` - Biophysical calculations
- `jinja2` - Template rendering (optional)

Install all with:
```bash
pip install -r requirements.txt
```

---

## 🎓 Learning Path

**Beginner**:
1. Read this README
2. Run quick start: `python -m src.main P69905`
3. Check [README_COMPREHENSIVE.md](README_COMPREHENSIVE.md#quick-start)

**Intermediate**:
1. Follow [REPORT_GENERATION_GUIDE.md](REPORT_GENERATION_GUIDE.md)
2. Try Python API examples
3. Generate batch reports

**Advanced**:
1. Study [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)
2. Customize templates and modules
3. Extend with new features

---

## 🤝 Contributing

Areas for enhancement:
- Additional PTM prediction models
- Performance optimization
- Additional visualization options
- Unit tests

---

## 📝 License

MIT License - Free to use and modify

---

## 🆘 Support

**Need help?**
1. Check the relevant documentation guide above
2. Review [REPORT_GENERATION_GUIDE.md#troubleshooting](REPORT_GENERATION_GUIDE.md#troubleshooting)
3. Check Python API in [README_COMPREHENSIVE.md](README_COMPREHENSIVE.md#api-reference)

**Still stuck?**
- Verify installation: `python3 -c "from src.config import *; print('✓')"`
- Check internet connection for UniProt API access
- Ensure Python 3.8+ is installed

---

## ⚡ Quick Commands

```bash
# Generate single report
python -m src.main P69905

# Generate multiple reports
python -m src.main P69905 P04637 P01308

# Check stored proteins
sqlite3 protein_reports.db "SELECT accession_id, protein_name FROM proteins;"

# Verify setup
python3 -c "from src.config import *; from src.api_client import APIClient; print('✓')"
```

---

**Happy analyzing! 🧬**

Start with: [README_COMPREHENSIVE.md](README_COMPREHENSIVE.md)
