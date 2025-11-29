# Report Generation Guide

Complete guide for generating comprehensive protein analysis reports using UniProt Lookup.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Command Line Usage](#command-line-usage)
3. [Python API Usage](#python-api-usage)
4. [Report Contents](#report-contents)
5. [Customization](#customization)
6. [Advanced Usage](#advanced-usage)
7. [Troubleshooting](#troubleshooting)

---

## Quick Start

### Fastest Way to Generate a Report (30 seconds)

```bash
# 1. Navigate to project directory
cd uniprot_lookup

# 2. Run the tool with a protein accession ID
python -m src.main P69905

# 3. Open the report
open report_P69905.html  # macOS
# or
xdg-open report_P69905.html  # Linux
# or
start report_P69905.html  # Windows
```

Done! You now have a comprehensive protein analysis report.

---

## Command Line Usage

### Single Protein Report

#### Basic Syntax

```bash
python -m src.main <ACCESSION_ID>
```

#### Examples

**Hemoglobin alpha**
```bash
python -m src.main P69905
```

**p53 (tumor suppressor)**
```bash
python -m src.main P04637
```

**Insulin**
```bash
python -m src.main P01308
```

**Output Location**: `report_P69905.html` (in current directory)

### Batch Processing (Multiple Proteins)

#### Generate Reports for Multiple Proteins

```bash
python -m src.main P69905 P04637 P01308
```

This generates:
- `report_P69905.html` (Hemoglobin)
- `report_P04637.html` (p53)
- `report_P01308.html` (Insulin)

#### Processing Order
- Sequential processing (one at a time)
- Each report appears in the output as it completes
- All data stored in database

### Viewing Generated Reports

#### Open in Web Browser

```bash
# macOS
open report_P69905.html

# Linux
firefox report_P69905.html
xdg-open report_P69905.html

# Windows
start report_P69905.html

# WSL (Windows Subsystem for Linux)
explorer.exe report_P69905.html
```

#### Check Files

```bash
# List all generated reports
ls -lh report_*.html

# Check file size
du -h report_P69905.html

# Check creation time
ls -lT report_P69905.html
```

---

## Python API Usage

### Method 1: Using ProteinReportGenerator (Recommended)

```python
from src.main import ProteinReportGenerator

# Initialize
generator = ProteinReportGenerator()

# Generate single report
report_path = generator.generate_report('P69905')
print(f"Report saved to: {report_path}")
```

### Method 2: Direct Enrichment Service

```python
from src.services.protein_enrichment_service import enrich_protein

# Get enrichment profile
profile = enrich_protein('P69905')

# Export to JSON
json_data = profile.to_json()
print(json_data)

# Or save to file
profile.to_json('enriched_P69905.json')
```

### Method 3: Complete Workflow

```python
from src.api_client import APIClient
from src.protein_analyzer import ProteinAnalyzer
from src.data_manager import DataManager
from src.services.protein_enrichment_service import enrich_protein

accession = 'P69905'

# Step 1: Fetch data
client = APIClient()
all_data = client.fetch_all_protein_data(accession)
print(f"✓ Fetched sequence: {len(all_data['sequence'])} aa")

# Step 2: Analyze
analyzer = ProteinAnalyzer(all_data['sequence'])
properties = analyzer.get_all_properties()
print(f"✓ Calculated {len(properties)} properties")

# Step 3: Enrich
profile = enrich_protein(accession)
print(f"✓ Created enrichment profile")

# Step 4: Store
dm = DataManager()
dm.save_protein({
    'accession_id': accession,
    'sequence': all_data['sequence'],
    'properties': properties
})
print(f"✓ Stored in database")

# Step 5: Export
profile.to_json(f'report_{accession}.json')
print(f"✓ Exported JSON report")
```

### Method 4: Batch Processing with Python

```python
from src.main import ProteinReportGenerator

generator = ProteinReportGenerator()

# Batch process
proteins = ['P69905', 'P04637', 'P01308']
generator.generate_batch_reports(proteins)

# Or manually
for protein_id in proteins:
    try:
        report_path = generator.generate_report(protein_id)
        print(f"✓ {protein_id}: {report_path}")
    except Exception as e:
        print(f"✗ {protein_id}: {e}")
```

### Method 5: Custom Report Generation with Control

```python
from src.main import ProteinReportGenerator
import json

class CustomReportGenerator(ProteinReportGenerator):
    def generate_report(self, accession_id: str, custom_config=None) -> bool:
        """Generate report with custom configuration."""
        try:
            # Fetch data
            protein_data = self.api_client.fetch_all_protein_data(accession_id)
            if not protein_data['sequence']:
                print(f"❌ No sequence for {accession_id}")
                return False

            # Enrich
            profile = self.enrichment_service.enrich_protein(accession_id)

            # Save to database
            self.data_manager.save_protein({
                'accession_id': accession_id,
                'sequence': protein_data['sequence'],
                'features': protein_data['features'],
            })

            # Custom output
            if custom_config:
                if custom_config.get('save_json'):
                    profile.to_json(f"enriched_{accession_id}.json")
                if custom_config.get('print_summary'):
                    print(f"\n{accession_id}:")
                    print(f"  Name: {profile.general.protein_name}")
                    print(f"  MW: {profile.biophysical.molecular_weight_kda:.1f} kDa")
                    print(f"  pI: {profile.biophysical.isoelectric_point_pi:.2f}")

            return True
        except Exception as e:
            print(f"Error: {e}")
            return False

# Use it
generator = CustomReportGenerator()
generator.generate_report('P69905', custom_config={
    'save_json': True,
    'print_summary': True
})
```

---

## Report Contents

### What's Included in Each Report

#### 1. Header Information
- Protein name
- Accession ID
- Gene name
- Organism
- Sequence length
- PDB ID (if available)

#### 2. Sequence Viewer
- Interactive protein sequence display
- Nightingale protein viewer integration
- PTM annotations directly on sequence
- Disulfide bonds highlighted
- Click to navigate

#### 3. Biophysical Properties
- **Molecular Weight**: Da and kDa
- **Isoelectric Point (pI)**: Charge at zero
- **GRAVY**: Hydropathy score
- **Aromaticity**: Aromatic amino acid content
- **Instability Index**: Protein stability
- **Cysteine Count**: For disulfide bond potential
- **Extinction Coefficients**: UV absorption (280 nm)

#### 4. Secondary Structure
- **α-Helix**: Percentage
- **β-Sheet**: Percentage
- **Coil**: Percentage
- Based on Chou-Fasman algorithm

#### 5. Charge Profile
- **Interactive Graph**: Charge vs pH
- **pH Range**: 2.0 - 12.0
- **Isoelectric Point**: Marked on graph
- **Charge at pH 7.4**: Physiological pH

#### 6. Post-Translational Modifications
- **Total PTM Count**: Number of modification sites
- **PTM Types**: MOD_RES, CARBOHYD, DISULFID, etc.
- **Interactive Visualization**: Click for details
- **Position Information**: Exact residue locations

#### 7. UniProt Metadata
- Gene names and synonyms
- Protein name variants
- Functional annotations
- Cross-references to other databases
- Gene Ontology annotations
- Keywords and topics
- Protein interactions
- Literature references

---

## Customization

### Changing Report Output Directory

#### Option 1: Edit Configuration

```python
# In src/config.py
REPORT = {
    'template': 'templates/offline_template.html',
    'output_dir': '/path/to/reports',  # Change this
    'prefix': 'report_',
}
```

#### Option 2: Programmatically

```python
from src.main import ProteinReportGenerator
from src.config import REPORT

# Change before generating
REPORT['output_dir'] = '/tmp/my_reports'

generator = ProteinReportGenerator()
report_path = generator.generate_report('P69905')
```

### Changing Report Filename

```python
from src.main import ProteinReportGenerator
from src.config import REPORT

# Change prefix
REPORT['prefix'] = 'analysis_'  # Results in: analysis_P69905.html

generator = ProteinReportGenerator()
generator.generate_report('P69905')
```

### Adding Custom CSS/Styling

Edit `templates/offline_template.html` to customize appearance:

```html
<!-- Add custom CSS -->
<style>
    .protein-info {
        background-color: #f0f0f0;
        padding: 20px;
        border-radius: 5px;
    }

    .property-value {
        font-weight: bold;
        color: #0066cc;
    }
</style>
```

### Modifying Report Template

The template is located at `templates/offline_template.html`.

Key template variables available:
- `protein_name`
- `accession_id`
- `sequence_length`
- `molecular_weight`
- `isoelectric_point`
- `secondary_structure`
- `ptm_summary`
- `charge_profile`

---

## Advanced Usage

### Processing with Error Handling

```python
from src.main import ProteinReportGenerator
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)

# Create generator
generator = ProteinReportGenerator()

# Process with error handling
proteins = ['P69905', 'INVALID_ID', 'P04637']
results = {
    'success': [],
    'failed': []
}

for protein_id in proteins:
    try:
        report_path = generator.generate_report(protein_id)
        results['success'].append({
            'id': protein_id,
            'report': report_path
        })
        print(f"✓ {protein_id}")
    except Exception as e:
        results['failed'].append({
            'id': protein_id,
            'error': str(e)
        })
        print(f"✗ {protein_id}: {e}")

# Summary
print(f"\nSummary:")
print(f"  Success: {len(results['success'])}")
print(f"  Failed: {len(results['failed'])}")

if results['failed']:
    print("\nFailed proteins:")
    for item in results['failed']:
        print(f"  - {item['id']}: {item['error']}")
```

### Database Integration

```python
from src.main import ProteinReportGenerator
from src.data_manager import DataManager
import sqlite3

# Generate reports
generator = ProteinReportGenerator()
proteins = ['P69905', 'P04637']

for protein_id in proteins:
    generator.generate_report(protein_id)

# Query database
dm = DataManager()
conn = sqlite3.connect('protein_reports.db')
cursor = conn.cursor()

# Get all stored proteins
cursor.execute("""
    SELECT accession_id, protein_name, molecular_weight, isoelectric_point
    FROM proteins
""")

for row in cursor.fetchall():
    acc_id, name, mw, pi = row
    print(f"{acc_id}: {name}")
    print(f"  MW: {mw/1000:.1f} kDa, pI: {pi:.2f}")

conn.close()
```

### Export to Different Formats

```python
from src.services.protein_enrichment_service import enrich_protein
import json
import csv

protein_id = 'P69905'
profile = enrich_protein(protein_id)

# Export as JSON
profile.to_json(f'{protein_id}_enriched.json')

# Export key properties as CSV
with open(f'{protein_id}_properties.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Property', 'Value'])
    writer.writerow(['Protein Name', profile.general.protein_name])
    writer.writerow(['Accession ID', profile.general.uniprot_id])
    writer.writerow(['Sequence Length', profile.sequence.length])
    writer.writerow(['Molecular Weight (Da)', profile.biophysical.molecular_weight_da])
    writer.writerow(['Isoelectric Point', profile.biophysical.isoelectric_point_pi])
    writer.writerow(['GRAVY Score', profile.biophysical.gravy])
    writer.writerow(['Instability Index', profile.biophysical.instability_index])

# Export to plain text
with open(f'{protein_id}_report.txt', 'w') as f:
    f.write(f"Protein: {profile.general.protein_name}\n")
    f.write(f"Accession: {profile.general.uniprot_id}\n")
    f.write(f"Sequence Length: {profile.sequence.length}\n")
    f.write(f"MW: {profile.biophysical.molecular_weight_kda:.1f} kDa\n")
    f.write(f"pI: {profile.biophysical.isoelectric_point_pi:.2f}\n")
    f.write(f"Stability: {'Stable' if profile.biophysical.instability_index < 40 else 'Unstable'}\n")
```

### Scheduled/Automated Reporting

```python
from src.main import ProteinReportGenerator
from datetime import datetime
import schedule
import time

# Define job
def generate_daily_reports():
    proteins = ['P69905', 'P04637', 'P01308']
    generator = ProteinReportGenerator()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"[{timestamp}] Starting report generation...")

    for protein_id in proteins:
        try:
            generator.generate_report(protein_id)
            print(f"[{timestamp}] ✓ {protein_id}")
        except Exception as e:
            print(f"[{timestamp}] ✗ {protein_id}: {e}")

    print(f"[{timestamp}] Completed")

# Schedule (requires: pip install schedule)
schedule.every().day.at("09:00").do(generate_daily_reports)

# Run scheduler
while True:
    schedule.run_pending()
    time.sleep(60)
```

---

## Troubleshooting

### Report Not Generated

**Issue**: No HTML file created

**Solutions**:
```bash
# 1. Check if accession ID is valid
curl "https://rest.uniprot.org/uniprotkb/P69905.json" | head

# 2. Verify permissions
ls -la report_*.html

# 3. Try verbose output
python -m src.main P69905 -v

# 4. Check database
sqlite3 protein_reports.db "SELECT COUNT(*) FROM proteins;"
```

### "Connection timeout" Error

**Issue**: API request takes too long

**Solution**:
```python
from src.api_client import APIClient

# Increase timeout
client = APIClient(timeout=30)  # 30 seconds instead of 10

sequence = client.fetch_protein_sequence('P69905')
```

### "Template not found" Error

**Issue**: HTML template file missing

**Solutions**:
```bash
# 1. Verify file exists
ls -la templates/offline_template.html

# 2. Check from correct directory
pwd  # Should be: .../uniprot_lookup

# 3. Verify path in config
grep -n "template" src/config.py
```

### Report File Size Very Large

**Issue**: HTML file is 300+ MB

**Causes**:
- Very large protein sequences
- Excessive PTM data
- Template includes large libraries

**Solution**:
```python
# Use JSON export instead
profile = enrich_protein('P69905')
profile.to_json('report.json')  # Much smaller

# Or compress HTML
import gzip
with open('report.html', 'rb') as f_in:
    with gzip.open('report.html.gz', 'wb') as f_out:
        f_out.write(f_in.read())
```

### BioPython Not Installed

**Issue**: "No module named 'Bio'"

**Solution**:
```bash
pip install biopython
pip install biopython --upgrade

# Verify
python -c "from Bio.SeqUtils.ProtParam import ProteinAnalysis; print('✓')"
```

### Database Locked

**Issue**: "Database is locked" error

**Solutions**:
```bash
# 1. Close other connections
# Check if other processes are using it
lsof protein_reports.db

# 2. Delete and recreate
rm protein_reports.db
python -m src.main P69905

# 3. Use temporary database
python -c "from src.services.protein_database import ProteinDatabase; db = ProteinDatabase(':memory:'); print('✓')"
```

### Invalid Accession ID

**Issue**: "Failed to fetch" error

**Solution**:
```bash
# Verify accession ID format
# UniProt accessions are typically: P|Q followed by 5 digits

# Check online
curl "https://rest.uniprot.org/uniprotkb/P69905.json"

# Try alternative accession
python -m src.main P04637  # Known valid ID
```

### Out of Memory

**Issue**: Process killed during batch processing

**Solution**:
```python
# Process one at a time with cleanup
import gc
from src.main import ProteinReportGenerator

generator = ProteinReportGenerator()

for protein_id in proteins:
    generator.generate_report(protein_id)
    gc.collect()  # Force garbage collection
```

---

## Performance Tips

### Speed Up Report Generation

```python
# Use batch API calls if available
from src.api_client import APIClient

client = APIClient()
all_proteins = ['P69905', 'P04637', 'P01308']

# Parallel processing (requires: pip install joblib)
from joblib import Parallel, delayed

results = Parallel(n_jobs=4)(
    delayed(client.fetch_all_protein_data)(pid)
    for pid in all_proteins
)
```

### Reduce Database Size

```bash
# Clean up old reports
sqlite3 protein_reports.db "DELETE FROM generated_reports WHERE generated_at < date('now', '-30 days');"

# Vacuum database
sqlite3 protein_reports.db "VACUUM;"

# Check size
du -h protein_reports.db
```

### Cache Frequently Used Data

```python
from functools import lru_cache
from src.api_client import APIClient

# Create cached version
@lru_cache(maxsize=100)
def cached_fetch(accession):
    client = APIClient()
    return client.fetch_uniprot_entry(accession)

# Use cached version
data = cached_fetch('P69905')
```

---

## Best Practices

1. **Always use correct project directory**
   ```bash
   cd uniprot_lookup  # Must be in this directory
   python -m src.main P69905
   ```

2. **Batch process with error handling**
   ```python
   for protein in proteins:
       try:
           generate_report(protein)
       except Exception as e:
           log.error(f"Failed: {protein}: {e}")
           continue
   ```

3. **Verify database before querying**
   ```python
   from src.data_manager import DataManager
   dm = DataManager()
   proteins = dm.list_proteins()
   if proteins:
       # Process...
   ```

4. **Keep accession IDs in separate file**
   ```bash
   # proteins.txt
   P69905
   P04637
   P01308

   # Read and process
   python -c "
   from src.main import ProteinReportGenerator
   gen = ProteinReportGenerator()
   with open('proteins.txt') as f:
       for pid in f:
           gen.generate_report(pid.strip())
   "
   ```

5. **Backup database regularly**
   ```bash
   cp protein_reports.db protein_reports.db.backup
   ```

---

## Summary

**Simplest Way**:
```bash
python -m src.main P69905
open report_P69905.html
```

**Batch Way**:
```bash
python -m src.main P69905 P04637 P01308
```

**Programmatic Way**:
```python
from src.main import ProteinReportGenerator
ProteinReportGenerator().generate_report('P69905')
```

**Database Query Way**:
```bash
sqlite3 protein_reports.db "SELECT accession_id, protein_name FROM proteins;"
```

---

**Need help?** Check [README_COMPREHENSIVE.md](README_COMPREHENSIVE.md) for full API reference.
