-- ============================================================================
-- UniProt Protein Database Schema
-- Stores all protein data used to generate reports
-- ============================================================================

-- Main proteins table
CREATE TABLE IF NOT EXISTS proteins (
    protein_id INTEGER PRIMARY KEY AUTOINCREMENT,
    accession_id TEXT NOT NULL UNIQUE,
    protein_name TEXT,
    gene_name TEXT,
    organism TEXT,
    subcellular_location TEXT,
    sequence TEXT NOT NULL,
    sequence_length INTEGER NOT NULL,
    molecular_weight REAL,
    isoelectric_point REAL,
    gravy_score REAL,
    aromaticity REAL,
    instability_index REAL,
    cysteine_count INTEGER,
    extinction_coeff_reduced INTEGER,
    extinction_coeff_oxidized INTEGER,
    pdb_id TEXT,
    has_structure BOOLEAN DEFAULT 0,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Secondary structure data (predicted vs experimental)
CREATE TABLE IF NOT EXISTS secondary_structure (
    structure_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    helix_predicted REAL,
    turn_predicted REAL,
    sheet_predicted REAL,
    helix_experimental REAL,
    turn_experimental REAL,
    sheet_experimental REAL,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Charge profile at different pH values
CREATE TABLE IF NOT EXISTS charge_profile (
    charge_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    ph_value REAL NOT NULL,
    charge REAL NOT NULL,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Post-translational modifications
CREATE TABLE IF NOT EXISTS ptm_summary (
    ptm_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    total_ptms INTEGER,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Individual PTM details
CREATE TABLE IF NOT EXISTS ptm_details (
    ptm_detail_id INTEGER PRIMARY KEY AUTOINCREMENT,
    ptm_id INTEGER NOT NULL,
    position INTEGER NOT NULL,
    ptm_type TEXT NOT NULL,
    description TEXT,
    feature_type TEXT,
    FOREIGN KEY (ptm_id) REFERENCES ptm_summary(ptm_id) ON DELETE CASCADE
);

-- PTM type counts
CREATE TABLE IF NOT EXISTS ptm_type_counts (
    ptm_type_id INTEGER PRIMARY KEY AUTOINCREMENT,
    ptm_id INTEGER NOT NULL,
    ptm_type TEXT NOT NULL,
    count INTEGER,
    FOREIGN KEY (ptm_id) REFERENCES ptm_summary(ptm_id) ON DELETE CASCADE
);

-- All protein features (domains, regions, sites, etc.)
CREATE TABLE IF NOT EXISTS protein_features (
    feature_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    feature_type TEXT NOT NULL,
    category TEXT,
    description TEXT,
    start_position INTEGER,
    end_position INTEGER,
    feature_length INTEGER,
    ft_id TEXT,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Disulfide bonds
CREATE TABLE IF NOT EXISTS disulfide_bonds (
    disulfide_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    position_1 INTEGER,
    position_2 INTEGER,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Reports generated for each protein
CREATE TABLE IF NOT EXISTS generated_reports (
    report_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    report_filename TEXT NOT NULL,
    report_path TEXT,
    file_size_kb INTEGER,
    generated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- ============================================================================
-- UniProt Metadata Tables (comprehensive data collection)
-- ============================================================================

-- Secondary accessions (alternative accession IDs)
CREATE TABLE IF NOT EXISTS protein_secondary_accessions (
    secondary_accession_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    secondary_accession TEXT,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Gene names (including synonyms and different gene types)
CREATE TABLE IF NOT EXISTS protein_genes (
    gene_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    gene_name TEXT,
    gene_type TEXT,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Alternative protein names (short names, submitted names, etc.)
CREATE TABLE IF NOT EXISTS protein_names (
    name_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    name_type TEXT,
    name_value TEXT,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Protein comments (function, catalytic activity, subcellular location, etc.)
CREATE TABLE IF NOT EXISTS protein_comments (
    comment_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    comment_type TEXT,
    comment_text TEXT,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Cross-references to other databases (PDB, SMART, etc.)
CREATE TABLE IF NOT EXISTS protein_cross_references (
    xref_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    database_name TEXT,
    database_id TEXT,
    properties TEXT,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Keywords from UniProt
CREATE TABLE IF NOT EXISTS protein_keywords (
    keyword_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    keyword_id_code TEXT,
    keyword_name TEXT,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Gene Ontology (GO) terms
CREATE TABLE IF NOT EXISTS protein_go_terms (
    go_term_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    go_id TEXT,
    go_term TEXT,
    go_aspect TEXT,
    evidence_code TEXT,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Protein-protein interactions
CREATE TABLE IF NOT EXISTS protein_interactions (
    interaction_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    interactant_accession TEXT,
    interactant_name TEXT,
    interaction_type TEXT,
    experiments INTEGER,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Publications/References
CREATE TABLE IF NOT EXISTS protein_references (
    reference_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protein_id INTEGER NOT NULL,
    pubmed_id TEXT,
    doi TEXT,
    title TEXT,
    authors TEXT,
    journal TEXT,
    publication_date TEXT,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id) ON DELETE CASCADE
);

-- Create indexes for common queries
CREATE INDEX IF NOT EXISTS idx_protein_accession ON proteins(accession_id);
CREATE INDEX IF NOT EXISTS idx_protein_name ON proteins(protein_name);
CREATE INDEX IF NOT EXISTS idx_protein_gene ON proteins(gene_name);
CREATE INDEX IF NOT EXISTS idx_created_at ON proteins(created_at);
CREATE INDEX IF NOT EXISTS idx_report_protein ON generated_reports(protein_id);
CREATE INDEX IF NOT EXISTS idx_secondary_accession ON protein_secondary_accessions(secondary_accession);
CREATE INDEX IF NOT EXISTS idx_gene_name ON protein_genes(gene_name);
CREATE INDEX IF NOT EXISTS idx_go_id ON protein_go_terms(go_id);
CREATE INDEX IF NOT EXISTS idx_xref_database ON protein_cross_references(database_name);
CREATE INDEX IF NOT EXISTS idx_comment_type ON protein_comments(comment_type);
