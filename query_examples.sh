#!/bin/bash
# Common queries to explore the UniProt database

echo "=== Database Overview ==="
sqlite3 uniprot.db "
SELECT
    'Proteins' as table_name, COUNT(*) as count FROM proteins
UNION ALL SELECT 'Features', COUNT(*) FROM protein_features
UNION ALL SELECT 'GO Terms', COUNT(*) FROM protein_go_terms
UNION ALL SELECT 'References', COUNT(*) FROM protein_references
UNION ALL SELECT 'Cross-refs', COUNT(*) FROM protein_cross_references
UNION ALL SELECT 'Keywords', COUNT(*) FROM protein_keywords;
"

echo ""
echo "=== All Proteins ==="
sqlite3 uniprot.db -header -column "
SELECT
    accession,
    protein_recommended_name as name,
    organism_scientific_name as organism,
    gene_name,
    sequence_length as length
FROM proteins;
"

echo ""
echo "=== Feature Types Distribution ==="
sqlite3 uniprot.db -header -column "
SELECT
    feature_type,
    COUNT(*) as count
FROM protein_features
GROUP BY feature_type
ORDER BY count DESC;
"

echo ""
echo "=== GO Terms by Aspect ==="
sqlite3 uniprot.db -header -column "
SELECT
    go_aspect,
    COUNT(*) as count
FROM protein_go_terms
GROUP BY go_aspect;
"

echo ""
echo "=== Top Databases in Cross-references ==="
sqlite3 uniprot.db -header -column "
SELECT
    database_name,
    COUNT(*) as count
FROM protein_cross_references
GROUP BY database_name
ORDER BY count DESC
LIMIT 10;
"
