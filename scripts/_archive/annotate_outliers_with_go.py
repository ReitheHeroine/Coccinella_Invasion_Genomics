#!/usr/bin/env python3
"""
annotate_outliers_with_go.py

Transfer GO annotations from Tribolium castaneum to Coccinella septempunctata
outlier genes using BLAST-based orthology.

Strategy:
1. For genes with direct Tribolium BLAST hits: use GO terms directly
2. For genes with Harmonia hits: BLAST Harmonia protein against Tribolium,
   then transfer GO terms from best Tribolium hit

Author: Generated for BIOL624 Project
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
from collections import defaultdict
import csv

# Paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_DIR = SCRIPT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
GO_DIR = DATA_DIR / "go_annotations"
RESULTS_DIR = PROJECT_DIR / "results"

# BLAST database location
BLAST_DB_DIR = Path("/tmp/blast_db")

# Files
TRIBOLIUM_GENE2GO = GO_DIR / "tribolium_gene2go.tsv"
TRIBOLIUM_PROTEIN2GENE = GO_DIR / "tribolium_protein2gene.tsv"
HARMONIA_PROTEIN2GENE = GO_DIR / "harmonia_protein2gene.tsv"
TRIBOLIUM_PROTEINS = DATA_DIR / "tribolium_castaneum_ncbi_dataset/ncbi_dataset/data/GCF_031307605.1/protein.faa"
HARMONIA_PROTEINS = DATA_DIR / "harmonia_axyridis_ncbi_dataset/ncbi_dataset/data/GCF_914767665.1/protein.faa"


def load_gene2go(filepath):
    """Load NCBI gene2go file into a dictionary mapping GeneID -> list of GO terms."""
    gene2go = defaultdict(list)
    go_descriptions = {}

    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene_id = row['GeneID']
            go_id = row['GO_ID']
            go_term = row['GO_term']
            category = row['Category']

            gene2go[gene_id].append({
                'go_id': go_id,
                'go_term': go_term,
                'category': category
            })
            go_descriptions[go_id] = {'term': go_term, 'category': category}

    print(f"Loaded GO annotations for {len(gene2go)} Tribolium genes")
    return gene2go, go_descriptions


def load_protein2gene(filepath):
    """Load protein ID to gene ID mapping."""
    mapping = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                protein_id, gene_id = parts
                mapping[protein_id] = gene_id
    print(f"Loaded {len(mapping)} protein-to-gene mappings from {filepath.name}")
    return mapping


def load_protein_sequences(fasta_path):
    """Load protein sequences from FASTA file."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                # Extract protein ID from header
                header = line[1:].strip()
                current_id = header.split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def blast_protein_against_tribolium(protein_id, protein_seq, db_path, evalue=1e-10):
    """BLAST a single protein sequence against Tribolium database."""
    import tempfile
    from Bio.Blast import NCBIXML

    # Write query sequence to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(f">{protein_id}\n{protein_seq}\n")
        query_file = f.name

    # Run BLAST
    output_file = tempfile.mktemp(suffix='.xml')

    cmd = [
        'blastp',
        '-query', query_file,
        '-db', str(db_path),
        '-out', output_file,
        '-outfmt', '5',  # XML format
        '-evalue', str(evalue),
        '-max_target_seqs', '1',
        '-num_threads', '4'
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)

        # Parse results
        with open(output_file, 'r') as f:
            blast_records = list(NCBIXML.parse(f))

        if blast_records and blast_records[0].alignments:
            best_hit = blast_records[0].alignments[0]
            best_hsp = best_hit.hsps[0]

            # Extract protein ID from hit_id (not hit_def)
            # Format is: ref|XP_001807778.2|
            raw_hit_id = best_hit.hit_id
            if raw_hit_id.startswith('ref|'):
                hit_id = raw_hit_id.split('|')[1]
            else:
                hit_id = raw_hit_id.split()[0]

            identity = (best_hsp.identities / best_hsp.align_length) * 100

            return {
                'hit_id': hit_id,
                'identity': identity,
                'evalue': best_hsp.expect
            }
    except Exception as e:
        print(f"  BLAST error for {protein_id}: {e}", file=sys.stderr)
    finally:
        # Clean up temp files
        if os.path.exists(query_file):
            os.remove(query_file)
        if os.path.exists(output_file):
            os.remove(output_file)

    return None


def extract_protein_id_from_description(description):
    """Extract protein ID from BLAST hit description."""
    # Format: ref|XP_045464621.1| description [species]
    if description.startswith('ref|'):
        parts = description.split('|')
        if len(parts) >= 2:
            return parts[1]
    # Try extracting XP_ or NP_ pattern
    import re
    match = re.search(r'([XN]P_\d+\.\d+)', description)
    if match:
        return match.group(1)
    return None


def process_outliers(input_file, output_file, gene2go, go_descriptions,
                     trib_protein2gene, harm_protein2gene, harm_sequences):
    """Process outlier SNPs and annotate with GO terms."""
    from Bio.Blast import NCBIXML

    tribolium_db = BLAST_DB_DIR / "tribolium_castaneum"

    # Check if Tribolium BLAST db exists
    if not Path(str(tribolium_db) + ".phr").exists():
        print("Creating Tribolium BLAST database...")
        cmd = [
            'makeblastdb',
            '-in', str(TRIBOLIUM_PROTEINS),
            '-dbtype', 'prot',
            '-out', str(tribolium_db),
            '-parse_seqids'
        ]
        subprocess.run(cmd, check=True)
        print("Tribolium BLAST database created.")

    # Read input file
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        rows = list(reader)

    print(f"Processing {len(rows)} outlier SNPs...")

    # Add GO columns
    new_fieldnames = fieldnames + ['GO_IDs', 'GO_Terms', 'GO_Categories', 'GO_Source']

    results = []
    genes_processed = set()
    go_gene_list = []  # For enrichment analysis

    for i, row in enumerate(rows):
        gene_symbol = row.get('Gene_Symbol', '')
        description = row.get('Gene_Description', '')
        source = row.get('Info_Source', '')
        ortholog_species = row.get('Ortholog_Species', '')

        go_ids = []
        go_source = 'None'

        # Skip if already processed this gene (for enrichment, we only need unique genes)
        gene_key = row.get('Gene_ID', '')

        # Try to get GO terms
        if ortholog_species == 'Tribolium_castaneum' or 'Tribolium' in description:
            # Direct Tribolium hit
            protein_id = extract_protein_id_from_description(description)
            if protein_id and protein_id in trib_protein2gene:
                gene_id = trib_protein2gene[protein_id]
                if gene_id in gene2go:
                    go_ids = [g['go_id'] for g in gene2go[gene_id]]
                    go_source = 'Tribolium_direct'
                    if gene_key not in genes_processed:
                        go_gene_list.append((gene_key, gene_id, go_ids))
                        genes_processed.add(gene_key)

        elif ortholog_species == 'Harmonia_axyridis' or 'Harmonia' in description:
            # Harmonia hit - need to BLAST against Tribolium
            protein_id = extract_protein_id_from_description(description)

            if protein_id and protein_id in harm_sequences and gene_key not in genes_processed:
                print(f"  [{i+1}/{len(rows)}] BLASTing {protein_id} against Tribolium...")

                blast_result = blast_protein_against_tribolium(
                    protein_id,
                    harm_sequences[protein_id],
                    tribolium_db
                )

                if blast_result:
                    trib_protein = blast_result['hit_id']
                    print(f"    -> Hit: {trib_protein} ({blast_result['identity']:.1f}% identity)")

                    if blast_result['identity'] >= 30:  # Lowered threshold for distant orthologs
                        if trib_protein in trib_protein2gene:
                            gene_id = trib_protein2gene[trib_protein]
                            print(f"    -> Tribolium GeneID: {gene_id}")
                            if gene_id in gene2go:
                                go_ids = [g['go_id'] for g in gene2go[gene_id]]
                                go_source = f"Tribolium_via_Harmonia_{blast_result['identity']:.1f}%"
                                go_gene_list.append((gene_key, gene_id, go_ids))
                                print(f"    -> Found {len(go_ids)} GO terms")

                genes_processed.add(gene_key)

        elif source == 'NCBI' and gene_key not in genes_processed:
            # Coccinella gene - try BLAST against Tribolium
            # We'd need to get the protein sequence from Coccinella
            # For now, mark as needing annotation
            go_source = 'Coccinella_no_ortholog'
            genes_processed.add(gene_key)

        # Format GO annotations for output
        if go_ids:
            go_terms = [go_descriptions[gid]['term'] for gid in go_ids if gid in go_descriptions]
            go_cats = [go_descriptions[gid]['category'] for gid in go_ids if gid in go_descriptions]
            row['GO_IDs'] = '; '.join(go_ids)
            row['GO_Terms'] = '; '.join(go_terms[:5])  # Limit to first 5 terms
            row['GO_Categories'] = '; '.join(set(go_cats))
        else:
            row['GO_IDs'] = ''
            row['GO_Terms'] = ''
            row['GO_Categories'] = ''

        row['GO_Source'] = go_source
        results.append(row)

    # Write output file
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=new_fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(results)

    print(f"\nWrote annotated outliers to {output_file}")

    # Write gene-to-GO mapping for enrichment analysis
    go_mapping_file = output_file.parent / "outlier_gene2go.tsv"
    with open(go_mapping_file, 'w') as f:
        f.write("Gene_ID\tTribolium_GeneID\tGO_IDs\n")
        for gene_key, gene_id, gos in go_gene_list:
            f.write(f"{gene_key}\t{gene_id}\t{';'.join(gos)}\n")

    print(f"Wrote gene-to-GO mapping to {go_mapping_file}")

    # Summary
    n_with_go = sum(1 for r in results if r['GO_IDs'])
    print(f"\nSummary:")
    print(f"  Total SNPs: {len(results)}")
    print(f"  SNPs with GO annotations: {n_with_go} ({100*n_with_go/len(results):.1f}%)")
    print(f"  Unique genes with GO: {len(go_gene_list)}")

    return go_gene_list


def main():
    parser = argparse.ArgumentParser(description='Annotate outlier SNPs with GO terms')
    parser.add_argument('input', help='Input enriched outliers TSV file')
    parser.add_argument('-o', '--output', help='Output file (default: input with _go suffix)')
    args = parser.parse_args()

    input_file = Path(args.input)
    if args.output:
        output_file = Path(args.output)
    else:
        output_file = input_file.parent / (input_file.stem + '_go.tsv')

    # Load data
    print("Loading GO annotations...")
    gene2go, go_descriptions = load_gene2go(TRIBOLIUM_GENE2GO)

    print("\nLoading protein-to-gene mappings...")
    trib_protein2gene = load_protein2gene(TRIBOLIUM_PROTEIN2GENE)
    harm_protein2gene = load_protein2gene(HARMONIA_PROTEIN2GENE)

    print("\nLoading Harmonia protein sequences for BLAST...")
    harm_sequences = load_protein_sequences(HARMONIA_PROTEINS)
    print(f"Loaded {len(harm_sequences)} Harmonia protein sequences")

    # Process outliers
    print("\nProcessing outliers...")
    go_gene_list = process_outliers(
        input_file, output_file, gene2go, go_descriptions,
        trib_protein2gene, harm_protein2gene, harm_sequences
    )

    print("\nDone!")


if __name__ == '__main__':
    main()
