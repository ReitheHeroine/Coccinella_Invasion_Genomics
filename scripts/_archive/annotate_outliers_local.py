#!/usr/bin/env python3
################################################################################
# title: annotate_outliers_local.py
# project: BIOL624 Final Project
# author: Reina Hastings
# date created: 2025-12-23
#
# purpose:
#   Complete LOCAL annotation pipeline for selection scan outliers:
#   1. Map SNPs to nearest genes using Coccinella septempunctata GFF
#   2. Extract protein sequences from local protein.faa
#   3. BLAST against Harmonia axyridis (closest relative)
#   4. BLAST Harmonia hits against Tribolium castaneum (for GO terms)
#   5. Transfer GO annotations from Tribolium
#
# workflow:
#   Coccinella SNP → Gene (GFF) → Protein (local .faa) →
#   BLAST → Harmonia → BLAST → Tribolium → GO terms
#
# usage:
#   # Setup BLAST databases (run once)
#   python annotate_outliers_local.py --setup-db
#
#   # Run full annotation pipeline
#   python annotate_outliers_local.py \
#       --outliers ../results/whole_pop_results/combined_outliers/all_outliers_top0.5pct.tsv \
#       --outdir ../results/whole_pop_results/combined_outliers
################################################################################

import argparse
import sys
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

# Check for required packages
try:
    from Bio import SeqIO
    from Bio.Blast import NCBIXML
except ImportError:
    print("ERROR: Biopython not installed.")
    print("Please run: pip install biopython")
    sys.exit(1)

################################################################################
# CONFIGURATION - Paths to local data
################################################################################

SCRIPT_DIR = Path(__file__).parent.absolute()
DATA_DIR = SCRIPT_DIR.parent / "data"

# Use /tmp for BLAST databases (avoids path issues with spaces)
BLAST_DB_DIR = Path("/tmp/blast_db_biol624")

# Coccinella septempunctata (main genome)
COCSEPT_DIR = DATA_DIR / "coccinella_septempunctata_/ncbi_dataset/data/GCF_907165205.1"
COCSEPT_GFF = COCSEPT_DIR / "genomic.gff"
COCSEPT_PROTEIN = COCSEPT_DIR / "protein.faa"

# Beetle databases for BLAST
BEETLE_DATABASES = {
    'harmonia': {
        'name': 'Harmonia axyridis',
        'common': 'harlequin ladybird',
        'protein_faa': DATA_DIR / "harmonia_axyridis_ncbi_dataset/ncbi_dataset/data/GCF_914767665.1/protein.faa",
        'gff': DATA_DIR / "harmonia_axyridis_ncbi_dataset/ncbi_dataset/data/GCF_914767665.1/genomic.gff",
        'db_path': BLAST_DB_DIR / "harmonia_axyridis",
        'min_identity': 50.0,
        'priority': 1
    },
    'tribolium': {
        'name': 'Tribolium castaneum',
        'common': 'red flour beetle',
        'protein_faa': DATA_DIR / "tribolium_castaneum_ncbi_dataset/ncbi_dataset/data/GCF_031307605.1/protein.faa",
        'gff': DATA_DIR / "tribolium_castaneum_ncbi_dataset/ncbi_dataset/data/GCF_031307605.1/genomic.gff",
        'db_path': BLAST_DB_DIR / "tribolium_castaneum",
        'min_identity': 40.0,
        'priority': 2
    },
    'leptinotarsa': {
        'name': 'Leptinotarsa decemlineata',
        'common': 'Colorado potato beetle',
        'protein_faa': DATA_DIR / "leptinotarsa_decemlineata_ncbi_dataset/ncbi_dataset/data/GCF_000500325.1/protein.faa",
        'gff': DATA_DIR / "leptinotarsa_decemlineata_ncbi_dataset/ncbi_dataset/data/GCF_000500325.1/genomic.gff",
        'db_path': BLAST_DB_DIR / "leptinotarsa_decemlineata",
        'min_identity': 40.0,
        'priority': 3
    }
}

# GO annotation file
GO_ANNOTATION_FILE = DATA_DIR / "go_annotations/tribolium_gene2go.tsv"

################################################################################
# BLAST DATABASE SETUP
################################################################################

def check_blast_installed():
    """Check if BLAST+ is installed."""
    try:
        result = subprocess.run(['blastp', '-version'], capture_output=True, text=True)
        if result.returncode == 0:
            version = result.stdout.split('\n')[0]
            print(f"  BLAST+ found: {version}")
            return True
    except FileNotFoundError:
        pass

    print("  BLAST+ not found!")
    print("\n  To install BLAST+:")
    print("    macOS:   brew install blast")
    print("    Ubuntu:  sudo apt-get install ncbi-blast+")
    print("    conda:   conda install -c bioconda blast")
    return False

def setup_blast_databases():
    """Create BLAST databases from protein FASTA files."""
    import shutil

    print("\n" + "="*80)
    print("SETTING UP LOCAL BLAST DATABASES")
    print("="*80 + "\n")

    if not check_blast_installed():
        sys.exit(1)

    BLAST_DB_DIR.mkdir(parents=True, exist_ok=True)
    print(f"BLAST database directory: {BLAST_DB_DIR}\n")

    for key, info in BEETLE_DATABASES.items():
        print(f"Creating database for {info['name']}...")

        faa_path = info['protein_faa']
        db_path = info['db_path']

        if not faa_path.exists():
            print(f"  Protein file not found: {faa_path}")
            continue

        # Copy file to /tmp to avoid path issues with spaces
        tmp_faa = BLAST_DB_DIR / f"{key}_protein.faa"
        print(f"  Copying to {tmp_faa}...")
        shutil.copy(str(faa_path), str(tmp_faa))

        # Count sequences
        seq_count = sum(1 for _ in SeqIO.parse(str(tmp_faa), "fasta"))
        print(f"  Input: {faa_path.name} ({seq_count:,} sequences)")

        # Create BLAST database
        cmd = [
            'makeblastdb',
            '-in', str(tmp_faa),
            '-dbtype', 'prot',
            '-out', str(db_path),
            '-title', info['name'],
            '-parse_seqids'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"  Database created: {db_path}")
        else:
            print(f"  Error: {result.stderr}")

    print("\n" + "="*80)
    print("DATABASE SETUP COMPLETE")
    print("="*80 + "\n")

def check_databases_exist():
    """Check if BLAST databases exist."""
    missing = []
    for key, info in BEETLE_DATABASES.items():
        db_files = list(BLAST_DB_DIR.glob(f"{info['db_path'].name}.*"))
        if not db_files:
            missing.append(info['name'])
    return missing

################################################################################
# GFF PARSING AND GENE MAPPING
################################################################################

def parse_gff_attributes(attr_string):
    """Parse GFF9 attribute string into dictionary."""
    attrs = {}
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attrs[key] = unquote(value)
    return attrs

def load_genes_from_gff(gff_file):
    """Load gene information from GFF file."""
    print(f"Loading genes from {gff_file}...")

    genes_by_chrom = defaultdict(list)
    gene_info = {}

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = parse_gff_attributes(fields[8])

            if feature == 'gene':
                gene_id = attrs.get('ID', '')
                if gene_id:
                    genes_by_chrom[chrom].append({
                        'gene_id': gene_id,
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'name': attrs.get('Name', gene_id.replace('gene-', '')),
                        'biotype': attrs.get('gene_biotype', 'unknown'),
                        'description': attrs.get('description', '')
                    })
                    gene_info[gene_id] = genes_by_chrom[chrom][-1]

    # Sort genes by position
    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(key=lambda x: x['start'])

    total = sum(len(g) for g in genes_by_chrom.values())
    print(f"  Loaded {total} genes across {len(genes_by_chrom)} chromosomes")

    return genes_by_chrom, gene_info

def build_gene_to_protein_map(gff_file):
    """Build mapping from gene IDs to protein IDs using GFF."""
    print(f"Building gene-to-protein mapping...")

    gene_to_protein = {}
    mrna_to_gene = {}

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature = fields[2]
            attrs = parse_gff_attributes(fields[8])

            # Map mRNA to gene
            if feature == 'mRNA':
                mrna_id = attrs.get('ID', '')
                parent = attrs.get('Parent', '')
                if mrna_id and parent:
                    mrna_to_gene[mrna_id] = parent

            # Map CDS to protein (via mRNA parent)
            elif feature == 'CDS':
                parent = attrs.get('Parent', '')
                protein_id = attrs.get('protein_id', '')

                if parent and protein_id:
                    # Find gene from mRNA
                    if parent in mrna_to_gene:
                        gene_id = mrna_to_gene[parent]
                        if gene_id not in gene_to_protein:
                            gene_to_protein[gene_id] = protein_id

    print(f"  Mapped {len(gene_to_protein)} genes to proteins")
    return gene_to_protein

def load_protein_sequences(faa_file):
    """Load protein sequences from FASTA file."""
    print(f"Loading protein sequences from {faa_file}...")

    proteins = {}
    for record in SeqIO.parse(str(faa_file), "fasta"):
        proteins[record.id] = str(record.seq)

    print(f"  Loaded {len(proteins)} protein sequences")
    return proteins

################################################################################
# SNP ANNOTATION
################################################################################

def find_nearest_gene(snp_chrom, snp_pos, genes_by_chrom):
    """Find the nearest gene(s) to a SNP position."""
    if snp_chrom not in genes_by_chrom:
        return ('no_genes', None, None)

    genes = genes_by_chrom[snp_chrom]

    # Check if SNP is within any gene (genic)
    for gene in genes:
        if gene['start'] <= snp_pos <= gene['end']:
            return ('genic', gene, 0)

    # Find nearest upstream and downstream genes
    upstream_gene = None
    downstream_gene = None
    upstream_dist = float('inf')
    downstream_dist = float('inf')

    for gene in genes:
        if gene['end'] < snp_pos:
            dist = snp_pos - gene['end']
            if dist < upstream_dist:
                upstream_dist = dist
                upstream_gene = gene
        elif gene['start'] > snp_pos:
            dist = gene['start'] - snp_pos
            if dist < downstream_dist:
                downstream_dist = dist
                downstream_gene = gene

    if upstream_dist < downstream_dist:
        return ('intergenic', upstream_gene, upstream_dist)
    elif downstream_dist < upstream_dist:
        return ('intergenic', downstream_gene, downstream_dist)
    elif upstream_gene:
        return ('intergenic', upstream_gene, upstream_dist)
    else:
        return ('no_nearby_genes', None, None)

def annotate_snps_with_genes(outlier_file, genes_by_chrom, gene_info):
    """Annotate each SNP with nearest gene."""
    print(f"\nAnnotating SNPs from {outlier_file}...")

    annotated = []
    stats = {'genic': 0, 'intergenic': 0, 'no_genes': 0}

    with open(outlier_file, 'r') as f:
        header = f.readline().strip().split('\t')

        # Find column indices
        try:
            locusname_idx = header.index('LocusName')
            chrom_idx = header.index('Chromosome')
        except ValueError as e:
            print(f"ERROR: Required column not found: {e}")
            sys.exit(1)

        for line in f:
            if not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < len(header):
                continue

            row = dict(zip(header, fields))
            snp_id = row['LocusName']
            snp_chrom = row['Chromosome']

            # Extract position
            try:
                if ':' in snp_id:
                    snp_pos = int(snp_id.split(':')[-1])
                else:
                    snp_pos = int(snp_id.split('_')[-1])
            except:
                continue

            # Find nearest gene
            location_type, gene, distance = find_nearest_gene(snp_chrom, snp_pos, genes_by_chrom)

            if location_type == 'genic':
                row['Location'] = 'Genic'
                row['Gene_ID'] = gene['gene_id']
                row['Gene_Name'] = gene['name']
                row['Gene_Biotype'] = gene['biotype']
                row['Distance'] = 0
                stats['genic'] += 1
            elif location_type == 'intergenic' and gene:
                row['Location'] = 'Intergenic'
                row['Gene_ID'] = gene['gene_id']
                row['Gene_Name'] = gene['name']
                row['Gene_Biotype'] = gene['biotype']
                row['Distance'] = distance
                stats['intergenic'] += 1
            else:
                row['Location'] = 'No_genes'
                row['Gene_ID'] = 'NA'
                row['Gene_Name'] = 'NA'
                row['Gene_Biotype'] = 'NA'
                row['Distance'] = 'NA'
                stats['no_genes'] += 1

            annotated.append(row)

    print(f"  Annotated {len(annotated)} SNPs")
    print(f"    Genic: {stats['genic']} ({100*stats['genic']/len(annotated):.1f}%)")
    print(f"    Intergenic: {stats['intergenic']} ({100*stats['intergenic']/len(annotated):.1f}%)")
    print(f"    No genes: {stats['no_genes']}")

    return annotated

################################################################################
# BLAST FUNCTIONS
################################################################################

def run_blastp(query_seq, db_path, min_identity=40.0, evalue=1e-10):
    """Run BLASTP and return best hit."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
        tmp.write(f">query\n{query_seq}\n")
        query_file = tmp.name

    output_file = tempfile.NamedTemporaryFile(suffix='.xml', delete=False).name

    try:
        cmd = [
            'blastp',
            '-query', query_file,
            '-db', str(db_path),
            '-out', output_file,
            '-outfmt', '5',
            '-evalue', str(evalue),
            '-max_target_seqs', '1',
            '-num_threads', '4'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            return None

        with open(output_file, 'r') as f:
            records = list(NCBIXML.parse(f))

        for record in records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    identity_pct = (hsp.identities / hsp.align_length) * 100

                    if identity_pct >= min_identity:
                        # Extract protein ID from hit (format: ref|XP_12345.1|)
                        hit_id = alignment.hit_id
                        if '|' in hit_id:
                            parts = hit_id.split('|')
                            # Get the accession (middle part, e.g., XP_12345.1)
                            hit_id = parts[1] if len(parts) > 1 else parts[0]

                        return {
                            'hit_id': hit_id,
                            'hit_def': alignment.hit_def,
                            'identity': identity_pct,
                            'evalue': hsp.expect,
                            'score': hsp.score
                        }
        return None

    finally:
        os.unlink(query_file)
        if os.path.exists(output_file):
            os.unlink(output_file)

def batch_blastp(query_seqs, db_path, min_identity=40.0, evalue=1e-10):
    """Run BLASTP on multiple sequences at once for efficiency."""
    if not query_seqs:
        return {}

    # Write all queries to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
        for qid, seq in query_seqs.items():
            tmp.write(f">{qid}\n{seq}\n")
        query_file = tmp.name

    output_file = tempfile.NamedTemporaryFile(suffix='.xml', delete=False).name

    results = {}

    try:
        cmd = [
            'blastp',
            '-query', query_file,
            '-db', str(db_path),
            '-out', output_file,
            '-outfmt', '5',
            '-evalue', str(evalue),
            '-max_target_seqs', '1',
            '-num_threads', '4'
        ]

        subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        with open(output_file, 'r') as f:
            records = list(NCBIXML.parse(f))

        for record in records:
            query_id = record.query.split()[0]

            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    identity_pct = (hsp.identities / hsp.align_length) * 100

                    if identity_pct >= min_identity:
                        # Extract protein ID from hit (format: ref|XP_12345.1|)
                        hit_id = alignment.hit_id
                        if '|' in hit_id:
                            parts = hit_id.split('|')
                            hit_id = parts[1] if len(parts) > 1 else parts[0]

                        results[query_id] = {
                            'hit_id': hit_id,
                            'hit_def': alignment.hit_def,
                            'identity': identity_pct,
                            'evalue': hsp.expect,
                            'score': hsp.score
                        }
                        break
                break

        return results

    except subprocess.TimeoutExpired:
        print("  BLAST timeout - processing partial results")
        return results
    finally:
        os.unlink(query_file)
        if os.path.exists(output_file):
            os.unlink(output_file)

################################################################################
# GO TERM ANNOTATION
################################################################################

def load_tribolium_go_annotations(go_file):
    """Load GO annotations for Tribolium."""
    print(f"Loading GO annotations from {go_file}...")

    gene_go = defaultdict(list)

    with open(go_file, 'r') as f:
        header = f.readline()
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 6:
                gene_id = fields[1]
                go_id = fields[2]
                go_term = fields[5]
                category = fields[7] if len(fields) > 7 else 'Unknown'

                gene_go[gene_id].append({
                    'go_id': go_id,
                    'go_term': go_term,
                    'category': category
                })

    print(f"  Loaded GO terms for {len(gene_go)} genes")
    return gene_go

def extract_gene_id_from_hit(hit_def):
    """Extract gene ID (LOC number) from BLAST hit definition."""
    # Pattern: "protein name LOC123456789 ..." or "[Species]"
    match = re.search(r'LOC(\d+)', hit_def)
    if match:
        return match.group(1)
    return None

def build_tribolium_protein_to_gene(gff_file):
    """Build mapping from Tribolium protein IDs to gene IDs."""
    print(f"Building Tribolium protein-to-gene mapping...")

    protein_to_gene = {}

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            if fields[2] == 'CDS':
                attrs = parse_gff_attributes(fields[8])
                protein_id = attrs.get('protein_id', '')

                # Get gene from Dbxref
                dbxref = attrs.get('Dbxref', '')
                gene_match = re.search(r'GeneID:(\d+)', dbxref)

                if protein_id and gene_match:
                    protein_to_gene[protein_id] = gene_match.group(1)

    print(f"  Mapped {len(protein_to_gene)} proteins to genes")
    return protein_to_gene

################################################################################
# MAIN PIPELINE
################################################################################

def run_annotation_pipeline(outlier_file, outdir):
    """Run the complete annotation pipeline."""

    print("\n" + "="*80)
    print("LOCAL OUTLIER ANNOTATION PIPELINE")
    print("Coccinella septempunctata -> Harmonia axyridis -> Tribolium castaneum -> GO")
    print("="*80 + "\n")

    # Check BLAST databases
    missing = check_databases_exist()
    if missing:
        print(f"ERROR: Missing BLAST databases: {', '.join(missing)}")
        print("Run with --setup-db first")
        sys.exit(1)

    # Step 1: Load reference data
    print("STEP 1: Loading reference data...")
    genes_by_chrom, gene_info = load_genes_from_gff(COCSEPT_GFF)
    gene_to_protein = build_gene_to_protein_map(COCSEPT_GFF)
    proteins = load_protein_sequences(COCSEPT_PROTEIN)
    tribolium_go = load_tribolium_go_annotations(GO_ANNOTATION_FILE)
    tribolium_protein_to_gene = build_tribolium_protein_to_gene(BEETLE_DATABASES['tribolium']['gff'])

    # Step 2: Annotate SNPs with genes
    print("\nSTEP 2: Annotating SNPs with nearest genes...")
    annotated_snps = annotate_snps_with_genes(outlier_file, genes_by_chrom, gene_info)

    # Step 3: Get unique genes for BLAST
    print("\nSTEP 3: Preparing genes for BLAST...")
    unique_genes = set()
    for snp in annotated_snps:
        gene_id = snp.get('Gene_ID', 'NA')
        if gene_id != 'NA' and gene_id in gene_to_protein:
            unique_genes.add(gene_id)

    print(f"  Found {len(unique_genes)} unique genes with protein sequences")

    # Prepare query sequences
    query_seqs = {}
    for gene_id in unique_genes:
        protein_id = gene_to_protein.get(gene_id)
        if protein_id and protein_id in proteins:
            query_seqs[gene_id] = proteins[protein_id]

    print(f"  Prepared {len(query_seqs)} protein sequences for BLAST")

    # Step 4: BLAST against Harmonia
    print("\nSTEP 4: BLASTing against Harmonia axyridis...")
    harmonia_db = BEETLE_DATABASES['harmonia']['db_path']
    harmonia_hits = batch_blastp(query_seqs, harmonia_db, min_identity=50.0)
    print(f"  Found {len(harmonia_hits)} Harmonia hits")

    # Step 5: BLAST Harmonia hits against Tribolium
    print("\nSTEP 5: BLASTing Harmonia hits against Tribolium castaneum...")

    # Load Harmonia protein sequences
    harmonia_proteins = load_protein_sequences(BEETLE_DATABASES['harmonia']['protein_faa'])

    # Prepare Harmonia sequences for second BLAST
    harmonia_query_seqs = {}
    gene_to_harmonia = {}
    for gene_id, hit in harmonia_hits.items():
        hit_id = hit['hit_id']
        if hit_id in harmonia_proteins:
            harmonia_query_seqs[gene_id] = harmonia_proteins[hit_id]
            gene_to_harmonia[gene_id] = hit

    print(f"  Prepared {len(harmonia_query_seqs)} Harmonia sequences for second BLAST")

    tribolium_db = BEETLE_DATABASES['tribolium']['db_path']
    tribolium_hits = batch_blastp(harmonia_query_seqs, tribolium_db, min_identity=40.0)
    print(f"  Found {len(tribolium_hits)} Tribolium hits")

    # Step 6: Map to GO terms
    print("\nSTEP 6: Mapping to GO terms...")
    gene_annotations = {}
    go_mapped = 0

    for gene_id in unique_genes:
        annotation = {
            'harmonia_hit': None,
            'harmonia_identity': 0,
            'tribolium_hit': None,
            'tribolium_identity': 0,
            'tribolium_gene_id': None,
            'go_terms': [],
            'go_ids': []
        }

        # Harmonia hit
        if gene_id in harmonia_hits:
            hit = harmonia_hits[gene_id]
            annotation['harmonia_hit'] = hit['hit_def'][:100]
            annotation['harmonia_identity'] = hit['identity']

        # Tribolium hit and GO terms
        if gene_id in tribolium_hits:
            hit = tribolium_hits[gene_id]
            annotation['tribolium_hit'] = hit['hit_def'][:100]
            annotation['tribolium_identity'] = hit['identity']

            # Get Tribolium gene ID
            protein_id = hit['hit_id']
            if protein_id in tribolium_protein_to_gene:
                trib_gene_id = tribolium_protein_to_gene[protein_id]
                annotation['tribolium_gene_id'] = trib_gene_id

                # Get GO terms
                if trib_gene_id in tribolium_go:
                    go_entries = tribolium_go[trib_gene_id]
                    annotation['go_terms'] = [g['go_term'] for g in go_entries[:5]]
                    annotation['go_ids'] = [g['go_id'] for g in go_entries[:5]]
                    go_mapped += 1

        gene_annotations[gene_id] = annotation

    print(f"  Mapped GO terms for {go_mapped} genes")

    # Step 7: Create output
    print("\nSTEP 7: Writing output files...")

    # Determine output filename
    basename = os.path.basename(outlier_file)
    threshold_match = re.search(r'(top\d+(?:\.\d+)?pct)', basename)
    threshold = threshold_match.group(1) if threshold_match else 'annotated'

    output_file = os.path.join(outdir, f'all_outliers_annotated_{threshold}.tsv')

    # Write annotated SNPs
    with open(output_file, 'w') as f:
        # Header
        base_cols = ['LocusName', 'Chromosome', 'Position', 'FST', 'Method',
                     'Comparison', 'Location', 'Gene_ID', 'Gene_Name', 'Gene_Biotype', 'Distance']
        blast_cols = ['Harmonia_Hit', 'Harmonia_Identity', 'Tribolium_Hit',
                      'Tribolium_Identity', 'Tribolium_GeneID', 'GO_IDs', 'GO_Terms']
        f.write('\t'.join(base_cols + blast_cols) + '\n')

        for snp in annotated_snps:
            gene_id = snp.get('Gene_ID', 'NA')
            annot = gene_annotations.get(gene_id, {})

            row = [
                snp.get('LocusName', 'NA'),
                snp.get('Chromosome', 'NA'),
                snp.get('Position', 'NA'),
                snp.get('FST', 'NA'),
                snp.get('Method', 'NA'),
                snp.get('Comparison', 'NA'),
                snp.get('Location', 'NA'),
                gene_id,
                snp.get('Gene_Name', 'NA'),
                snp.get('Gene_Biotype', 'NA'),
                str(snp.get('Distance', 'NA')),
                annot.get('harmonia_hit', 'NA') or 'NA',
                f"{annot.get('harmonia_identity', 0):.1f}" if annot.get('harmonia_identity') else 'NA',
                annot.get('tribolium_hit', 'NA') or 'NA',
                f"{annot.get('tribolium_identity', 0):.1f}" if annot.get('tribolium_identity') else 'NA',
                annot.get('tribolium_gene_id', 'NA') or 'NA',
                ';'.join(annot.get('go_ids', [])) or 'NA',
                ';'.join(annot.get('go_terms', [])) or 'NA'
            ]
            f.write('\t'.join(row) + '\n')

    print(f"  Wrote {len(annotated_snps)} annotated SNPs to {output_file}")

    # Write gene summary
    gene_summary_file = os.path.join(outdir, f'gene_summary_{threshold}.tsv')
    with open(gene_summary_file, 'w') as f:
        f.write('\t'.join(['Gene_ID', 'Gene_Name', 'SNP_Count', 'Harmonia_Identity',
                          'Tribolium_Identity', 'Tribolium_GeneID', 'GO_IDs', 'GO_Terms']) + '\n')

        # Count SNPs per gene
        gene_snp_counts = defaultdict(int)
        gene_names = {}
        for snp in annotated_snps:
            gene_id = snp.get('Gene_ID', 'NA')
            if gene_id != 'NA':
                gene_snp_counts[gene_id] += 1
                gene_names[gene_id] = snp.get('Gene_Name', 'NA')

        for gene_id, count in sorted(gene_snp_counts.items(), key=lambda x: -x[1]):
            annot = gene_annotations.get(gene_id, {})
            row = [
                gene_id,
                gene_names.get(gene_id, 'NA'),
                str(count),
                f"{annot.get('harmonia_identity', 0):.1f}" if annot.get('harmonia_identity') else 'NA',
                f"{annot.get('tribolium_identity', 0):.1f}" if annot.get('tribolium_identity') else 'NA',
                annot.get('tribolium_gene_id', 'NA') or 'NA',
                ';'.join(annot.get('go_ids', [])) or 'NA',
                ';'.join(annot.get('go_terms', [])) or 'NA'
            ]
            f.write('\t'.join(row) + '\n')

    print(f"  Wrote gene summary to {gene_summary_file}")

    # Print summary
    print("\n" + "="*80)
    print("ANNOTATION SUMMARY")
    print("="*80)
    print(f"Total SNPs annotated: {len(annotated_snps)}")
    print(f"Unique genes: {len(unique_genes)}")
    print(f"Genes with Harmonia orthologs: {len(harmonia_hits)}")
    print(f"Genes with Tribolium orthologs: {len(tribolium_hits)}")
    print(f"Genes with GO annotations: {go_mapped}")
    print("="*80 + "\n")

    return output_file

################################################################################
# MAIN
################################################################################

def main():
    parser = argparse.ArgumentParser(
        description='Complete local annotation pipeline for selection scan outliers',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Setup BLAST databases (run once)
  python annotate_outliers_local.py --setup-db

  # Annotate outliers
  python annotate_outliers_local.py \\
      --outliers ../results/whole_pop_results/combined_outliers/all_outliers_top0.5pct.tsv \\
      --outdir ../results/whole_pop_results/combined_outliers
        '''
    )

    parser.add_argument('--setup-db', action='store_true',
                       help='Create BLAST databases from protein files')
    parser.add_argument('--outliers',
                       help='Outlier SNP file (all_outliers_topXpct.tsv)')
    parser.add_argument('--outdir', default='.',
                       help='Output directory')

    args = parser.parse_args()

    if args.setup_db:
        setup_blast_databases()
        return

    if not args.outliers:
        parser.error("--outliers is required (or use --setup-db)")

    if not os.path.exists(args.outliers):
        print(f"ERROR: File not found: {args.outliers}")
        sys.exit(1)

    # Check BLAST+
    print("Checking BLAST+ installation...")
    if not check_blast_installed():
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    run_annotation_pipeline(args.outliers, args.outdir)

    print("ANNOTATION COMPLETE")

if __name__ == '__main__':
    main()
