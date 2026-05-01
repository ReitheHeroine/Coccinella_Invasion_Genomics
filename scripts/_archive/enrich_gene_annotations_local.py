#!/usr/bin/env python3
################################################################################
# title: enrich_gene_annotations_local.py
# project: BIOL624 Final Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2025-12-23
# last modified: 2025-12-23
#
# purpose:
#   Enrich gene annotations using LOCAL BLAST against beetle protein databases:
#   1. Query NCBI Gene database for Coccinella septempunctata info
#   2. If needed, LOCAL BLAST against downloaded protein databases:
#      - Harmonia axyridis (harlequin ladybird - closest relative)
#      - Tribolium castaneum (red flour beetle - well-annotated)
#      - Leptinotarsa decemlineata (Colorado potato beetle)
#
# inputs:
#   - Annotated SNP file from annotate_snps.py
#   - Local protein FASTA files for target species
#
# outputs:
#   - Enriched annotation file with gene symbols, descriptions, orthologs
#   - Summary statistics including BLAST hit quality
#
# requirements:
#   - BLAST+ (install with: brew install blast)
#   - Biopython: pip install biopython
#
# usage:
#   # First time: create BLAST databases
#   python enrich_gene_annotations_local.py --setup-db
#
#   # Run enrichment
#   python enrich_gene_annotations_local.py \
#       --input ../results/whole_pop_results/combined_outliers/top10_snps_annotated_top0.5pct.tsv \
#       --email your@email.com \
#       --outdir ../results/whole_pop_results/combined_outliers
################################################################################

import argparse
import sys
import os
import re
import subprocess
import tempfile
from collections import defaultdict, Counter
from pathlib import Path

# Check for required packages
try:
    from Bio import Entrez, SeqIO
    from Bio.Blast import NCBIXML
except ImportError:
    print("ERROR: Biopython not installed.")
    print("Please run: pip install biopython --break-system-packages")
    sys.exit(1)

################################################################################
# CONFIGURATION - Paths to local data
################################################################################

# Base data directory (relative to script location)
SCRIPT_DIR = Path(__file__).parent.absolute()
DATA_DIR = SCRIPT_DIR.parent / "data"

# BLAST database directory (use /tmp to avoid path issues with spaces)
BLAST_DB_DIR = Path("/tmp/blast_db")

# Protein database paths
BEETLE_DATABASES = {
    'Harmonia_axyridis': {
        'name': 'Harmonia axyridis',
        'common': 'harlequin ladybird',
        'protein_faa': DATA_DIR / "harmonia_axyridis_ncbi_dataset/ncbi_dataset/data/GCF_914767665.1/protein.faa",
        'db_path': BLAST_DB_DIR / "harmonia_axyridis",
        'min_identity': 60.0,
        'priority': 1  # Closest relative
    },
    'Tribolium_castaneum': {
        'name': 'Tribolium castaneum',
        'common': 'red flour beetle',
        'protein_faa': DATA_DIR / "tribolium_castaneum_ncbi_dataset/ncbi_dataset/data/GCF_031307605.1/protein.faa",
        'db_path': BLAST_DB_DIR / "tribolium_castaneum",
        'min_identity': 50.0,
        'priority': 2
    },
    'Leptinotarsa_decemlineata': {
        'name': 'Leptinotarsa decemlineata',
        'common': 'Colorado potato beetle',
        'protein_faa': DATA_DIR / "leptinotarsa_decemlineata_ncbi_dataset/ncbi_dataset/data/GCF_000500325.1/protein.faa",
        'db_path': BLAST_DB_DIR / "leptinotarsa_decemlineata",
        'min_identity': 50.0,
        'priority': 3
    }
}

# C. septempunctata reference
COCSEPT_GFF = DATA_DIR / "cocsept_ncbi_dataset/GCF_907165205.1/genomic.gff"

################################################################################
# BLAST DATABASE SETUP
################################################################################

def check_blast_installed():
    """Check if BLAST+ is installed."""
    try:
        result = subprocess.run(['blastp', '-version'], capture_output=True, text=True)
        if result.returncode == 0:
            version = result.stdout.split('\n')[0]
            print(f"  ✓ BLAST+ found: {version}")
            return True
    except FileNotFoundError:
        pass

    print("  ✗ BLAST+ not found!")
    print("\n  To install BLAST+:")
    print("    macOS:   brew install blast")
    print("    Ubuntu:  sudo apt-get install ncbi-blast+")
    print("    conda:   conda install -c bioconda blast")
    return False

def setup_blast_databases():
    """Create BLAST databases from protein FASTA files."""
    print("\n" + "="*80)
    print("SETTING UP LOCAL BLAST DATABASES")
    print("="*80 + "\n")

    # Check BLAST+ installation
    print("Checking BLAST+ installation...")
    if not check_blast_installed():
        sys.exit(1)

    # Create blast_db directory
    db_dir = DATA_DIR / "blast_db"
    db_dir.mkdir(exist_ok=True)
    print(f"\nBLAST database directory: {db_dir}\n")

    for species_key, info in BEETLE_DATABASES.items():
        print(f"Creating database for {info['name']}...")

        faa_path = info['protein_faa']
        db_path = info['db_path']

        if not faa_path.exists():
            print(f"  ✗ Protein file not found: {faa_path}")
            continue

        # Count sequences
        seq_count = sum(1 for _ in SeqIO.parse(str(faa_path), "fasta"))
        print(f"  Input: {faa_path.name} ({seq_count:,} sequences)")

        # Create BLAST database
        cmd = [
            'makeblastdb',
            '-in', str(faa_path),
            '-dbtype', 'prot',
            '-out', str(db_path),
            '-title', info['name'],
            '-parse_seqids'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"  ✓ Database created: {db_path}")
        else:
            print(f"  ✗ Error creating database: {result.stderr}")

    print("\n" + "="*80)
    print("DATABASE SETUP COMPLETE")
    print("="*80 + "\n")

def check_databases_exist():
    """Check if BLAST databases exist."""
    missing = []
    for species_key, info in BEETLE_DATABASES.items():
        db_files = list(Path(info['db_path']).parent.glob(f"{info['db_path'].name}.*"))
        if not db_files:
            missing.append(info['name'])
    return missing

################################################################################
# PROTEIN SEQUENCE RETRIEVAL
################################################################################

def get_protein_from_gff(gene_id, gff_file):
    """
    Extract protein accession from GFF file for a gene.
    Returns the protein ID if found.
    """
    clean_id = gene_id.replace('gene-', '')

    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                if clean_id in line:
                    # Look for protein_id in attributes
                    if 'protein_id=' in line:
                        match = re.search(r'protein_id=([^;]+)', line)
                        if match:
                            return match.group(1)
    except Exception as e:
        pass

    return None

def fetch_protein_sequence_ncbi(gene_id, email):
    """
    Retrieve protein sequence from NCBI for a C. septempunctata gene.
    """
    Entrez.email = email
    clean_id = gene_id.replace('gene-', '')

    try:
        # Search for protein entries for this gene
        search_term = f"{clean_id}[Gene Name] AND Coccinella septempunctata[Organism]"

        handle = Entrez.esearch(db="protein", term=search_term, retmax=1)
        protein_record = Entrez.read(handle)
        handle.close()

        if protein_record['IdList']:
            protein_id = protein_record['IdList'][0]

            # Fetch protein sequence in FASTA format
            handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
            protein_seq = handle.read()
            handle.close()

            return protein_seq

        return None

    except Exception as e:
        return None

################################################################################
# LOCAL BLAST
################################################################################

def run_local_blast(query_seq, db_path, species_name, min_identity=50.0):
    """
    Run local BLASTP against a beetle protein database.
    Returns best hit information or None.
    """
    # Create temporary query file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
        tmp.write(query_seq)
        query_file = tmp.name

    # Create temporary output file
    output_file = tempfile.NamedTemporaryFile(suffix='.xml', delete=False).name

    try:
        # Run BLASTP
        cmd = [
            'blastp',
            '-query', query_file,
            '-db', str(db_path),
            '-out', output_file,
            '-outfmt', '5',  # XML format
            '-evalue', '1e-5',
            '-max_target_seqs', '3',
            '-num_threads', '2'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            return None

        # Parse results
        with open(output_file, 'r') as f:
            blast_records = list(NCBIXML.parse(f))

        best_hit = None

        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    identity_pct = (hsp.identities / hsp.align_length) * 100

                    if identity_pct >= min_identity and hsp.expect < 1e-5:
                        # Extract info from hit
                        hit_title = alignment.title
                        gene_symbol = extract_gene_symbol(hit_title)

                        hit_info = {
                            'gene_symbol': gene_symbol,
                            'description': hit_title,
                            'identity': identity_pct,
                            'e_value': hsp.expect,
                            'species': species_name,
                            'align_length': hsp.align_length,
                            'score': hsp.score
                        }

                        if best_hit is None or identity_pct > best_hit['identity']:
                            best_hit = hit_info

                        break  # Only use first HSP

        return best_hit

    finally:
        # Cleanup temp files
        os.unlink(query_file)
        if os.path.exists(output_file):
            os.unlink(output_file)

def extract_gene_symbol(blast_title):
    """Extract gene symbol from BLAST hit title."""
    # Remove species bracket
    if '[' in blast_title:
        desc_part = blast_title.split('[')[0].strip()
    else:
        desc_part = blast_title

    # Remove database identifiers
    if '|' in desc_part:
        parts = desc_part.split('|')
        desc_part = parts[-1].strip()

    # Get first few words as gene symbol
    words = desc_part.split()
    if words:
        first_word = words[0]
        if len(first_word) <= 15 and any(c.isalpha() for c in first_word):
            return first_word
        return ' '.join(words[:4])

    return 'unknown'

def blast_against_local_databases(protein_seq, gene_id):
    """
    BLAST protein sequence against all local beetle databases.
    Returns best hit across all species.
    """
    best_overall = None

    # Sort databases by priority (closest relatives first)
    sorted_dbs = sorted(BEETLE_DATABASES.items(), key=lambda x: x[1]['priority'])

    for species_key, info in sorted_dbs:
        db_path = info['db_path']

        # Check if database exists
        db_files = list(Path(db_path).parent.glob(f"{db_path.name}.*"))
        if not db_files:
            print(f"      → {info['name']}: Database not found, skipping")
            continue

        print(f"      → {info['name']}...", end=' ', flush=True)

        hit = run_local_blast(protein_seq, db_path, info['name'], info['min_identity'])

        if hit:
            print(f"{hit['identity']:.1f}% identity - {hit['gene_symbol']}")

            if best_overall is None or hit['identity'] > best_overall['identity']:
                best_overall = hit
                best_overall['common_name'] = info['common']

            # If excellent hit, can stop early
            if hit['identity'] >= 85:
                break
        else:
            print("No significant hits")

    return best_overall

################################################################################
# NCBI GENE DATABASE SEARCH
################################################################################

def fetch_gene_info_from_ncbi(gene_id, email):
    """Fetch gene information from NCBI Gene database."""
    Entrez.email = email
    clean_id = gene_id.replace('gene-', '')

    try:
        search_term = f"{clean_id}[Gene Name] AND Coccinella septempunctata[Organism]"

        handle = Entrez.esearch(db="gene", term=search_term, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if record['IdList']:
            gene_uid = record['IdList'][0]

            handle = Entrez.efetch(db="gene", id=gene_uid, retmode="xml")
            gene_data = Entrez.read(handle)
            handle.close()

            if gene_data and len(gene_data) > 0:
                gene_record = gene_data[0]

                gene_symbol = clean_id
                if 'Entrezgene_gene' in gene_record:
                    gene_ref = gene_record['Entrezgene_gene'].get('Gene-ref', {})
                    if 'Gene-ref_locus' in gene_ref:
                        gene_symbol = gene_ref['Gene-ref_locus']

                description = 'No description available'
                if 'Entrezgene_summary' in gene_record:
                    description = gene_record['Entrezgene_summary']
                elif 'Entrezgene_gene' in gene_record:
                    gene_ref = gene_record['Entrezgene_gene'].get('Gene-ref', {})
                    if 'Gene-ref_desc' in gene_ref:
                        description = gene_ref['Gene-ref_desc']

                return {
                    'Gene_Symbol': gene_symbol,
                    'Gene_Type': 'protein-coding',
                    'Description': description,
                    'NCBI_GeneID': gene_uid,
                    'Source': 'NCBI',
                    'Ortholog_Species': 'Coccinella_septempunctata',
                    'Ortholog_Identity': '100.0'
                }

        return None

    except Exception as e:
        return None

################################################################################
# MAIN ENRICHMENT WORKFLOW
################################################################################

def enrich_annotations(input_file, output_file, email, outdir):
    """Main enrichment workflow using local BLAST."""

    print(f'\n{"="*80}')
    print('GENE ENRICHMENT PIPELINE (LOCAL BLAST)')
    print('Workflow: NCBI Gene Search → Local BLAST (Coleoptera)')
    print(f'{"="*80}\n')

    # Check for missing databases
    missing = check_databases_exist()
    if missing:
        print(f"WARNING: Missing BLAST databases for: {', '.join(missing)}")
        print("Run with --setup-db to create databases\n")

    print(f'Reading annotations from {input_file}...\n')

    # Read input
    annotations = []
    genes_to_query = set()

    with open(input_file, 'r') as f:
        header = f.readline().strip().split('\t')

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= len(header):
                annotation = dict(zip(header, fields))
                annotations.append(annotation)

                gene_id = annotation.get('Gene_ID', 'NA')
                if gene_id != 'NA' and '/' not in gene_id:
                    genes_to_query.add(gene_id)

    print(f'  Found {len(annotations)} SNPs')
    print(f'  Found {len(genes_to_query)} unique genes to query\n')

    # Enrich each gene
    gene_info_cache = {}
    blast_count = 0

    print(f'{"="*80}')
    print(f'ENRICHING {len(genes_to_query)} GENES')
    print(f'{"="*80}\n')

    for i, gene_id in enumerate(sorted(genes_to_query), 1):
        print(f'[{i}/{len(genes_to_query)}] {gene_id}')

        # STEP 1: Try NCBI Gene database
        print(f'    Step 1: NCBI Gene database...', end=' ', flush=True)
        info = fetch_gene_info_from_ncbi(gene_id, email)

        if info:
            print(f'✓ Found')
            print(f'      Symbol: {info["Gene_Symbol"]}')
            desc_preview = info["Description"][:60] if len(info["Description"]) > 60 else info["Description"]
            print(f'      Description: {desc_preview}...')
            gene_info_cache[gene_id] = info
        else:
            print(f'✗ Not found')

            # STEP 2: Get protein sequence and BLAST locally
            print(f'    Step 2: Retrieving protein sequence...', end=' ', flush=True)
            protein_seq = fetch_protein_sequence_ncbi(gene_id, email)

            if protein_seq:
                seq_len = len(protein_seq.split('\n', 1)[1].replace('\n', ''))
                print(f'✓ Got sequence ({seq_len} aa)')
                print(f'    Step 3: Local BLAST against Coleoptera...')

                blast_info = blast_against_local_databases(protein_seq, gene_id)
                blast_count += 1

                if blast_info:
                    print(f'      ✓ Best hit:')
                    print(f'        Symbol: {blast_info["gene_symbol"]}')
                    print(f'        Species: {blast_info["species"]}')
                    print(f'        Identity: {blast_info["identity"]:.1f}%')

                    gene_info_cache[gene_id] = {
                        'Gene_Symbol': blast_info['gene_symbol'],
                        'Gene_Type': 'protein-coding',
                        'Description': blast_info['description'],
                        'NCBI_GeneID': 'NA',
                        'Source': 'Local_BLAST',
                        'Ortholog_Species': blast_info['species'].replace(' ', '_'),
                        'Ortholog_Identity': f"{blast_info['identity']:.1f}"
                    }
                else:
                    print(f'      ✗ No significant hits found')
                    gene_info_cache[gene_id] = {
                        'Gene_Symbol': gene_id.replace('gene-', ''),
                        'Gene_Type': 'unknown',
                        'Description': 'No orthologs found in Coleoptera',
                        'NCBI_GeneID': 'NA',
                        'Source': 'Not_found',
                        'Ortholog_Species': 'NA',
                        'Ortholog_Identity': '0.0'
                    }
            else:
                print(f'✗ No protein sequence available')
                gene_info_cache[gene_id] = {
                    'Gene_Symbol': gene_id.replace('gene-', ''),
                    'Gene_Type': 'unknown',
                    'Description': 'No protein sequence in database',
                    'NCBI_GeneID': 'NA',
                    'Source': 'No_sequence',
                    'Ortholog_Species': 'NA',
                    'Ortholog_Identity': '0.0'
                }

        print()  # Blank line between genes

    print(f'{"="*80}')
    print(f'Performed local BLAST for {blast_count} genes')
    print(f'{"="*80}\n')

    # Enrich SNP annotations
    print(f'Enriching {len(annotations)} SNP annotations...\n')

    enriched_annotations = []

    for annotation in annotations:
        gene_id = annotation.get('Gene_ID', 'NA')

        if gene_id in gene_info_cache:
            info = gene_info_cache[gene_id]
            enriched = {
                **annotation,
                'Gene_Symbol': info['Gene_Symbol'],
                'Gene_Type': info['Gene_Type'],
                'Gene_Description': info['Description'],
                'NCBI_GeneID': info['NCBI_GeneID'],
                'Info_Source': info['Source'],
                'Ortholog_Species': info['Ortholog_Species'],
                'Ortholog_Identity': info['Ortholog_Identity']
            }
        else:
            enriched = {
                **annotation,
                'Gene_Symbol': annotation.get('Gene_Name', 'NA'),
                'Gene_Type': 'unknown',
                'Gene_Description': 'No information available',
                'NCBI_GeneID': 'NA',
                'Info_Source': 'None',
                'Ortholog_Species': 'NA',
                'Ortholog_Identity': '0.0'
            }

        enriched_annotations.append(enriched)

    # Write output
    print(f'Writing enriched annotations to {output_file}...')

    column_order = [
        'SNP_ID', 'Comparison', 'Method', 'Chromosome', 'Position', 'Location',
        'Gene_ID', 'Gene_Symbol', 'Gene_Type', 'Gene_Description',
        'Distance', 'Strand', 'Context',
        'NCBI_GeneID', 'Info_Source', 'Ortholog_Species', 'Ortholog_Identity'
    ]

    with open(output_file, 'w') as f:
        f.write('\t'.join(column_order) + '\n')

        for annotation in enriched_annotations:
            values = [str(annotation.get(col, 'NA')) for col in column_order]
            f.write('\t'.join(values) + '\n')

    print(f'  ✓ Wrote {len(enriched_annotations)} enriched SNP annotations\n')

    # Print summary
    print_summary(enriched_annotations, gene_info_cache)

def print_summary(enriched_snps, gene_info_cache):
    """Print comprehensive enrichment summary."""

    print(f'{"="*80}')
    print('ENRICHMENT SUMMARY')
    print(f'{"="*80}\n')

    # Information sources
    print('Information Sources:')
    source_counts = Counter(snp['Info_Source'] for snp in enriched_snps)
    for source, count in source_counts.most_common():
        pct = 100 * count / len(enriched_snps)
        print(f'  {source:20s}: {count:3d} ({pct:5.1f}%)')

    # Gene types
    print('\nGene Types:')
    type_counts = Counter(snp['Gene_Type'] for snp in enriched_snps)
    for gene_type, count in type_counts.most_common():
        print(f'  {gene_type:20s}: {count:3d}')

    # Ortholog species (for BLAST hits)
    print('\nOrtholog Sources (Local BLAST hits):')
    blast_snps = [s for s in enriched_snps if s['Info_Source'] == 'Local_BLAST']
    if blast_snps:
        species_counts = Counter(snp['Ortholog_Species'] for snp in blast_snps)
        for species, count in species_counts.most_common():
            print(f'  {species:35s}: {count:3d}')
    else:
        print('  (No BLAST hits)')

    # Genes with multiple SNPs
    print('\nGenes with Multiple SNPs (HIGH PRIORITY):')
    gene_snp_counts = Counter(snp['Gene_Symbol'] for snp in enriched_snps
                              if snp['Gene_Symbol'] != 'NA')
    multi_snp_genes = {g: c for g, c in gene_snp_counts.items() if c > 1}

    if multi_snp_genes:
        for gene, count in sorted(multi_snp_genes.items(), key=lambda x: -x[1])[:10]:
            for snp in enriched_snps:
                if snp['Gene_Symbol'] == gene:
                    desc = snp['Gene_Description'][:40]
                    source = snp['Info_Source']
                    break
            print(f'  {gene:20s}: {count} SNPs - {desc}... [{source}]')
    else:
        print('  (None found)')

    print(f'\nTotal SNPs: {len(enriched_snps)}')
    print(f'Unique genes queried: {len(gene_info_cache)}')
    ncbi_count = sum(1 for g in gene_info_cache.values() if g["Source"] == "NCBI")
    blast_count = sum(1 for g in gene_info_cache.values() if g["Source"] == "Local_BLAST")
    print(f'Successfully annotated: {ncbi_count + blast_count}')
    print(f'  - From NCBI: {ncbi_count}')
    print(f'  - From Local BLAST: {blast_count}')
    print()

################################################################################
# MAIN
################################################################################

def main():
    parser = argparse.ArgumentParser(
        description='Enrich gene annotations using LOCAL BLAST',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Setup BLAST databases (run once)
  python enrich_gene_annotations_local.py --setup-db

  # Run enrichment
  python enrich_gene_annotations_local.py \\
      --input top10_snps_annotated_top0.5pct.tsv \\
      --email your@email.com \\
      --outdir ../results/annotations

Requirements:
  - BLAST+ (brew install blast)
  - Biopython (pip install biopython)
        '''
    )

    parser.add_argument('--setup-db', action='store_true',
                       help='Create BLAST databases from protein FASTA files')
    parser.add_argument('--input', required=False,
                       help='Annotated SNP file from annotate_snps.py')
    parser.add_argument('--email', required=False,
                       help='Your email (required by NCBI)')
    parser.add_argument('--output', required=False,
                       help='Output file path (optional - auto-generated)')
    parser.add_argument('--outdir', required=False, default='.',
                       help='Output directory (default: current directory)')

    args = parser.parse_args()

    # Database setup mode
    if args.setup_db:
        setup_blast_databases()
        return

    # Enrichment mode - require input and email
    if not args.input:
        parser.error("--input is required for enrichment (or use --setup-db)")
    if not args.email:
        parser.error("--email is required for enrichment")

    # Validate email
    if '@' not in args.email:
        print('ERROR: Please provide a valid email address')
        sys.exit(1)

    # Check BLAST+ installation
    print("Checking BLAST+ installation...")
    if not check_blast_installed():
        sys.exit(1)

    # Create output directory
    if args.outdir != '.':
        os.makedirs(args.outdir, exist_ok=True)

    # Auto-generate output filename
    if args.output is None:
        input_basename = os.path.basename(args.input)
        input_basename = input_basename.replace('_annotated_', '_enriched_')

        threshold_match = re.search(r'(top\d+(?:\.\d+)?pct)', input_basename)

        if threshold_match:
            threshold = threshold_match.group(1)
            output_filename = f'top10_snps_enriched_{threshold}.tsv'
        else:
            output_filename = 'top10_snps_enriched.tsv'

        args.output = os.path.join(args.outdir, output_filename)
    else:
        if not os.path.dirname(args.output):
            args.output = os.path.join(args.outdir, args.output)

    print(f'\nOutput file: {args.output}')

    try:
        enrich_annotations(args.input, args.output, args.email, args.outdir)

        print(f'{"="*80}')
        print('✓ ENRICHMENT COMPLETE')
        print(f'{"="*80}\n')

    except KeyboardInterrupt:
        print('\n\nInterrupted by user')
        sys.exit(1)
    except Exception as e:
        print(f'\nERROR: {str(e)}')
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
