#!/usr/bin/env python3
"""
Enrich SNP annotations using gene descriptions from GFF file (fully offline).
Extracts product names, gene types, and other info directly from GFF annotations.
"""

import argparse
import os
import sys
from collections import defaultdict
from urllib.parse import unquote


def parse_gff_attributes(attr_string):
    """Parse GFF9 attribute string into dictionary."""
    attrs = {}
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attrs[key] = unquote(value)  # Handle URL encoding like %2C
    return attrs

def load_gene_info_from_gff(gff_file):
    """Load gene descriptions from GFF file."""
    print(f'Loading gene information from {gff_file}...')

    gene_info = {}
    gene_products = defaultdict(list)

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attrs = parse_gff_attributes(fields[8])

            # Get gene-level info
            if feature_type == 'gene':
                gene_id = attrs.get('ID', '')
                if gene_id:
                    dbxref = attrs.get('Dbxref', '')
                    ncbi_geneid = dbxref.replace('GeneID:', '') if 'GeneID:' in dbxref else ''
                    gene_info[gene_id] = {
                        'gene_id': gene_id,
                        'name': attrs.get('Name', gene_id.replace('gene-', '')),
                        'gene_biotype': attrs.get('gene_biotype', 'unknown'),
                        'ncbi_geneid': ncbi_geneid,
                        'description': attrs.get('description', ''),
                        'products': []
                    }

            # Get product descriptions from mRNA/CDS features
            elif feature_type in ['mRNA', 'CDS']:
                parent = attrs.get('Parent', '')
                product = attrs.get('product', '')

                # Find the gene parent
                if parent.startswith('gene-'):
                    gene_parent = parent
                elif parent.startswith('rna-'):
                    # Need to look up the gene from mRNA
                    continue
                else:
                    gene_parent = parent

                if gene_parent and product:
                    gene_products[gene_parent].append(product)

    # Merge product info into gene_info
    for gene_id, products in gene_products.items():
        if gene_id in gene_info:
            # Use the most common/first product description
            gene_info[gene_id]['products'] = list(set(products))
            if products:
                gene_info[gene_id]['description'] = products[0]

    # Also handle mRNA parents
    with open(gff_file) as f:
        mrna_to_gene = {}
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] == 'mRNA':
                attrs = parse_gff_attributes(fields[8])
                mrna_id = attrs.get('ID', '')
                parent = attrs.get('Parent', '')
                product = attrs.get('product', '')
                if mrna_id and parent:
                    mrna_to_gene[mrna_id] = parent
                if parent in gene_info and product and not gene_info[parent]['description']:
                    gene_info[parent]['description'] = product

    print(f'  ✓ Loaded info for {len(gene_info)} genes')
    genes_with_desc = sum(1 for g in gene_info.values() if g['description'])
    print(f'  ✓ {genes_with_desc} genes have product descriptions')

    return gene_info

def enrich_annotations(input_file, output_file, gene_info):
    """Enrich SNP annotations with gene descriptions."""
    print(f'\nEnriching annotations from {input_file}...')

    # Read input
    with open(input_file) as f:
        header = f.readline().strip().split('\t')
        rows = [line.strip().split('\t') for line in f if line.strip()]

    print(f'  Found {len(rows)} SNPs')

    # Find gene_id column
    gene_id_idx = header.index('Gene_ID') if 'Gene_ID' in header else None
    if gene_id_idx is None:
        print("ERROR: Gene_ID column not found")
        sys.exit(1)

    # Add new columns
    new_columns = ['Gene_Symbol', 'Gene_Type', 'Gene_Description', 'NCBI_GeneID', 'Info_Source']
    output_header = header + new_columns

    enriched_rows = []
    stats = {'found': 0, 'not_found': 0}

    for row in rows:
        gene_id = row[gene_id_idx]

        # Handle multiple genes (intergenic)
        if '/' in gene_id:
            gene_ids = [g.strip() for g in gene_id.split('/')]
            gene_id = gene_ids[0]  # Use first gene

        if gene_id in gene_info:
            info = gene_info[gene_id]
            new_values = [
                info['name'],
                info['gene_biotype'],
                info['description'] if info['description'] else 'No description in GFF',
                info['ncbi_geneid'] if info['ncbi_geneid'] else 'NA',
                'GFF'
            ]
            stats['found'] += 1
        else:
            # Try without 'gene-' prefix
            alt_id = (f"gene-{gene_id}" if not gene_id.startswith('gene-')
                      else gene_id.replace('gene-', ''))
            if alt_id in gene_info:
                info = gene_info[alt_id]
                new_values = [
                    info['name'],
                    info['gene_biotype'],
                    info['description'] if info['description'] else 'No description in GFF',
                    info['ncbi_geneid'] if info['ncbi_geneid'] else 'NA',
                    'GFF'
                ]
                stats['found'] += 1
            else:
                new_values = [
                    gene_id.replace('gene-', ''),
                    'unknown',
                    'Gene not found in GFF',
                    'NA',
                    'Not_found'
                ]
                stats['not_found'] += 1

        enriched_rows.append(row + new_values)

    # Write output
    print(f'\nWriting enriched annotations to {output_file}...')
    with open(output_file, 'w') as f:
        f.write('\t'.join(output_header) + '\n')
        for row in enriched_rows:
            f.write('\t'.join(row) + '\n')

    print(f'  ✓ Wrote {len(enriched_rows)} enriched annotations')
    print('\nSummary:')
    print(f'  Found in GFF: {stats["found"]} ({100*stats["found"]/len(rows):.1f}%)')
    print(f'  Not found: {stats["not_found"]} ({100*stats["not_found"]/len(rows):.1f}%)')

def main():
    parser = argparse.ArgumentParser(
        description='Enrich SNP annotations using GFF file (fully offline)'
    )
    parser.add_argument('--input', required=True, help='Annotated SNP file')
    parser.add_argument('--gff', required=True, help='Reference GFF file')
    parser.add_argument('--output', help='Output file (auto-generated if not specified)')
    parser.add_argument('--outdir', default='.', help='Output directory')

    args = parser.parse_args()

    # Auto-generate output filename
    if args.output is None:
        basename = os.path.basename(args.input)
        basename = basename.replace('_annotated_', '_enriched_')
        args.output = os.path.join(args.outdir, basename)

    print(f'Input: {args.input}')
    print(f'GFF: {args.gff}')
    print(f'Output: {args.output}')

    # Load gene info and enrich
    gene_info = load_gene_info_from_gff(args.gff)
    enrich_annotations(args.input, args.output, gene_info)

    print('\n' + '='*60)
    print('✓ ENRICHMENT COMPLETE (offline mode)')
    print('='*60)

if __name__ == '__main__':
    main()
