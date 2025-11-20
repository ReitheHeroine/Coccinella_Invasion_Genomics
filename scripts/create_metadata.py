#!/usr/bin/env python3
"""
title: create_metadata.py
project: BIOL624 Final Project — Selection Detection in Lady Beetles
author: Reina Hastings
contact: reinahastings13@gmail.com
date created: 2025-11-17
last modified: 2025-11-20

purpose:
    Parse Excel metadata, extract city/country/latitude/longitude, and
    generate standardized population files for downstream analyses.

inputs:
    - Excel metadata file (e.g., C7_samples_locations_updated.xlsx)
    - Required columns: Sample, Location Info, Continent (CHI/WEU/EEU/USA)

outputs:
    - samples.tsv
    - popmap.txt
    - pop_<POP>.txt for each population in the chosen grouping
    - parse_report.tsv
    - missing_rows.tsv
    - create_metadata.log

required libraries:
    - pandas, argparse, pathlib, re
    
usage example:
    python create_metadata.py \
        --excel ../metadata/C7_samples_locations_updated.xlsx  \
        --outdir ../metadata \
        --group-by continent
"""

import re
import argparse
from pathlib import Path
import pandas as pd

NUM = r'[+-]?\d+(?:\.\d+)?'
DEG = r'(?:°|º)?'
HEMI = r'[NnSsEeWw]?'
COORD_TOKEN_RE = re.compile(rf'\s*({NUM})\s*{DEG}\s*({HEMI})\s*$')


def to_float_with_hemi(tok: str):
    tok = tok.strip()
    m = COORD_TOKEN_RE.match(tok)
    if not m:
        m2 = re.search(NUM, tok)
        return float(m2.group(0)) if m2 else None
    val = float(m.group(1))
    hemi = (m.group(2) or "").upper()
    if hemi == 'S':
        val = -abs(val)
    if hemi == 'W':
        val = -abs(val)
    return val


def looks_like_coord_pair(text: str) -> bool:
    parts = [p.strip() for p in text.split(',') if p.strip()]
    if len(parts) != 2:
        return False
    a_ok = re.search(NUM, parts[0]) is not None
    b_ok = re.search(NUM, parts[1]) is not None
    return a_ok and b_ok


def parse_location_info(raw: str):
    if pd.isna(raw):
        return {'city': None, 'country': None, 'lat': None, 'lon': None, 'raw': raw}

    s = str(raw)

    # find all (...) groups, pick last one that looks like coordinates
    paren_groups = re.findall(r'\(([^()]*)\)', s)
    coord_text = None
    for g in paren_groups[::-1]:
        if looks_like_coord_pair(g):
            coord_text = g
            break

    if coord_text:
        s_wo_coords = re.sub(r'\(' + re.escape(coord_text) + r'\)\s*', '', s)
    else:
        s_wo_coords = s

    # remove any remaining non-coordinate parentheses
    s_head = re.sub(r'\([^()]*\)', '', s_wo_coords).strip().rstrip(',')

    city = country = None
    if ',' in s_head:
        parts = [p.strip() for p in s_head.split(',') if p.strip()]
        if len(parts) >= 2:
            city = parts[0]
            country = parts[-1]
        else:
            city = s_head
    else:
        toks = s_head.split()
        if len(toks) >= 2:
            city = ' '.join(toks[:-1])
            country = toks[-1]
        else:
            city = s_head

    lat = lon = None
    if coord_text:
        p = [p.strip() for p in coord_text.split(',') if p.strip()]
        if len(p) == 2:
            lat = to_float_with_hemi(p[0])
            lon = to_float_with_hemi(p[1])

    return {'city': city, 'country': country, 'lat': lat, 'lon': lon, 'raw': raw}


def main():
    ap = argparse.ArgumentParser(description='Create samples.tsv and population files from Excel metadata.')
    ap.add_argument('--excel', required=True, help='Path to Excel metadata file')
    ap.add_argument('--outdir', default='', help='Output directory (default: current working dir)')
    ap.add_argument('--sample-col', default='Sample', help='Column with sample IDs')
    ap.add_argument('--locinfo-col', default='Location Info',
                    help='Column containing "City, Country (lat, lon)"')
    ap.add_argument('--continent-col', default='Continent',
                    help='Column containing continent code (e.g., CHI/WEU/EEU/USA)')
    ap.add_argument('--group-by', default='continent', choices=['city', 'country', 'continent'],
                    help='Label grouping for pop files (default: continent)')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # set up log file
    log_path = outdir / 'create_metadata.log'
    log = open(log_path, 'w', encoding='utf-8')

    def log_print(msg: str):
        print(msg)
        log.write(msg + '\n')

    log_print(f"Running create_metadata.py")
    log_print(f"  Excel file   : {args.excel}")
    log_print(f"  Output dir   : {outdir}")
    log_print(f"  group_by     : {args.group_by}")
    log_print(f"  sample_col   : {args.sample_col}")
    log_print(f"  locinfo_col  : {args.locinfo_col}")
    log_print(f"  continent_col: {args.continent_col}")

    # read first sheet
    try:
        df = pd.read_excel(args.excel, sheet_name=0)
    except Exception as e:
        log_print(f"ERROR: Failed to read Excel file: {e}")
        log.close()
        raise SystemExit(f"ERROR: Failed to read Excel file '{args.excel}': {e}")

    log_print(f"Read {len(df)} rows from Excel.")

    # hard column checks
    for col in (args.sample_col, args.locinfo_col, args.continent_col):
        if col not in df.columns:
            log_print(f"ERROR: expected column '{col}' not found. Available: {list(df.columns)}")
            log.close()
            raise SystemExit(f"ERROR: expected column '{col}' not found in Excel sheet.")

    parsed = df[args.locinfo_col].apply(parse_location_info).apply(pd.Series)

    samples = pd.DataFrame({
        'sample_id': df[args.sample_col].astype(str),
        'country': parsed['country'],
        'city': parsed['city'],
        'latitude': parsed['lat'],
        'longitude': parsed['lon'],
        'continent': df[args.continent_col]
    })

    # choose grouping
    if args.group_by == 'city':
        pop = samples['city']
    elif args.group_by == 'country':
        pop = samples['country']
    else:
        pop = samples['continent']

    population = (
        pop.astype(str)
           .str.strip()
           .str.replace(' ', '_', regex=False)
           .str.replace('[(),]', '', regex=True)
    )

    samples.insert(1, 'population', population)

    # outputs
    samples_path = outdir / 'samples.tsv'
    popmap_path = outdir / 'popmap.txt'
    parse_report_path = outdir / 'parse_report.tsv'
    missing_path = outdir / 'missing_rows.tsv'

    samples[['sample_id', 'population', 'country', 'city', 'latitude', 'longitude', 'continent']] \
        .to_csv(samples_path, sep='\t', index=False)

    with open(popmap_path, 'w', encoding='utf-8') as f:
        for sid, poplbl in samples[['sample_id', 'population']].itertuples(index=False):
            f.write(f"{sid}\t{poplbl}\n")

    # log population summary
    pop_counts = samples['population'].value_counts(dropna=False)
    log_print("Population label summary (population -> n_samples):")
    for lbl, count in pop_counts.items():
        log_print(f"  {lbl!r}: {count}")

    # make population-specific pop files
    created_pop_files = []
    for grp in sorted(samples['population'].unique()):
        if pd.isna(grp):
            continue
        popfile = outdir / f"pop_{grp}.txt"
        subset = samples.loc[samples['population'] == grp, 'sample_id']
        subset.to_csv(popfile, index=False, header=False)
        created_pop_files.append(popfile)
        log_print(f"  Wrote pop file for group '{grp}': {popfile} (n={len(subset)})")

    if not created_pop_files:
        log_print("WARNING: No population-specific pop_*.txt files were created (no non-NA groups?).")

    # parse report
    pr = pd.DataFrame({
        'sample_id': samples['sample_id'],
        'location_info': df[args.locinfo_col],
        'parsed_city': samples['city'],
        'parsed_country': samples['country'],
        'parsed_latitude': samples['latitude'],
        'parsed_longitude': samples['longitude'],
    })
    pr.to_csv(parse_report_path, sep='\t', index=False)

    # missing rows
    required = ['sample_id', 'population', 'country', 'city', 'latitude', 'longitude', 'continent']
    missing_mask = samples[required].isna().any(axis=1)
    missing_rows = samples.loc[missing_mask]
    if len(missing_rows) > 0:
        missing_rows.to_csv(missing_path, sep='\t', index=False)
        log_print(f"Wrote rows with missing required fields: {missing_path} (n={len(missing_rows)})")
    else:
        log_print("No rows with missing required fields.")

    log_print(f"Wrote: {samples_path}")
    log_print(f"Wrote: {popmap_path}")
    log_print(f"Wrote: {parse_report_path}")
    log_print(f"Log written to: {log_path}")
    log_print("Done.")

    log.close()


if __name__ == '__main__':
    main()