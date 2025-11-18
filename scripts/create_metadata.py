#!/usr/bin/env python3
"""
title: 'create_metadata.py'
# project: 'BIOL624 Final Project: Selection Detection in Lady Beetles'
# author: 'Reina Hastings'
# contact: 'reinahastings13@gmail.com'
# date created: 11/12/2025
# last modified: 11/17/2025
# purpose: 'WIP'
# inputs:
# outputs:
# required libraries: pandas 
# notes:
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
    hemi = m.group(2).upper() if m.group(2) else None
    if hemi == 'S':
        val = -abs(val)
    if hemi == 'W':
        val = -abs(val)
    return val

def looks_like_coord_pair(text: str) -> bool:
    """Return True if 'a, b' looks like coordinates (numbers with optional N/S/E/W)."""
    parts = [p.strip() for p in text.split(',') if p.strip()]
    if len(parts) != 2:
        return False
    a_ok = re.search(NUM, parts[0]) is not None
    b_ok = re.search(NUM, parts[1]) is not None
    return a_ok and b_ok

def parse_location_info(raw: str):
    """
    Parse 'City (Neighborhood), COUNTRY (lat, lon)' into:
    city, country, lat, lon.
    Chooses the parenthetical group that actually contains coordinates.
    """
    if pd.isna(raw):
        return {'city': None, 'country': None, 'lat': None, 'lon': None, 'raw': raw}

    s = str(raw)

    # Find all (...) groups; pick the one that looks like coordinates
    paren_groups = re.findall(r'\(([^()]*)\)', s)
    coord_text = None
    for g in paren_groups[::-1]:  # check from the end; coords usually last
        if looks_like_coord_pair(g):
            coord_text = g
            break

    # Remove the chosen coordinate group from the string (if present)
    if coord_text is not None:
        s_wo_coords = re.sub(r'\(' + re.escape(coord_text) + r'\)\s*', '', s)
    else:
        s_wo_coords = s

    # Also remove any remaining non-coordinate parentheses (e.g., '(Latsida)')
    s_head = re.sub(r'\([^()]*\)', '', s_wo_coords).strip().rstrip(',')

    # Now split head into City, Country
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

    # Parse coordinates if we found them
    lat = lon = None
    if coord_text:
        p = [p.strip() for p in coord_text.split(',') if p.strip()]
        if len(p) == 2:
            lat = to_float_with_hemi(p[0])
            lon = to_float_with_hemi(p[1])

    return {'city': city, 'country': country, 'lat': lat, 'lon': lon, 'raw': raw}



def main():
    ap = argparse.ArgumentParser(description='Create samples.tsv and popmap.txt from Excel metadata.')
    ap.add_argument('--excel', required=True, help='Path to Excel (e.g., C7_samples_locations.xlsx)')
    ap.add_argument('--sheet', default=0, help='Sheet name or index (default 0)')
    ap.add_argument('--outdir', default='', help='Output directory (default: current directory)')
    # Fixed column names expected in the Excel:
    ap.add_argument('--sample-col', default='Sample', help='Excel column with sample IDs (default: Sample)')
    ap.add_argument('--locinfo-col', default='Location Info', help='Excel column with "City, Country (lat, lon)" (default: Location Info)')
    ap.add_argument('--continent-col', default='Continent', help='Excel column with continent (default: Continent)')
    # Population labeling rule:
    ap.add_argument('--group-by', default='continent', choices=['city', 'country', 'continent'],
                    help='Choose how to label populations in popmap (default: city)')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(args.excel, sheet_name=args.sheet)

    # hard checks for expected columns
    for col in (args.sample_col, args.locinfo_col, args.continent_col):
        if col not in df.columns:
            raise SystemExit(f"ERROR: expected column '{col}' not found. Available: {list(df.columns)}")

    parsed = df[args.locinfo_col].apply(parse_location_info).apply(pd.Series)

    samples = pd.DataFrame({
        'sample_id': df[args.sample_col].astype(str),
        'country': parsed['country'],
        'city': parsed['city'],
        'latitude': parsed['lat'],
        'longitude': parsed['lon'],
        'continent': df[args.continent_col]
    })

    # population labels
    if args.group_by == 'city':
        pop = samples['city']
    elif args.group_by == 'country':
        pop = samples['country']
    else:
        pop = samples['continent']

    population = (pop.astype(str)
                    .str.strip()
                    .str.replace(' ', '_', regex=False)
                    .str.replace('[(),]', '', regex=True))

    samples.insert(1, 'population', population)

    # write outputs
    samples_path = outdir / 'samples.tsv'
    popmap_path = outdir / 'popmap.txt'
    parse_report_path = outdir / 'parse_report.tsv'
    missing_path = outdir / 'missing_rows.tsv'

    samples[['sample_id','population','country','city','latitude','longitude','continent']] \
        .to_csv(samples_path, sep='\t', index=False)

    with open(popmap_path, 'w', encoding='utf-8') as f:
        for sid, poplbl in samples[['sample_id','population']].itertuples(index=False):
            f.write(f'{sid}\t{poplbl}\n')

    # debugging report (shows original location info and parsed pieces)
    pr = pd.DataFrame({
        'sample_id': samples['sample_id'],
        'location_info': df[args.locinfo_col],
        'parsed_city': samples['city'],
        'parsed_country': samples['country'],
        'parsed_latitude': samples['latitude'],
        'parsed_longitude': samples['longitude'],
    })
    pr.to_csv(parse_report_path, sep='\t', index=False)

    # flag rows with ANY missing core fields
    required = ['sample_id','population','country','city','latitude','longitude','continent']
    missing_mask = samples[required].isna().any(axis=1) | (samples['population'].astype(str) == 'nan')
    missing_rows = samples.loc[missing_mask]
    if len(missing_rows):
        missing_rows.to_csv(missing_path, sep='\t', index=False)
        print(f'Wrote: {missing_path} (n={len(missing_rows)})')
    else:
        print('No rows with missing required fields.')

    print(f'Wrote: {samples_path}')
    print(f'Wrote: {popmap_path}')
    print(f'Wrote: {parse_report_path}')

if __name__ == '__main__':
    main()