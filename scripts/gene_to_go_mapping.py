#!/usr/bin/env python3
################################################################################
# title: gene_to_go_mapping.py
# project: BIOL624 Final Project -- Population genomics of Coccinella septempunctata
# author: Reina Hastings (reinahastings13@gmail.com)
# date created: 2026-04-14
# last modified: 2026-04-23 (ensure_blast_db: stage both -in and -out through
#                            /tmp to work around makeblastdb's inability to
#                            handle spaces in -out paths; handoff §4.1 fix.)
#
# purpose:
#   Task 4.2: Map foreground and background C. septempunctata genes to GO terms
#   via two routes -- (A) direct lookup in NCBI gene2go for taxid 41139 and
#   (B) two-step orthology transfer (Csep -> Harmonia -> Tribolium) using local
#   BLAST. Both routes are merged into a unified per-gene GO table tagged with
#   evidence source. Feeds Task 4.3 (clusterProfiler ORA, pattern-stratified).
#
# inputs:
#   --foreground   results/enrichment/foreground_genes.tsv
#                  (cols include: gene_id, ncbi_gene_id, ...)
#   --background   results/enrichment/background_genes.tsv
#                  (cols include: gene_id, ncbi_gene_id, ...)
#   --gene2go      data/go_annotations/gene2go (NCBI; ~9.9 GB)
#   --tribolium-go data/go_annotations/tribolium_gene2go.tsv
#   --tribolium-gff data/tribolium_castaneum_ncbi_dataset/.../genomic.gff
#   --csep-gff     data/coccinella_septempunctata_/.../genomic.gff
#   --csep-faa     data/coccinella_septempunctata_/.../protein.faa
#   --harmonia-faa data/harmonia_axyridis_ncbi_dataset/.../protein.faa
#   --tribolium-faa data/tribolium_castaneum_ncbi_dataset/.../protein.faa
#
# outputs (all in results/enrichment/):
#   gene2go_foreground.tsv     gene_id, go_id, go_term, go_ontology,
#                              evidence_source, ortholog_id, blast_identity
#   gene2go_background.tsv     same schema, all background genes
#   gene_go_cache.tsv          persistent cross-run cache (same schema +
#                              cache_date, plus a "no_annotation" sentinel row
#                              for genes confirmed to have no GO terms so we
#                              don't re-BLAST them on the next run)
#   gene_go_cache_meta.txt     MD5 hashes of the inputs that built the cache,
#                              plus build date / parameters. Used for
#                              invalidation.
#   go_mapping_summary.txt     coverage by route, ontology, gene set,
#                              cache hit/miss stats
#
# usage:
#   python scripts/gene_to_go_mapping.py \
#     --enrichment-dir results/enrichment \
#     --data-dir data \
#     --threads 8 [--force-rebuild]
#
# notes:
#   - Csep taxid 41139. Foreground/background ncbi_gene_id is the numeric NCBI
#     GeneID; this is what NCBI gene2go is keyed on.
#   - Cache is keyed on gene_id (LOC...) alone. Foreground/background tables
#     are always regenerated from the cache by filtering to the current gene
#     lists, so changing flanking window or foreground filter does NOT
#     invalidate the cache.
#   - Cache invalidation: any change in the MD5 of gene2go, tribolium_gene2go,
#     or the BLAST input proteomes triggers a full rebuild. --force-rebuild
#     bypasses the cache regardless.
#   - Route B is run only on genes lacking Route A annotation AND not present
#     in the cache, since foreground is a strict subset of background
#     (Task 4.1 contract).
#   - BLAST identity thresholds match prior Dec 2025 run: 50% Csep -> Harmonia,
#     30% Harmonia -> Tribolium (lowered from 40% to retain coverage at the
#     larger phylogenetic distance).
#   - BLAST DB persists to data/blast_db/ to make reruns cheap.
################################################################################

import argparse
import hashlib
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from urllib.parse import unquote

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("ERROR: Biopython required. pip install biopython")

# --- Constants -------------------------------------------------------------- #
CSEP_TAXID = "41139"
TRIB_TAXID = "7070"
ID_THRESH_HARMONIA = 50.0   # Csep protein -> Harmonia
ID_THRESH_TRIBOLIUM = 30.0  # Harmonia protein -> Tribolium
EVALUE = 1e-10

CATEGORY_TO_ONTOLOGY = {
    "Process": "BP", "Function": "MF", "Component": "CC",
    "BP": "BP", "MF": "MF", "CC": "CC",
}


# --- GFF parsing ------------------------------------------------------------ #
def parse_gff_attrs(s):
    out = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = unquote(v)
    return out


def build_gene_to_protein(gff_path):
    """Csep / Harmonia / Tribolium GFF -> {ncbi_gene_id (str): protein_id (str)}.

    Picks the first protein-coding CDS per gene (longest-isoform selection
    is overkill for ortholog assignment at the GO-term level).
    """
    mrna_to_gene_numeric = {}
    gene_to_protein = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            attrs = parse_gff_attrs(fields[8])
            if feature == "mRNA":
                mrna_id = attrs.get("ID", "")
                dbxref = attrs.get("Dbxref", "")
                m = re.search(r"GeneID:(\d+)", dbxref)
                if mrna_id and m:
                    mrna_to_gene_numeric[mrna_id] = m.group(1)
            elif feature == "CDS":
                parent = attrs.get("Parent", "")
                protein_id = attrs.get("protein_id", "")
                if not protein_id:
                    continue
                gene_num = mrna_to_gene_numeric.get(parent)
                if gene_num is None:
                    dbxref = attrs.get("Dbxref", "")
                    m = re.search(r"GeneID:(\d+)", dbxref)
                    gene_num = m.group(1) if m else None
                if gene_num and gene_num not in gene_to_protein:
                    gene_to_protein[gene_num] = protein_id
    return gene_to_protein


def load_protein_seqs(faa_path):
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(faa_path), "fasta")}


# --- Route A: NCBI gene2go direct lookup ------------------------------------ #
def load_gene2go_for_taxon(gene2go_path, taxid, target_genes):
    """Stream NCBI gene2go (very large), keep rows for taxid where GeneID is in
    target_genes. Returns {gene_id (str): [(go_id, go_term, ontology, evidence)]}.
    """
    out = defaultdict(list)
    target = set(str(g) for g in target_genes)
    print(f"  Streaming NCBI gene2go for taxid {taxid} (this is the big file)...",
          flush=True)
    t0 = time.time()
    n_lines = 0
    n_kept = 0
    with open(gene2go_path) as f:
        # Skip header
        f.readline()
        # Header columns: #tax_id GeneID GO_ID Evidence Qualifier GO_term PubMed Category
        for line in f:
            n_lines += 1
            if n_lines % 5_000_000 == 0:
                print(f"    {n_lines:,} lines scanned, {n_kept:,} kept "
                      f"({time.time()-t0:.0f}s)", flush=True)
            # Cheap prefix filter before split
            if not line.startswith(taxid + "\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            gene_id = fields[1]
            if gene_id not in target:
                continue
            go_id = fields[2]
            evidence = fields[3]
            go_term = fields[5]
            category = fields[7]
            ont = CATEGORY_TO_ONTOLOGY.get(category, category)
            out[gene_id].append((go_id, go_term, ont, evidence))
            n_kept += 1
    print(f"    done: {n_lines:,} lines, {n_kept:,} kept, "
          f"{len(out):,} genes annotated ({time.time()-t0:.0f}s)", flush=True)
    return out


def load_tribolium_gene2go(path):
    """Pre-filtered Tribolium gene2go (small)."""
    out = defaultdict(list)
    with open(path) as f:
        f.readline()
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            gene_id = fields[1]
            go_id = fields[2]
            evidence = fields[3]
            go_term = fields[5]
            ont = CATEGORY_TO_ONTOLOGY.get(fields[7], fields[7])
            out[gene_id].append((go_id, go_term, ont, evidence))
    return out


# --- BLAST ------------------------------------------------------------------ #
def ensure_blast_db(faa_path, db_path, label):
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    sentinel = Path(str(db_path) + ".phr")
    if sentinel.exists():
        print(f"  BLAST DB already present: {label}", flush=True)
        return
    print(f"  Building BLAST DB for {label}...", flush=True)
    # makeblastdb fails on paths with spaces for both -in AND -out (see handoff
    # Section 4.1 resolved question 2026-04-23). Stage BOTH through /tmp, then
    # move the built DB sidecars (*.phr/.pin/.psq/.pdb/.pog/.pos/.pot/.ptf/.pto)
    # back to the target directory.
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        staged_in  = tdp / Path(faa_path).name
        staged_out = tdp / db_path.name
        shutil.copy(str(faa_path), str(staged_in))
        cmd = ["makeblastdb", "-in", str(staged_in), "-dbtype", "prot",
               "-out", str(staged_out), "-title", label, "-parse_seqids"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"makeblastdb failed: {r.stderr}")
        moved = 0
        for produced in tdp.glob(db_path.name + ".*"):
            shutil.move(str(produced), str(db_path.parent / produced.name))
            moved += 1
        if moved == 0:
            sys.exit(f"makeblastdb produced no output sidecars for {label}")
    print(f"  Built: {db_path} ({moved} sidecars)", flush=True)


def batch_blastp_tabular(query_seqs, db_path, min_identity, threads, label):
    """Run blastp once for all queries; outfmt 6 with selected columns.
    Returns {qid: {'sseqid', 'pident', 'evalue', 'bitscore'}} for best hit
    above identity threshold.
    """
    if not query_seqs:
        return {}
    print(f"  BLAST {label}: {len(query_seqs):,} queries (min_id "
          f"{min_identity}%)...", flush=True)
    t0 = time.time()
    with tempfile.TemporaryDirectory() as td:
        qfile = Path(td) / "q.faa"
        with open(qfile, "w") as fh:
            for qid, seq in query_seqs.items():
                fh.write(f">{qid}\n{seq}\n")
        ofile = Path(td) / "out.tsv"
        cmd = [
            "blastp", "-query", str(qfile), "-db", str(db_path),
            "-out", str(ofile), "-outfmt",
            "6 qseqid sseqid pident length evalue bitscore",
            "-evalue", str(EVALUE),
            "-max_target_seqs", "1", "-max_hsps", "1",
            "-num_threads", str(threads),
        ]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"blastp failed ({label}): {r.stderr}")
        hits = {}
        with open(ofile) as fh:
            for line in fh:
                qid, sseqid, pident, _length, evalue, bitscore = line.rstrip(
                    "\n").split("\t")
                pident_f = float(pident)
                if pident_f < min_identity:
                    continue
                if qid in hits:
                    continue  # first hit is best (single max_target_seqs)
                # Strip "ref|XP_xxx.x|" wrappers
                if "|" in sseqid:
                    parts = [p for p in sseqid.split("|") if p]
                    sseqid = parts[1] if len(parts) > 1 else parts[0]
                hits[qid] = {
                    "sseqid": sseqid, "pident": pident_f,
                    "evalue": float(evalue), "bitscore": float(bitscore),
                }
    print(f"    {len(hits):,} hits >= {min_identity}% in "
          f"{time.time()-t0:.0f}s", flush=True)
    return hits


# --- Cache I/O -------------------------------------------------------------- #
CACHE_COLS = ["gene_id", "go_id", "go_term", "go_ontology",
              "evidence_source", "ortholog_id", "blast_identity",
              "cache_date"]
NO_ANNOT_SENTINEL = "no_annotation"  # written as evidence_source for genes
                                     # confirmed annotation-less, to short-
                                     # circuit re-BLAST on later runs


def md5_of_file(path, block=4 * 1024 * 1024):
    """Streaming MD5 (handles the 9.9 GB gene2go file without loading it)."""
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            buf = fh.read(block)
            if not buf:
                break
            h.update(buf)
    return h.hexdigest()


def compute_input_hashes(paths):
    """paths: dict label -> Path. Returns dict label -> md5 hex."""
    out = {}
    for label, p in paths.items():
        print(f"  hashing {label} ({p.name})...", flush=True)
        out[label] = md5_of_file(p)
    return out


def read_cache_meta(meta_path):
    """Returns dict of key->value, or None if absent/unparseable."""
    if not meta_path.exists():
        return None
    out = {}
    try:
        with open(meta_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                if "\t" not in line:
                    continue
                k, v = line.split("\t", 1)
                out[k] = v
        return out
    except Exception:
        return None


def write_cache_meta(meta_path, hashes, params):
    with open(meta_path, "w") as fh:
        fh.write("# gene_go_cache_meta.txt -- input fingerprints for cache "
                 "invalidation\n")
        fh.write(f"build_date\t{datetime.now().isoformat(timespec='seconds')}\n")
        for k, v in params.items():
            fh.write(f"param_{k}\t{v}\n")
        for label, h in hashes.items():
            fh.write(f"md5_{label}\t{h}\n")


def cache_is_valid(meta, hashes, params):
    """Compare current input hashes & params against stored meta."""
    if meta is None:
        return False, "no cache meta on disk"
    for label, h in hashes.items():
        stored = meta.get(f"md5_{label}")
        if stored != h:
            return False, f"input changed: {label}"
    for k, v in params.items():
        stored = meta.get(f"param_{k}")
        if stored is not None and stored != str(v):
            return False, f"parameter changed: {k} ({stored} -> {v})"
    return True, "ok"


def read_cache(cache_path):
    """Returns {gene_id: list of cache rows (dicts)}.
    Includes rows with evidence_source == 'no_annotation' (sentinel).
    """
    out: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    if not cache_path.exists():
        return out
    with open(cache_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                f += [""] * (len(header) - len(f))
            row = {c: f[idx[c]] for c in CACHE_COLS if c in idx}
            out[row["gene_id"]].append(row)
    return out


def write_cache(cache_path, cache_rows_by_gene):
    """cache_rows_by_gene: {gene_id: [row dicts]}. Writes header + rows."""
    with open(cache_path, "w") as fh:
        fh.write("\t".join(CACHE_COLS) + "\n")
        for g in sorted(cache_rows_by_gene):
            for r in cache_rows_by_gene[g]:
                fh.write("\t".join(str(r.get(c, "")) for c in CACHE_COLS)
                         + "\n")


def cache_rows_from_route(gene_ids, route_a, route_b, ortholog_meta,
                          today_str):
    """Build cache rows for the given gene_ids, including the no_annotation
    sentinel for genes annotated by neither route.
    """
    out = defaultdict(list)
    for g in gene_ids:
        a = route_a.get(g, [])
        b = route_b.get(g, [])
        if not a and not b:
            out[g].append({
                "gene_id": g, "go_id": "", "go_term": "", "go_ontology": "",
                "evidence_source": NO_ANNOT_SENTINEL,
                "ortholog_id": "", "blast_identity": "",
                "cache_date": today_str,
            })
            continue
        seen_a = set()
        for go_id, go_term, ont, _ev in a:
            key = (go_id, ont)
            if key in seen_a:
                continue
            seen_a.add(key)
            out[g].append({
                "gene_id": g, "go_id": go_id, "go_term": go_term,
                "go_ontology": ont, "evidence_source": "direct",
                "ortholog_id": "", "blast_identity": "",
                "cache_date": today_str,
            })
        seen_b = set()
        meta = ortholog_meta.get(g, {})
        for go_id, go_term, ont, _ev in b:
            key = (go_id, ont)
            if key in seen_b:
                continue
            seen_b.add(key)
            out[g].append({
                "gene_id": g, "go_id": go_id, "go_term": go_term,
                "go_ontology": ont,
                "evidence_source": "orthology_Tcas_via_Haxy",
                "ortholog_id": meta.get("ortholog_id", ""),
                "blast_identity": (f"{meta['pident']:.1f}"
                                   if "pident" in meta else ""),
                "cache_date": today_str,
            })
    return out


def split_cache_into_routes(cache_rows):
    """Inverse of cache_rows_from_route: rebuild route_a / route_b /
    ortholog_meta from cache rows so the rest of the pipeline can consume them
    uniformly. Sentinel rows produce empty entries (gene IS in cache, just
    has no terms).
    """
    route_a = defaultdict(list)
    route_b = defaultdict(list)
    ortholog_meta = {}
    cached_genes = set()
    for g, rows in cache_rows.items():
        cached_genes.add(g)
        for r in rows:
            src = r.get("evidence_source", "")
            if src == NO_ANNOT_SENTINEL:
                continue
            tup = (r["go_id"], r["go_term"], r["go_ontology"], "")
            if src == "direct":
                route_a[g].append(tup)
            elif src.startswith("orthology"):
                route_b[g].append(tup)
                if r.get("ortholog_id"):
                    pid = r.get("blast_identity", "")
                    ortholog_meta[g] = {
                        "ortholog_id": r["ortholog_id"],
                        "pident": float(pid) if pid else 0.0,
                    }
    return route_a, route_b, ortholog_meta, cached_genes


# --- Output writing --------------------------------------------------------- #
def write_gene2go_table(path, rows):
    """rows: list of dicts with the standard 7 columns."""
    cols = ["gene_id", "go_id", "go_term", "go_ontology",
            "evidence_source", "ortholog_id", "blast_identity"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def assemble_rows(gene_ids, route_a, route_b, ortholog_meta):
    """Build the long-format gene2go table for a given gene set.

    route_a: {gene_id: [(go_id, go_term, ont, evidence), ...]} -- direct.
    route_b: {gene_id: [(go_id, go_term, ont, evidence), ...]} -- transferred.
    ortholog_meta: {gene_id: {'ortholog_id': str, 'pident': float}} for B-genes.
    Per-gene union: a (go_id, ontology) pair seen in both routes is emitted
    twice (once per source) so evidence provenance is preserved per the
    handoff acceptance criteria.
    """
    rows = []
    for g in gene_ids:
        seen_pairs_a = set()
        for go_id, go_term, ont, _evidence in route_a.get(g, []):
            key = (go_id, ont)
            if key in seen_pairs_a:
                continue
            seen_pairs_a.add(key)
            rows.append({
                "gene_id": g, "go_id": go_id, "go_term": go_term,
                "go_ontology": ont, "evidence_source": "direct",
                "ortholog_id": "", "blast_identity": "",
            })
        seen_pairs_b = set()
        meta = ortholog_meta.get(g, {})
        for go_id, go_term, ont, _evidence in route_b.get(g, []):
            key = (go_id, ont)
            if key in seen_pairs_b:
                continue
            seen_pairs_b.add(key)
            rows.append({
                "gene_id": g, "go_id": go_id, "go_term": go_term,
                "go_ontology": ont,
                "evidence_source": "orthology_Tcas_via_Haxy",
                "ortholog_id": meta.get("ortholog_id", ""),
                "blast_identity": (f"{meta['pident']:.1f}"
                                   if "pident" in meta else ""),
            })
    return rows


def coverage_summary(label, gene_ids, route_a, route_b):
    """Returns dict with counts/percentages."""
    n = len(gene_ids)
    a_only = b_only = both = neither = 0
    by_ont: dict[str, set[str]] = {"BP": set(), "MF": set(), "CC": set()}
    for g in gene_ids:
        in_a = g in route_a and len(route_a[g]) > 0
        in_b = g in route_b and len(route_b[g]) > 0
        if in_a and in_b:
            both += 1
        elif in_a:
            a_only += 1
        elif in_b:
            b_only += 1
        else:
            neither += 1
        for src in (route_a.get(g, []), route_b.get(g, [])):
            for _, _, ont, _ in src:
                if ont in by_ont:
                    by_ont[ont].add(g)
    annotated = a_only + b_only + both
    def pct(x):
        return f"{100.0*x/n:.1f}%" if n else "n/a"
    return {
        "label": label, "n_genes": n,
        "n_annotated": annotated, "pct_annotated": pct(annotated),
        "n_route_a_only": a_only, "n_route_b_only": b_only,
        "n_both": both, "n_neither": neither,
        "n_BP": len(by_ont["BP"]), "n_MF": len(by_ont["MF"]),
        "n_CC": len(by_ont["CC"]),
    }


def write_summary(path, fg_cov, bg_cov, n_unique_terms, args, blast_stats):
    with open(path, "w") as fh:
        fh.write("Task 4.2: Gene-to-GO Mapping Summary\n")
        fh.write("=" * 60 + "\n")
        fh.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write(f"Csep taxid: {CSEP_TAXID}\n")
        fh.write(f"BLAST identity thresholds: "
                 f"Csep->Harmonia {ID_THRESH_HARMONIA}%, "
                 f"Harmonia->Tribolium {ID_THRESH_TRIBOLIUM}%\n")
        fh.write(f"BLAST evalue: {EVALUE}\n\n")
        for cov in (fg_cov, bg_cov):
            fh.write(f"--- {cov['label']} ({cov['n_genes']} genes) ---\n")
            fh.write(f"  annotated: {cov['n_annotated']} "
                     f"({cov['pct_annotated']})\n")
            fh.write(f"  route A only (direct):     {cov['n_route_a_only']}\n")
            fh.write(f"  route B only (orthology):  {cov['n_route_b_only']}\n")
            fh.write(f"  both routes:               {cov['n_both']}\n")
            fh.write(f"  neither:                   {cov['n_neither']}\n")
            fh.write(f"  by ontology -- BP: {cov['n_BP']}, "
                     f"MF: {cov['n_MF']}, CC: {cov['n_CC']}\n\n")
        fh.write(f"Unique GO terms (foreground union): "
                 f"{n_unique_terms['foreground']}\n")
        fh.write(f"Unique GO terms (background union): "
                 f"{n_unique_terms['background']}\n\n")
        fh.write("--- Cache stats ---\n")
        fh.write(f"  cache_hits (genes served from cache):      "
                 f"{blast_stats.get('cache_hits', 0)}\n")
        fh.write(f"  cache_misses (newly annotated this run):   "
                 f"{blast_stats.get('newly_annotated', 0)}\n")
        fh.write("\n--- BLAST stage stats (this run only) ---\n")
        for k, v in blast_stats.items():
            if k in ("cache_hits", "newly_annotated"):
                continue
            fh.write(f"  {k}: {v}\n")
        fh.write("\n")
        # Acceptance flags
        fg_pct = float(fg_cov["pct_annotated"].rstrip("%"))
        if fg_pct < 40:
            fh.write("[FOLLOW-UP NEEDED] Foreground annotation coverage "
                     f"{fg_pct:.1f}% < 40% -- consider InterProScan as next "
                     "step.\n")
        elif fg_pct < 60:
            fh.write(f"NOTE: Foreground annotation coverage {fg_pct:.1f}% "
                     "below 60% target but above 40% red-flag threshold.\n")
        else:
            fh.write(f"OK: Foreground annotation coverage {fg_pct:.1f}% "
                     "meets >=60% target.\n")


# --- Input loading ---------------------------------------------------------- #
def load_gene_table(path):
    """Returns list of (gene_id_str, ncbi_gene_id_str)."""
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gid_idx = header.index("gene_id")
        ncbi_idx = header.index("ncbi_gene_id")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append((f[gid_idx], f[ncbi_idx]))
    return rows


# --- Main ------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description="Task 4.2: Csep gene -> GO mapping (direct + orthology).")
    ap.add_argument("--enrichment-dir", required=True,
                    help="Dir with foreground_genes.tsv / background_genes.tsv")
    ap.add_argument("--data-dir", required=True,
                    help="Project data/ root")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--skip-route-b", action="store_true",
                    help="Skip BLAST orthology transfer (Route A only)")
    ap.add_argument("--force-rebuild", action="store_true",
                    help="Ignore cache and re-annotate all genes from scratch")
    args = ap.parse_args()

    enr = Path(args.enrichment_dir)
    data = Path(args.data_dir)

    fg_path = enr / "foreground_genes.tsv"
    bg_path = enr / "background_genes.tsv"
    gene2go_path = data / "go_annotations" / "gene2go"
    trib_g2g_path = data / "go_annotations" / "tribolium_gene2go.tsv"
    csep_gff = (data / "coccinella_septempunctata_/ncbi_dataset/data/"
                "GCF_907165205.1/genomic.gff")
    csep_faa = (data / "coccinella_septempunctata_/ncbi_dataset/data/"
                "GCF_907165205.1/protein.faa")
    harmonia_faa = (data / "harmonia_axyridis_ncbi_dataset/ncbi_dataset/data/"
                    "GCF_914767665.1/protein.faa")
    tribolium_faa = (data / "tribolium_castaneum_ncbi_dataset/ncbi_dataset/"
                     "data/GCF_031307605.1/protein.faa")
    tribolium_gff = (data / "tribolium_castaneum_ncbi_dataset/ncbi_dataset/"
                     "data/GCF_031307605.1/genomic.gff")
    blast_db_dir = data / "blast_db"

    for p in (fg_path, bg_path, gene2go_path, trib_g2g_path, csep_gff,
              csep_faa, harmonia_faa, tribolium_faa, tribolium_gff):
        if not p.exists():
            sys.exit(f"Missing input: {p}")

    print("=" * 70)
    print("Task 4.2: Gene-to-GO Mapping (Routes A + B, with cache)")
    print("=" * 70, flush=True)

    cache_path = enr / "gene_go_cache.tsv"
    cache_meta_path = enr / "gene_go_cache_meta.txt"

    # --- Load gene sets --- #
    print("\n[1/7] Loading foreground and background gene tables...",
          flush=True)
    fg_rows = load_gene_table(fg_path)
    bg_rows = load_gene_table(bg_path)
    fg_gene_ids = [g for g, _ in fg_rows]
    bg_gene_ids = [g for g, _ in bg_rows]
    fg_to_ncbi = {g: n for g, n in fg_rows}
    bg_to_ncbi = {g: n for g, n in bg_rows}
    print(f"  foreground: {len(fg_rows)} genes, background: "
          f"{len(bg_rows)} genes", flush=True)
    if not set(fg_gene_ids).issubset(set(bg_gene_ids)):
        n_missing = len(set(fg_gene_ids) - set(bg_gene_ids))
        print(f"  WARNING: {n_missing} foreground genes not in background "
              "(contract violation from Task 4.1)", flush=True)

    # --- Cache validation --- #
    print("\n[2/7] Checking cache validity (input MD5s)...", flush=True)
    hash_inputs = {
        "gene2go": gene2go_path,
        "tribolium_gene2go": trib_g2g_path,
        "harmonia_proteins": harmonia_faa,
        "tribolium_proteins": tribolium_faa,
        "csep_proteins": csep_faa,
    }
    current_hashes = compute_input_hashes(hash_inputs)
    current_params = {
        "csep_taxid": CSEP_TAXID,
        "id_thresh_harmonia": ID_THRESH_HARMONIA,
        "id_thresh_tribolium": ID_THRESH_TRIBOLIUM,
        "evalue": EVALUE,
        "skip_route_b": int(args.skip_route_b),
    }
    stored_meta = read_cache_meta(cache_meta_path)
    valid, reason = cache_is_valid(stored_meta, current_hashes, current_params)
    if args.force_rebuild:
        valid, reason = False, "--force-rebuild flag set"
    if valid:
        print(f"  cache VALID ({reason}); loading {cache_path.name}...",
              flush=True)
        cache_rows = read_cache(cache_path)
        print(f"  loaded cache for {len(cache_rows):,} genes", flush=True)
    else:
        print(f"  cache INVALID or missing ({reason}); will rebuild "
              "from scratch", flush=True)
        cache_rows = defaultdict(list)

    _cached_route_a, _cached_route_b, _cached_meta, cached_genes = (
        split_cache_into_routes(cache_rows))

    # Genes still needing annotation = present in (foreground U background) but
    # not in the cache.
    needed_gene_ids = [g for g in bg_gene_ids if g not in cached_genes]
    needed_set = set(needed_gene_ids)
    # Foreground genes not in background (rare, contract violation) also need
    # annotation
    for g in fg_gene_ids:
        if g not in cached_genes and g not in needed_set:
            needed_gene_ids.append(g)
            needed_set.add(g)
    print(f"  cache hits: {len(cached_genes):,}; "
          f"to annotate: {len(needed_gene_ids):,}", flush=True)

    # NCBI ids for the genes we still need to annotate (Route A scope)
    needed_ncbi = set()
    for g in needed_gene_ids:
        n = bg_to_ncbi.get(g) or fg_to_ncbi.get(g)
        if n:
            needed_ncbi.add(n)

    # --- Route A: direct gene2go for Csep (only for genes not in cache) --- #
    new_route_a = defaultdict(list)
    blast_stats = {
        "cache_hits": len(cached_genes),
        "newly_annotated": 0,
    }
    if needed_gene_ids:
        print("\n[3/7] Route A: direct NCBI gene2go lookup (taxid 41139)...",
              flush=True)
        direct_by_ncbi = load_gene2go_for_taxon(gene2go_path, CSEP_TAXID,
                                                needed_ncbi)
        ncbi_to_genes = defaultdict(list)
        for g in needed_gene_ids:
            n = bg_to_ncbi.get(g) or fg_to_ncbi.get(g)
            if n:
                ncbi_to_genes[n].append(g)
        for n, terms in direct_by_ncbi.items():
            for g in ncbi_to_genes.get(n, []):
                new_route_a[g] = terms
        print(f"  Route A annotated {len(new_route_a):,} of "
              f"{len(needed_gene_ids):,} new genes", flush=True)
    else:
        print("\n[3/7] Route A: SKIPPED (all genes served from cache)",
              flush=True)

    # --- Route B: two-step BLAST orthology transfer (new genes only) --- #
    new_route_b = defaultdict(list)
    new_ortholog_meta = {}
    if args.skip_route_b:
        print("\n[4/7] Route B: SKIPPED (--skip-route-b)", flush=True)
    elif not needed_gene_ids:
        print("\n[4/7] Route B: SKIPPED (no new genes to annotate)",
              flush=True)
    else:
        print("\n[4/7] Route B: setting up BLAST databases...", flush=True)
        ensure_blast_db(harmonia_faa, blast_db_dir / "harmonia_axyridis",
                        "Harmonia axyridis")
        ensure_blast_db(tribolium_faa, blast_db_dir / "tribolium_castaneum",
                        "Tribolium castaneum")

        print("\n[5/7] Building gene -> protein maps and loading sequences...",
              flush=True)
        csep_gene2prot = build_gene_to_protein(csep_gff)
        trib_protein2gene = build_gene_to_protein(tribolium_gff)
        # Invert Tribolium so we can map Tribolium protein hit -> gene
        trib_protein2gene = {v: k for k, v in trib_protein2gene.items()}
        # Rebuild via dedicated parse to avoid losing alt isoforms:
        trib_protein2gene = {}
        with open(tribolium_gff) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "CDS":
                    continue
                attrs = parse_gff_attrs(f[8])
                pid = attrs.get("protein_id", "")
                m = re.search(r"GeneID:(\d+)", attrs.get("Dbxref", ""))
                if pid and m:
                    trib_protein2gene[pid] = m.group(1)
        print(f"  Csep gene->protein: {len(csep_gene2prot):,}; "
              f"Trib protein->gene: {len(trib_protein2gene):,}", flush=True)

        csep_proteins = load_protein_seqs(csep_faa)
        print(f"  Csep proteins loaded: {len(csep_proteins):,}", flush=True)

        # Route B candidates: needed genes that did not get Route A coverage.
        route_b_candidates = [g for g in needed_gene_ids
                              if g not in new_route_a]
        query_seqs = {}
        n_no_protein = 0
        for g in route_b_candidates:
            ncbi = bg_to_ncbi.get(g) or fg_to_ncbi.get(g)
            pid = csep_gene2prot.get(ncbi) if ncbi else None
            if pid and pid in csep_proteins:
                query_seqs[g] = csep_proteins[pid]
            else:
                n_no_protein += 1
        print(f"  Route B candidates: {len(route_b_candidates)} genes; "
              f"{len(query_seqs)} have protein, {n_no_protein} do not",
              flush=True)
        blast_stats["route_b_candidate_genes"] = len(route_b_candidates)
        blast_stats["route_b_with_protein_seq"] = len(query_seqs)

        print("\n[6/7] BLAST step 1: Csep -> Harmonia axyridis", flush=True)
        harm_hits = batch_blastp_tabular(
            query_seqs, blast_db_dir / "harmonia_axyridis",
            ID_THRESH_HARMONIA, args.threads, "Csep->Harmonia")
        blast_stats["harmonia_hits"] = len(harm_hits)

        print("\n[6/7] BLAST step 2: Harmonia -> Tribolium castaneum",
              flush=True)
        harmonia_proteins = load_protein_seqs(harmonia_faa)
        harm_query_seqs = {}
        n_haxy_seq_missing = 0
        for g, h in harm_hits.items():
            sseqid = h["sseqid"]
            if sseqid in harmonia_proteins:
                harm_query_seqs[g] = harmonia_proteins[sseqid]
            else:
                n_haxy_seq_missing += 1
        if n_haxy_seq_missing:
            print(f"    {n_haxy_seq_missing} Harmonia hits missing seq "
                  "(skipped)", flush=True)
        trib_hits = batch_blastp_tabular(
            harm_query_seqs, blast_db_dir / "tribolium_castaneum",
            ID_THRESH_TRIBOLIUM, args.threads, "Harmonia->Tribolium")
        blast_stats["tribolium_hits"] = len(trib_hits)

        print("\n[7/7] Inheriting Tribolium GO terms...", flush=True)
        trib_go = load_tribolium_gene2go(trib_g2g_path)
        n_with_terms = 0
        for g, h in trib_hits.items():
            trib_protein = h["sseqid"]
            trib_gene = trib_protein2gene.get(trib_protein)
            if not trib_gene:
                continue
            terms = trib_go.get(trib_gene, [])
            if not terms:
                continue
            new_route_b[g] = terms
            new_ortholog_meta[g] = {
                "ortholog_id": f"Tcas:{trib_gene}",
                "pident": h["pident"],
            }
            n_with_terms += 1
        blast_stats["genes_with_inherited_GO"] = n_with_terms

    blast_stats["newly_annotated"] = len(
        set(new_route_a) | set(new_route_b))

    # --- Update cache (append new genes; rewrite if invalid) --- #
    today_str = datetime.now().strftime("%Y-%m-%d")
    new_cache_rows = cache_rows_from_route(
        needed_gene_ids, new_route_a, new_route_b, new_ortholog_meta,
        today_str)
    # If cache was invalid we threw out cache_rows above; otherwise we merge.
    merged_cache = defaultdict(list)
    for g, rs in cache_rows.items():
        merged_cache[g] = rs
    for g, rs in new_cache_rows.items():
        merged_cache[g] = rs  # overwrites any prior sentinel
    write_cache(cache_path, merged_cache)
    write_cache_meta(cache_meta_path, current_hashes, current_params)
    print(f"\nCache updated: {len(merged_cache):,} genes total "
          f"({blast_stats['cache_hits']:,} hits, "
          f"{len(needed_gene_ids):,} new entries written)", flush=True)

    # --- Build per-gene-set route_a / route_b from full merged cache --- #
    full_a, full_b, full_meta, _ = split_cache_into_routes(merged_cache)
    route_a_fg = {g: full_a[g] for g in fg_gene_ids if g in full_a}
    route_a_bg = {g: full_a[g] for g in bg_gene_ids if g in full_a}
    route_b_fg = {g: full_b[g] for g in fg_gene_ids if g in full_b}
    route_b_bg = {g: full_b[g] for g in bg_gene_ids if g in full_b}

    # --- Assemble & write output tables (filter cache to current sets) --- #
    print("\nAssembling output tables from cache...", flush=True)
    fg_rows_out = assemble_rows(fg_gene_ids, route_a_fg, route_b_fg, full_meta)
    bg_rows_out = assemble_rows(bg_gene_ids, route_a_bg, route_b_bg, full_meta)
    write_gene2go_table(enr / "gene2go_foreground.tsv", fg_rows_out)
    write_gene2go_table(enr / "gene2go_background.tsv", bg_rows_out)
    print(f"  wrote gene2go_foreground.tsv ({len(fg_rows_out)} rows)",
          flush=True)
    print(f"  wrote gene2go_background.tsv ({len(bg_rows_out)} rows)",
          flush=True)

    # --- Coverage summary --- #
    fg_cov = coverage_summary("Foreground", fg_gene_ids, route_a_fg, route_b_fg)
    bg_cov = coverage_summary("Background", bg_gene_ids, route_a_bg, route_b_bg)
    n_unique = {
        "foreground": len({r["go_id"] for r in fg_rows_out}),
        "background": len({r["go_id"] for r in bg_rows_out}),
    }
    write_summary(enr / "go_mapping_summary.txt",
                  fg_cov, bg_cov, n_unique, args, blast_stats)
    print("  wrote go_mapping_summary.txt", flush=True)

    print("\nDone.", flush=True)
    print(f"Foreground coverage: {fg_cov['pct_annotated']} "
          f"({fg_cov['n_annotated']}/{fg_cov['n_genes']})", flush=True)
    print(f"Background coverage: {bg_cov['pct_annotated']} "
          f"({bg_cov['n_annotated']}/{bg_cov['n_genes']})", flush=True)


if __name__ == "__main__":
    main()
