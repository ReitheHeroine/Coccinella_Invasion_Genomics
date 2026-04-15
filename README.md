# Population Genomics of Invasion in *Coccinella septempunctata*

A population genomics pipeline that detects signatures of selection in the seven-spot ladybug (*Coccinella septempunctata*) across its native Eurasian range and its invasive North American range. The workflow takes per-chromosome VCFs of whole-genome variants for 88 individuals (China, Western Europe, Eastern Europe, USA), computes population structure and diversity statistics, runs multiple selection scans, maps candidate SNPs to genes via a multi-species ortholog chain, and tests for GO-term overrepresentation to identify biological processes associated with invasion success.

## Biological Context

*C. septempunctata* is a Eurasian coccinellid that was repeatedly introduced to North America beginning in the 1950s and is now broadly established across the continent. Successful biological invasion typically requires rapid adaptation to novel climate, host, and competitor regimes, and the resulting demographic pulses can leave detectable signatures of selection in the genome. This project asks which regions of the *C. septempunctata* genome show evidence of selection between the native and invasive ranges, and whether the implicated genes are enriched for functions plausibly tied to invasion biology (e.g., thermal tolerance, immunity, chemosensation). The framing matters beyond a single species: coccinellids are both biocontrol agents and, as invaders, documented drivers of native ladybug decline, so understanding the genetic basis of their establishment informs both applied pest management and invasion ecology theory.

## Pipeline Overview

```
   raw short reads (upstream, not in this repo)
            │
            ▼
   [ QC, alignment to GCF_907165205.1, variant calling ]
            │
            ▼
   per-chromosome raw VCFs (NC_058189.1 ... NC_058198.1)
            │                                    │
            ▼                                    ▼
   [ chr_filter.sh: vcftools          [ create_metadata.py ]
     MAF 0.025, missingness ≤12.5%,
     QUAL ≥50, mean DP 10-50,
     indels removed, per-sample DP 10-50 ]
            │                            popmap.txt, pop_*.txt
            ▼                                    │
   ┌────────┴────────┐                           │
   ▼                 ▼                           │
[ vcftools ]    [ PLINK BED/BIM/FAM ]            │
   │                 │                           │
   ├─ Tajima's D     └──► [ pcadapt scan ]       │
   ├─ π (pi)                  │                  │
   └─ FST (Weir &             ▼                  │
        Cockerham)      PC-based outliers        │
            │           (q < 0.05 / 0.01)        │
            ▼                 │                  │
     percentile FST           │                  │
     outliers (top 0.5/1/2%)  │                  │
            │                 │                  │
            └────────┬────────┘                  │
                     ▼                           │
         [ cross_method_concordance.R ]          │
                     │                           │
                     ▼                           │
     [ OutFLANK validation (df=2 fit) ]          │
                     │                           │
                     ▼                           │
              candidate SNPs ◄───────────────────┘
                     │
                     ▼
    [ snp_to_gene_mapping.R / annotate_outliers_local.py ]
        nearest gene from GFF (genic vs intergenic + distance)
                     │
                     ▼
        [ BLAST ortholog chain: Csep → Harmonia → Tribolium ]
                     │
                     ▼
       [ gene_to_go_mapping.py ]
         Route A: NCBI gene2go (taxid 41139, streamed)
         Route B: ortholog transfer via Tribolium gene2go
                     │
                     ▼
          [ go_enrichment_ora.R ]
         clusterProfiler ORA with custom background
                     │
                     ▼
       enriched GO terms + candidate gene table
```

Steps through the VCF-based scans are automated by the shell and R driver scripts in `scripts/`. SNP-to-gene mapping and GO enrichment are run as a sequenced set of scripts, not a single workflow manager; there is no Snakefile or Makefile. A separate environmental niche model (`scripts/Future_Climate_Projection`) projects current and future habitat suitability under climate scenarios using WorldClim 2.1 and GBIF occurrences; it is run independently of the genomic pipeline.

### Upstream processing (pre-VCF)

Short-read QC, alignment to the *C. septempunctata* reference assembly (GCF_907165205.1), and joint variant calling were performed upstream of this repository and are not reproduced here; the pipeline consumes per-chromosome raw VCFs as its entry point. The one upstream step that is tracked is site- and genotype-level filtering, implemented in `scripts/chr_filter.sh`. For each of the ten autosomal chromosome VCFs (NC_058189.1 through NC_058198.1), the script runs vcftools with the following filters and writes a bgzipped, tabix-indexed output to `FILTERED_CHR_VCFS/`:

- `--remove-indels` — retain biallelic SNPs only (selection-scan methods downstream assume SNPs)
- `--maf 0.025` — minor-allele-frequency floor, suppresses singletons and very rare sites where allele-frequency-based statistics (FST, pcadapt loadings) are unstable
- `--max-missing 0.875` — drop sites missing in more than 12.5% of samples
- `--minQ 50` — variant QUAL threshold
- `--min-meanDP 10 --max-meanDP 50` — per-site mean depth bounds, excluding low-coverage and collapsed/repetitive regions
- `--minDP 10 --maxDP 50` — per-genotype depth bounds applied with the same range

All downstream scripts expect the filtered, indexed per-chromosome VCFs produced by this step.

## Getting Started

### Requirements

External tools (install via conda or your package manager):

```
vcftools, plink (1.9), BLAST+ (makeblastdb, blastp), bcftools
```

Python (3.10+):

```
pip install biopython pandas openpyxl requests
```

R (4.3+):

```r
install.packages(c("pcadapt", "qvalue", "OutFLANK", "tidyverse",
                   "data.table", "ggplot2", "clusterProfiler",
                   "raster", "dismo", "spThin", "geosphere"))
```

### Configuration

Before running, the pipeline expects:

1. Raw per-chromosome VCFs in `data/VARIANTS_BY_CHR/` (one file per autosome, NC_058189.1 through NC_058198.1, gzipped; filtered VCFs produced by `chr_filter.sh` land in `data/VARIANTS_BY_CHR/FILTERED_CHR_VCFS/` and are indexed by tabix).
2. The *C. septempunctata* reference genome (NCBI assembly GCF_907165205.1) and its GFF3 in `data/coccinella_septempunctata_/`.
3. Reference proteomes and GFFs for the ortholog chain in `data/harmonia_axyridis_ncbi_dataset/`, `data/tribolium_castaneum_ncbi_dataset/`, and `data/leptinotarsa_decemlineata_ncbi_dataset/` (downloaded from NCBI Datasets).
4. NCBI `gene2go.gz` staged locally for streaming lookups (Route A of GO mapping).
5. The sample metadata spreadsheet `metadata/C7_samples_locations_updated.xlsx`.

None of these data files are tracked. Paths are set at the top of each driver script; update them to match your filesystem before running.

### Usage

```bash
# 1. Build population files from the metadata sheet
python metadata/create_metadata.py

# 2. Filter per-chromosome raw VCFs (site- and genotype-level)
bash scripts/chr_filter.sh

# 3. Per-chromosome diversity and pairwise FST
bash scripts/diversity.sh
bash scripts/fst.sh

# 4. Selection scans
Rscript scripts/pcadapt_selection_scan.R
Rscript scripts/unified_selection_scan.R
Rscript scripts/outflank_diagnostic.R
Rscript scripts/cross_method_concordance.R

# 5. SNP-to-gene mapping and functional annotation
Rscript scripts/snp_to_gene_mapping.R
python scripts/annotate_outliers_local.py
python scripts/gene_to_go_mapping.py

# 6. GO overrepresentation analysis
Rscript scripts/go_enrichment_ora.R

# 7. (Independent) Environmental niche modeling
Rscript scripts/Future_Climate_Projection
```

## Project Structure

```
.
├── README.md
├── metadata/
│   ├── create_metadata.py               # xlsx -> popmap + pop_*.txt, regex geocoding
│   ├── samples.tsv                      # sample, population, lat/lon, continent
│   ├── popmap.txt                       # sample -> population
│   ├── pop_{Asia,CHI,EEU,WEU,Europe,USA,North_America}.txt
│   └── C7_samples_locations*.xlsx       # source metadata (collection sites)
├── scripts/
│   ├── chr_filter.sh                    # vcftools: per-chr SNP/site/genotype filtering
│   ├── diversity.sh                     # vcftools: Tajima's D, windowed and site pi
│   ├── fst.sh                           # vcftools: Weir-Cockerham FST, 5 kb windows
│   ├── pcadapt_selection_scan.R         # PCA-based outliers, q-value FDR
│   ├── pcadapt_step2_3_checkin.R        # scree / K selection / missingness QC
│   ├── unified_selection_scan.R         # percentile FST outliers, 6 pairwise comps
│   ├── per_snp_selection_scan.R         # SNP-level selection metrics
│   ├── outflank_diagnostic.R            # OutFLANK goodness-of-fit check
│   ├── calculate_outflank_concordance.R # OutFLANK vs percentile FST overlap
│   ├── cross_method_concordance.R       # pcadapt vs FST outlier agreement
│   ├── merge_selection_metrics.R        # join per-SNP stats across methods
│   ├── snp_to_gene_mapping.R            # nearest-gene assignment from GFF
│   ├── annotate_outliers_local.py       # local BLAST ortholog chain + GO transfer
│   ├── annotate_outliers_with_go.py     # Tribolium-based GO transfer for outliers
│   ├── enrich_from_gff.py               # offline: pull product/biotype from GFF
│   ├── enrich_gene_annotations_local.py # NCBI Gene lookup with BLAST fallback
│   ├── gene_to_go_mapping.py            # two-route GO mapping with persistent cache
│   ├── go_enrichment_ora.R              # clusterProfiler ORA, custom background
│   ├── go_enrichment_analysis.R         # alternative enrichment formulation
│   ├── create_diversity_tables.R        # population-level diversity summary tables
│   ├── create_population_pi_table.R
│   ├── create_population_tajima_table.R
│   ├── summarize_fst.R
│   ├── plot_diversity.R                 # pi and Tajima's D visualizations
│   ├── plot_fst.r                       # FST Manhattan and summary plots
│   ├── create_metadata.py               # duplicate of metadata/create_metadata.py
│   └── Future_Climate_Projection        # MaxEnt ENM: WorldClim 2.1 + GBIF + spThin
├── data/                                # (not tracked) VCFs, reference genomes, BLAST DBs
├── results/                             # (not tracked) regenerable outputs
├── project_notes/                       # (not tracked) internal analysis logs
└── METHODS.md                           # (not tracked) working methods draft
```

## Key Design Decisions

**Percentile FST as the primary scan, OutFLANK reserved for validation.** OutFLANK fits a neutral chi-squared distribution to the FST spectrum and flags outliers exceeding its tail, which is well-suited for low-to-moderate differentiation but becomes overly conservative when many loci have genuinely high FST (as in among-invasive comparisons with strong founder effects). The pipeline therefore uses uniform top-percentile thresholds (0.5%, 1%, 2%) across all six pairwise comparisons for the primary scan, then uses OutFLANK goodness-of-fit diagnostics and df=2 fits as an orthogonal validation layer rather than as the outlier caller.

**PCA-based selection detection (pcadapt) alongside FST.** pcadapt treats allele-frequency variation along principal components as the null, so it tests for selection conditional on the structure actually present in the data rather than on a priori population labels. Running pcadapt and FST in parallel and reporting their intersection via `cross_method_concordance.R` guards against inference artifacts tied to how populations are defined.

**Local multi-species ortholog chain for GO transfer.** Direct NCBI gene2go coverage for *C. septempunctata* (taxid 41139) is thin, so the pipeline layers a fallback: Csep -> *Harmonia axyridis* (50% identity BLAST, close coccinellid relative) -> *Tribolium castaneum* (30-40% identity BLAST, well-annotated Coleoptera reference), with Tribolium gene2go carrying most of the transferred annotation. Doing the BLAST against locally staged proteomes rather than remote services keeps runs reproducible and removes network failure modes; a persistent MD5-indexed cache avoids re-streaming the 9.9 GB gene2go file on subsequent runs.

**Pattern-stratified GO enrichment instead of a single foreground set.** pcadapt PC loadings distinguish SNPs driving the native-vs-invasive axis from those driving within-invasive differentiation, so the ORA is run on stratified foreground sets (e.g., invasion-specific, post-invasion, shared) against the annotated background. Collapsing these into a single enrichment run would dilute signals that act on only one axis of structure and would overcount genes that separate populations for demographic rather than adaptive reasons.

## Tools and Libraries

vcftools (Tajima's D, windowed/site pi, Weir-Cockerham FST), PLINK 1.9 (LD pruning, BED/BIM/FAM conversion), BLAST+ (makeblastdb, blastp), bcftools, pcadapt (PCA-based selection scan), qvalue (Storey FDR), OutFLANK (FST neutrality test), clusterProfiler (GO ORA), tidyverse, data.table, ggplot2, Biopython (SeqIO, GFF parsing, BLAST XML parsing), pandas, openpyxl, requests, NCBI Datasets / Entrez (gene2go, gene records), raster, dismo (MaxEnt), spThin, geosphere, BiodiversityR.

## Author

Reina Hastings. GitHub: [@reinahastings](https://github.com/reinahastings). This is a personal portfolio project; no lab affiliation.