# EDTA-Audit
**Taxonomy-aware QC for EDTA transposable element (TE) annotations — with explicit technical failure detection.**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](#license)
[![Python](https://img.shields.io/badge/Python-%E2%89%A53.8-blue.svg)](#requirements)

**EDTA-Audit** is a lightweight, dependency-free post-processing QC tool for auditing *de novo* TE annotations produced by **EDTA**. It is designed for **large-scale pipelines** (HPC/cloud) where you need to quickly identify:

- **Technically incomplete EDTA outputs** (missing/empty key files or directories)
- **Biologically implausible TE profiles** using **conservative, taxonomy-aware heuristics**
- **Actionable diagnostics** (transparent tags; no black-box scoring)

> **Important**: EDTA-Audit is an **auditing/screening layer**, not a replacement for EDTA, and not a final arbiter of “biological truth”.

---

## What it checks

EDTA-Audit performs two layers of QC:

### 1) Technical completeness (execution integrity)
A sample is considered **TECH FAIL** if any required output is missing/empty:

- `*.EDTA.TEanno.sum`
- `*.EDTA.TEanno.gff3`
- `*.EDTA.TElib.fa`
- `*.EDTA.anno/` (directory; optional with a flag)

When technical completeness fails, **biological QC is skipped** and the record is still reported so you can triage reruns.

### 2) Taxonomy-aware biological sanity checks (conservative heuristics)
Given a profile (`plant`, `animal`, `fungi`), EDTA-Audit parses `TEanno.sum` and checks conservative lower bounds and “missing-major-class” style warnings.

Outputs are categorized as:

- **PASS**: metrics meet conservative expectations
- **SUSPECT**: borderline/atypical; manual review recommended
- **FAIL**: below lower bounds; likely incomplete run / collapse / severe under-annotation
- **TECH FAIL**: missing/empty required outputs (reported as `tech_status=FAIL`, `bio_status=NA`)

---

## Quick start

### Requirements
- Linux or macOS (Windows works if your EDTA outputs are accessible in the filesystem)
- Python **≥ 3.8**
- No external dependencies

### Install
```bash
git clone https://github.com/DiLiu-Lab/EDTA-Audit-Tool.git
cd EDTA-Audit-Tool
```

### Run (most common)
```bash
python3 edta_audit.py -d /path/to/EDTA_runs --profile plant
```

### Run (deeper directory structures)
```bash
python3 edta_audit.py -d /path/to/root -m 8 --profile plant
```

---

## Directory layout assumptions

EDTA-Audit scans for **sample directories** under `--dir`. By default it treats any directory whose **basename** matches:

```
^GC[AF]_\d+\.\d+.*_genomic$
```

Examples:
- `GCA_000001405.28_GRCh38.p13_genomic`
- `GCF_000001635.27_GRCm39_genomic`

---

## Output files

All outputs are written under `--out` (default: `audit_out/`):

- `all.tsv` — all samples
- `pass.tsv` — overall PASS
- `suspect.tsv` — overall SUSPECT
- `fail.tsv` — overall FAIL (includes technical failures)
- `summary.txt` — aggregated counts + most common diagnostic tags

### TSV schema (columns)
The TSV files are tab-delimited and include:

- identifiers: `profile`, `genome_id`, `sample_id`, `sample_dir`
- status: `overall_status`, `tech_status`, `bio_status`
- key metrics (from `TEanno.sum`): total TE %, LTR/TIR/Helitron/LINE/SINE %, bp-masked, etc.
- diagnostics: `tags`, `notes`, `missing`
- provenance: `teanno_sum`

Tip: Load `all.tsv` into R/Python/Excel and filter by `overall_status`, `tech_status`, or `tags`.

---

## Real-world batch example [153 plant genomes](benchmarking/demo.png)

Below is an example run on **153** genomes with `--profile plant`, scanning up to depth 6.

### Command
```bash
chmod +x edta_audit.py
./edta_audit.py -d /path/to/EDTA_root -o . -m 6 --profile plant
```

### Console output (excerpt)
```text
[2026-02-08 08:43:11] Scanning sample dirs: /path/to/EDTA_root (recursive=1, max_depth=6)
[2026-02-08 08:43:11] Profile: plant
[2026-02-08 08:43:11] sample_regex: ^GC[AF]_\d+\.\d+.*_genomic$
[2026-02-08 08:43:11] Found sample dirs: 153
[2026-02-08 08:43:11] Done.
[2026-02-08 08:43:11] Total samples: 153 | PASS: 81 | SUSPECT: 2 | FAIL: 70
```

### Aggregated summary (`summary.txt`)
```text
Counts
  TOTAL	153
  OVERALL_PASS	81
  OVERALL_SUSPECT	2
  OVERALL_FAIL	70
  TECH_PASS	104
  TECH_FAIL	49
  BIO_PASS	81
  BIO_SUSPECT	2
  BIO_FAIL	21
  BIO_NA	49
```

### Quick interpretation
- **Overall PASS**: 81 / 153 (**52.9%**)  
- **TECH_FAIL** (missing/empty required outputs): 49 / 153 (**32.0%**) → `bio_status=NA`  
- Among **technical PASS** samples (104), **BIO_FAIL** = 21 (**20.2% of TECH_PASS**)

Most frequent failure drivers [from `summary.txt`](benchmarking/summary.txt):
```text
FAIL reason tags (TECH_* + HARD_FAIL_*)
  TECH_MISSING_TEanno.sum	49
  TECH_MISSING_TEanno.gff3	49
  TECH_MISSING_EDTA.anno/	49
  TECH_MISSING_TElib.fa	44
  HARD_FAIL_LTR_TOTAL	21
  HARD_FAIL_TOTAL_TE	20
  HARD_FAIL_TIR	20
```

> This example illustrates why EDTA-Audit separates **technical integrity** from **biological plausibility**: many “FAIL” cases are actually **incomplete outputs** that should be rerun, while a smaller subset are **biologically implausible** given a plant profile.


## How metrics are computed (TEanno.sum parsing)

EDTA-Audit parses the **Repeat Classes** table inside `TEanno.sum` and extracts:

- total TE fraction: `TOTAL` (preferred) or `TOTAL_INTERSPERSED` as fallback
- LTR summary: sum of `LTR/*`
- LTR “known”: `LTR/Copia + LTR/Gypsy` (when present)
- LTR “unknown”: `LTR/unknown` (case-insensitive fallback)
- TIR summary: sum of `TIR/*`
- Helitron: `NONTIR/helitron` (case-insensitive fallback)
- LINE summary: sum of `LINE/*`
- SINE summary: sum of `SINE/*`

If a class is not reported in `TEanno.sum`, EDTA-Audit emits explicit `WARN_*_NOT_REPORTED` tags.

---

## Default thresholds (conservative)

These are **defaults** intended to catch extreme under-annotation and technical failures.
They are **not** meant to define ideal TE compositions for all lineages.

### Plant profile (default)
- Total TE %: PASS ≥ **20**, SUSPECT ≥ **10**
- LTR total %: PASS ≥ **5**, SUSPECT ≥ **1**
- TIR total %: PASS ≥ **2**, SUSPECT ≥ **0.5**
- Additional hints/warnings:
  - very low Helitron %
  - LINE/SINE not reported or zero
  - high LTR unknown fraction (relative to total LTR)

### Animal profile (default)
- Total TE %: PASS ≥ **10**, SUSPECT ≥ **3**
- Warnings (not hard fails):
  - LTR total % < **0.1**
  - TIR total % < **0.1**
  - non-LTR (LINE + SINE) < **0.2** or not reported

### Fungal profile (default)
- Total TE %: PASS ≥ **5**, SUSPECT ≥ **0.5**
- Warnings (permissive):
  - LTR total % < **0.1**
  - TIR total % < **0.1**
  - Helitron/LINE/SINE not reported

> **Customizing thresholds**: The thresholds are implemented as simple dataclasses in the script and can be edited for lineage-specific expectations (e.g., certain clades with unusually low TE fractions).

---

## Threshold rationale and representative references

EDTA-Audit uses **conservative, profile-based lower bounds** intended to catch *extreme under-annotation* and *technical incompleteness* at scale. The papers below are **representative background** on TE abundance, lineage-specific patterns, and TE annotation practices. They are **not** used as hard “ground truth” for any single genome.

### Plant background
- Genome size variation and TE contribution in plants. [R1]
- Long-term TE landscape insights in *Arabidopsis thaliana*. [R2]
- Benchmarking TE annotation and streamlined pipelines (EDTA context). [R3]

### Animal background
- Review of TE mechanisms and mammalian genomic applications. [R4]
- Near-complete *Drosophila melanogaster* reference assembly (high-quality baseline for TE discovery). [R5]
- Natural TE abundance variation across *Caenorhabditis elegans* species. [R6]
- Regulatory TE compendium in ENCODE (human). [R7]
- TE diversity and developmental expression in zebrafish. [R8]

### Fungal background
- Population-level TE expression dynamics in a fungal crop pathogen. [R9]
- Kingdom-scale analysis of horizontal transfer shaping fungal mobilomes. [R10]

---


## Interpreting statuses & tags

### Status fields
- `tech_status`:
  - `PASS` if required outputs exist and are non-empty
  - `FAIL` if any required output is missing/empty
- `bio_status`:
  - `PASS` / `SUSPECT` / `FAIL` when `tech_status=PASS`
  - `NA` when `tech_status=FAIL` (biological QC skipped)
- `overall_status`:
  - equals `bio_status` when `tech_status=PASS`
  - equals `FAIL` when `tech_status=FAIL`

### Common tag patterns
- Technical missing:
  - `TECH_MISSING_TEanno.sum`
  - `TECH_MISSING_TEanno.gff3`
  - `TECH_MISSING_TElib.fa`
  - `TECH_MISSING_EDTA.anno/`
- Hard biological fails (profile-dependent):
  - `HARD_FAIL_TOTAL_TE`
  - `HARD_FAIL_LTR_TOTAL`
  - `HARD_FAIL_TIR`
- Borderline:
  - `HARD_SUSPECT_TOTAL_TE`
  - `HARD_SUSPECT_LTR_TOTAL`
  - `HARD_SUSPECT_TIR`
- Informational:
  - `INFO_TOTAL_SOURCE_TOTAL` or `INFO_TOTAL_SOURCE_TOTAL_INTERSPERSED`

Use `summary.txt` to see which tags dominate failures in a batch.

---

## Reproducibility & best practices (academic-friendly)

- Keep a snapshot of:
  - EDTA-Audit version/tag
  - EDTA version and command line
  - reference libraries used by EDTA (if any)
  - assembly metadata (N50, contig count, size)
- Report:
  - number of samples screened
  - proportions of technical failures vs biological fails
  - the profile used (`plant`/`animal`/`fungi`) and any threshold edits
- Treat **SUSPECT** as “review required,” not automatically excluded.

---

## FAQ

**Q: Why is my sample overall FAIL but `bio_status` is NA?**  
A: That’s a technical failure. `tech_status=FAIL` and `bio_status=NA` means required outputs were missing/empty, so biological QC was skipped. Check `missing`, `tags`, and rerun EDTA.

**Q: Does FAIL mean the annotation is impossible?**  
A: No. FAIL means the output is **unlikely to be reliable** under conservative lower bounds, or it was technically incomplete. Edge cases exist; review `tags` and consider lineage-specific biology.

**Q: My organism has unusually low TE content. What should I do?**  
A: Edit thresholds (dataclasses in the script) or use a more permissive profile and document the change in your methods.

---

## Contributing

Issues and PRs are welcome:
- new taxa profiles (e.g., “protist”)
- improved TEanno.sum parsing robustness
- additional diagnostics (e.g., library coverage checks)

Please include:
- a minimal reproducible sample directory layout (or synthetic fixture)
- expected `tags`/statuses

---

## License

MIT — see `LICENSE`.

---

## Author

**Di Liu**  
GitHub: DiLiu-Lab

---

## Acknowledgements

- **EDTA**: Extensive de novo Transposable Element Annotator (upstream annotation pipeline)
- Community feedback from large-scale EDTA users on HPC and consortium genome projects



## References

- **[R1]** Lee, S.-I., & Kim, N.-S. (2014). *Transposable elements and genome size variations in plants.* **Genomics & Informatics**, 12(3), 87–97. https://doi.org/10.5808/GI.2014.12.3.87  
- **[R2]** Quesneville, H. (2020). *Twenty years of transposable element analysis in the Arabidopsis thaliana genome.* **Mobile DNA**, 11, 28. https://doi.org/10.1186/s13100-020-00223-x  
- **[R3]** Ou, S., Su, W., Liao, Y., Chougule, K., Ware, D., Peterson, T., Jiang, N., Hirsch, C. N., & Hufford, M. B. (2019). *Benchmarking transposable element annotation methods for creation of a streamlined, comprehensive pipeline.* **Genome Biology**, 20, 275. https://doi.org/10.1186/s13059-019-1905-y  
- **[R4]** Han, M., Perkins, M. H., Novaes, L. S., Xu, T., & Chang, H. (2023). *Advances in transposable elements: from mechanisms to applications in mammalian genomics.* **Frontiers in Genetics**, 14, 1290146. https://doi.org/10.3389/fgene.2023.1290146  
- **[R5]** Liu, Y.-N., Gao, J.-J., Zhuang, X.-L., Wu, D.-D., & Sun, Y.-B. (2025). *Near complete assembly of Drosophila melanogaster Canton S strain genome.* **Nature Communications**, 17, 329. https://doi.org/10.1038/s41467-025-67031-w  
- **[R6]** Laricchia, K. M., Zdraljevic, S., Cook, D. E., & Andersen, E. C. (2017). *Natural Variation in the Distribution and Abundance of Transposable Elements Across the Caenorhabditis elegans Species.* **Molecular Biology and Evolution**, 34(9), 2187–2202. https://doi.org/10.1093/molbev/msx155  
- **[R7]** Du, A. Y., Chobirko, J. D., Zhuo, X., Feschotte, C., Wang, T., et al. (2024). *Regulatory transposable elements in the encyclopedia of DNA elements.* **Nature Communications**, 15, 7594. https://doi.org/10.1038/s41467-024-51921-6  
- **[R8]** Chang, N.-C., Rovira, Q., Wells, J. N., & Feschotte, C. (2022). *Zebrafish transposable elements show extensive diversification in age, genomic distribution, and developmental expression.* **Genome Research**, 32(7), 1408–1423. https://doi.org/10.1101/gr.275655.121  
- **[R9]** Abraham, L. N., Oggenfuss, U., & Croll, D. (2024). *Population-level transposable element expression dynamics influence trait evolution in a fungal crop pathogen.* **mBio**, 15(3), e02840-23. https://doi.org/10.1128/mbio.02840-23  
- **[R10]** Romeijn, J., Bañales, I., Seidl, M. F., et al. (2025). *Extensive horizontal transfer of transposable elements shapes fungal mobilomes.* **Current Biology**. https://doi.org/10.1016/j.cub.2025.11.012  














