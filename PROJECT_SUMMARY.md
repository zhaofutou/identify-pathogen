# BioLab Project — AI Agent Work & Handoff Summary

**Working directory:** `C:\Users\rdpuser\bacterial_assemble`
**Date:** 2026-06-11
**Purpose:** Comprehensive handoff document for the next AI agent or developer.

---

## 1. Project Overview

A **Flask web application** (`webapp/app_docker.py`, running as `app.py` in the container) that serves as a clinical microbiology and genomics platform, featuring:

1. **Bacterial Genome Assembly**: Upload R1/R2 paired-end Illumina reads → FastQC quality control → SPAdes assembly → QUAST evaluation → Prokka annotation.
2. **Outbreak Investigation**: Metagenomic screening of clinical/environmental samples using a panel of 16 diarrheal pathogens via two modes:
   * **minimap2 Mode**: Reference-based read mapping against the 16 target pathogen genomes.
   * **Kraken2 Mode (Hybrid)**: Shotgun taxonomic classification against the 8 GB PLUSPF database, dynamically validated on-the-fly using `minimap2` alignments to resolve homology false positives.

---

## 2. Infrastructure & Environment

* **Host OS**: Windows 11 / Windows 10 LTSC.
* **WSL2 Distro**: Ubuntu (default user `root`, systemd enabled).
* **Docker Engine**: Installed natively inside WSL2 (not Docker Desktop).
* **Port Configuration**: App runs on **port 5050** (since port 5000 is occupied by Windows `svchost`).

### Key Environmental Gotchas & Fixes:
1. **WSL2 Idle Timeout VM Shutdown**: WSL2 automatically shuts down the VM when no interactive `wsl.exe` consoles are active on Windows. Since background Docker containers are not recognized as client sessions, WSL2 was shutting down mid-analysis, killing the container.
   * *Fix*: Started a persistent background keep-alive loop (`while true; do sleep 3600; done`) inside WSL2 via [run_biolab.bat](file:///C:/Users/rdpuser/bacterial_assemble/run_biolab.bat).
2. **Path Translation & Permissions**: Native docker-compose volumes are used for uploads, work, and results to avoid slow, buggy Windows-to-Linux path translation (`/mnt/c` bind mounts) and permissions errors.
3. **Reference Genomes & Kraken2 DB**:
   * Pre-downloaded fasta references from the WSL2 host (`/root/biolab/refs`) are mounted read-only into `/data/refs` to avoid slow, rate-limited NCBI downloads.
   * The 8 GB Kraken2 database is stored locally in the project directory (`C:\Users\rdpuser\bacterial_assemble\kraken2_db`) and mounted read-only into `/data/kraken2`.

---

## 3. Genomic Pipelines

### 3a. Genome Assembly (Dockerized)
All assembly tools are fully installed inside the Docker image:
* **FastQC** (v0.11.9)
* **SPAdes** (v3.13.1)
* **QUAST** (v5.2.0, installed via pip)
* **Prokka** (v1.14.6)
The pipeline is invoked in a background thread via `run_assembly()` inside [app_docker.py](file:///C:/Users/rdpuser/bacterial_assemble/webapp/app_docker.py).

### 3b. Outbreak Investigation & Hybrid Verification
* **minimap2 Mode**: Maps sample reads against the panel of 16 references. Positive thresholds are: `mapping_rate > 1.0%`, `avg_depth > 3x`, and `genome_coverage > 10%`.
* **Kraken2 Mode (Hybrid Pipeline)**:
  1. Classifies taxonomy using the 8 GB database.
  2. Automatically filters out human host DNA (taxid `9606`), human generic (`9605`), and common lab vectors/contaminants.
  3. **Alignment Validation**: If a known pathogen (Salmonella, Shigella, Campylobacter, Vibrio, Clostridioides, Norovirus, Rotavirus, Adenovirus, Astrovirus, Sapovirus, etc.) is detected at >0.5% reads, the app runs a quick `minimap2` alignment to calculate its actual **genome coverage** and **depth**.
  4. **Homology Filtering**: If genome coverage is $<10\%$, the hit is classified as `Normal Flora / Homology (Low Cov)` to eliminate false positives. If coverage is $>10\%$, it is marked positive.
  5. High-coverage pathogens are sorted to the top of the suspects table.

---

## 4. File Inventory

```
C:\Users\rdpuser\bacterial_assemble\
├── webapp/
│   ├── app_docker.py             ← Main Flask app (Docker Edition, hybrid validation, pathogen sorting)
│   └── templates/
│       ├── index.html            ← Landing home page
│       ├── assembly.html         ← Genome assembly page
│       ├── outbreak.html         ← Outbreak screening (dynamic headers, reads count, coverage)
│       └── results.html          ← Job results reporting
├── kraken2_db/                 ← The 8 GB Kraken2 database folder (ignored by Git)
├── Dockerfile                    ← Ubuntu 22.04 + minimap2 + samtools + bedtools + kraken2 + fastqc + spades + prokka + quast
├── docker-compose.yml            ← Compose config (ports, volume mounts for reads, refs, kraken2 db, and uploads)
├── entrypoint.sh                 ← Container entrypoint (web/shell/test modes)
├── run_biolab.bat                ← Start script (launches WSL2 keep-alive, docker-compose up, and downs container on close)
├── pathogen_panel.py             ← Standalone pathogen panel script
├── run_outbreak_test.py          ← Working direct-WSL2 test mapping script
├── enable_wsl.bat                ← Host setup script (WSL and VM platform enable)
└── PROJECT_SUMMARY.md            ← This document
```

---

## 5. Next Steps for the Next AI Agent

1. **Pathogen Panel Expansion**:
   * Add new accession numbers to `DIARRHEA_PANEL` and `PATHOGEN_PANEL` in [app_docker.py](file:///C:/Users/rdpuser/bacterial_assemble/webapp/app_docker.py) as needed.
   * When expanding, make sure to download the corresponding `.fasta` file to the host reference directory `/root/biolab/refs/` so it mounts automatically.
2. **Job Persistence**:
   * Replace the in-memory `jobs = {}` dictionary with a lightweight SQLite database so job history is preserved across container restarts.
3. **Genome Assembly Testing**:
   * Test the SPAdes -> Prokka pipeline inside Docker using the paired-end test files `2025SM060_{1,2}.clean.fq.gz` in the root directory.
