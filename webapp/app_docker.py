#!/usr/bin/env python3
"""
BioLab Service — Genome Assembly & Outbreak Investigation (Docker Edition)
==========================================================================
All bioinformatics tools run natively inside the container.
No WSL, no path translation, no permission issues.
"""

import os
import sys
import json
import time
import uuid
import shutil
import subprocess
import threading
from datetime import datetime
from pathlib import Path
from flask import Flask, render_template, request, jsonify, send_file, redirect, url_for

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024 * 1024  # 50GB max upload

# ─── Configuration (Docker paths) ────────────────────────────────
DATA_DIR = Path('/data')
UPLOAD_DIR = DATA_DIR / 'uploads'
RESULTS_DIR = DATA_DIR / 'results'
WORK_DIR = DATA_DIR / 'work'
REFS_DIR = DATA_DIR / 'refs'
KRAKEN2_DB = DATA_DIR / 'kraken2'

for d in [UPLOAD_DIR, RESULTS_DIR, WORK_DIR, REFS_DIR, KRAKEN2_DB]:
    d.mkdir(parents=True, exist_ok=True)

# Track running jobs
jobs = {}  # job_id -> {status, progress, result, ...}

# ─── Pathogen Reference Panel ───────────────────────────────────
PATHOGEN_PANEL = {
    'Salmonella_enterica': 'NC_003197.2',
    'Escherichia_coli_O157_H7': 'NC_002655.2',
    'Shigella_flexneri_2a': 'NC_008258.1',
    'Campylobacter_jejuni': 'NC_002163.1',
    'Vibrio_cholerae_O1': 'NC_002505.1',
    'Clostridioides_difficile': 'NC_009089.1',
    'Listeria_monocytogenes': 'NC_003210.1',
    'Aeromonas_hydrophila': 'NC_008570.1',
    'Norovirus_GII': 'NC_008228.1',
    'Rotavirus_A': 'NC_004526.1',
    'Adenovirus_40': 'NC_001454.1',
    'Astrovirus': 'NC_001943.1',
    'Cryptosporidium_parvum': 'NC_006982.2',
}

# Expanded diarrhea panel
DIARRHEA_PANEL = {
    # Bacteria
    'Salmonella_enterica': 'NC_003197.2',
    'Escherichia_coli_O157_H7': 'NC_002655.2',
    'Shigella_flexneri_2a': 'NC_008258.1',
    'Campylobacter_jejuni': 'NC_002163.1',
    'Vibrio_cholerae_O1': 'NC_002505.1',
    'Clostridioides_difficile': 'NC_009089.1',
    'Listeria_monocytogenes': 'NC_003210.1',
    'Yersinia_enterocolitica': 'NC_008800.1',
    'Staphylococcus_aureus': 'NC_007795.1',
    'Bacillus_cereus': 'NC_004722.1',
    'Aeromonas_hydrophila': 'NC_008570.1',
    # Viruses
    'Norovirus_GII': 'NC_008228.1',
    'Rotavirus_A': 'NC_004526.1',
    'Adenovirus_40': 'NC_001454.1',
    'Astrovirus': 'NC_001943.1',
    'Sapovirus': 'NC_006269.1',
}


# ─── Helper Functions ────────────────────────────────────────────
def run_cmd(cmd, timeout=3600):
    """Run a shell command and return output."""
    try:
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=timeout
        )
        return result.stdout + result.stderr, result.returncode
    except subprocess.TimeoutExpired:
        return 'Command timed out', 1
    except Exception as e:
        return str(e), 1


def update_job(job_id, **kwargs):
    """Update job status."""
    if job_id in jobs:
        jobs[job_id].update(kwargs)


def download_reference(accession):
    """Download reference genome from NCBI."""
    ref_file = REFS_DIR / f'{accession}.fasta'
    if ref_file.exists() and ref_file.stat().st_size > 100:
        return ref_file

    url = (f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
           f'db=nucleotide&id={accession}&rettype=fasta&retmode=text')
    out, rc = run_cmd(
        f'curl -sL -o {ref_file} "{url}" --connect-timeout 30 --max-time 120',
        timeout=130
    )
    if ref_file.exists() and ref_file.stat().st_size > 100:
        return ref_file
    return None


def is_kraken2_ready():
    """Check if Kraken2 database is available."""
    return (KRAKEN2_DB / 'hash.k2d').exists()


# ─── Genome Assembly Pipeline ────────────────────────────────────
def run_assembly(job_id, r1_path, r2_path, threads=16, memory=32):
    """Run the full genome assembly pipeline."""
    try:
        update_job(job_id, status='running', step='Preparing files')

        job_dir = RESULTS_DIR / job_id
        job_dir.mkdir(parents=True, exist_ok=True)

        links_dir = WORK_DIR / job_id
        links_dir.mkdir(parents=True, exist_ok=True)

        r1_link = links_dir / 'R1.fastq.gz'
        r2_link = links_dir / 'R2.fastq.gz'
        shutil.copy2(r1_path, r1_link)
        shutil.copy2(r2_path, r2_link)

        # Step 1: FastQC
        update_job(job_id, step='Running FastQC', progress=10)
        fastqc_dir = job_dir / 'fastqc'
        fastqc_dir.mkdir(exist_ok=True)
        run_cmd(f'fastqc -t {threads} -o {fastqc_dir} {r1_link} {r2_link}')

        # Step 2: SPAdes assembly
        update_job(job_id, step='Running SPAdes assembly', progress=25)
        spades_dir = job_dir / 'spades'
        spades_cmd = (
            f'spades.py --isolate '
            f'-1 {r1_link} -2 {r2_link} '
            f'-o {spades_dir} '
            f'-t {threads} -m {memory} '
            f'--cov-cutoff auto 2>&1'
        )
        run_cmd(spades_cmd, timeout=1800)

        # Step 3: QUAST evaluation
        update_job(job_id, step='Running QUAST evaluation', progress=60)
        quast_dir = job_dir / 'quast'
        quast_cmd = (
            f'python3 /usr/local/bin/quast.py '
            f'{spades_dir}/scaffolds.fasta '
            f'-o {quast_dir} -t {threads} 2>&1'
        )
        run_cmd(quast_cmd, timeout=300)

        # Step 4: Prokka annotation
        update_job(job_id, step='Running Prokka annotation', progress=70)
        prokka_dir = job_dir / 'prokka'
        prokka_cmd = (
            f'prokka --outdir {prokka_dir} '
            f'--genus Salmonella --species enterica '
            f'--cpus {threads} '
            f'--prefix assembly '
            f'{spades_dir}/scaffolds.fasta 2>&1'
        )
        run_cmd(prokka_cmd, timeout=600)

        # Step 5: Collect stats
        update_job(job_id, step='Generating report', progress=90)
        stats = collect_assembly_stats(spades_dir, quast_dir, prokka_dir)

        with open(job_dir / 'results.json', 'w') as f:
            json.dump(stats, f, indent=2)

        update_job(job_id, status='completed', progress=100, result=stats,
                   step='Complete')

    except Exception as e:
        update_job(job_id, status='error', error=str(e))


def collect_assembly_stats(spades_dir, quast_dir, prokka_dir):
    """Collect assembly statistics."""
    stats = {
        'timestamp': datetime.now().isoformat(),
        'assembly': {},
        'quast': {},
        'annotation': {},
    }

    scaffolds = spades_dir / 'scaffolds.fasta'
    if scaffolds.exists():
        quast_report = quast_dir / 'report.tsv'
        if quast_report.exists():
            quast_raw = {}
            with open(quast_report) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        key = parts[0].strip().lower()
                        val = parts[1].strip()
                        quast_raw[key] = val
            total = 0
            try: total = int(quast_raw.get('total length', '0').replace(',', ''))
            except: pass
            n50 = 0
            try: n50 = int(quast_raw.get('n50', '0').replace(',', ''))
            except: pass
            longest = 0
            try: longest = int(quast_raw.get('largest contig', '0').replace(',', ''))
            except: pass
            num_scaffolds = quast_raw.get('# contigs', '-')
            gc = quast_raw.get('gc (%)', '-')
            stats['assembly'] = {
                'total_length': total, 'num_scaffolds': num_scaffolds,
                'n50': n50, 'longest': longest, 'gc': gc,
            }

    quast_report = quast_dir / 'report.tsv'
    if quast_report.exists():
        with open(quast_report) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    key = parts[0].strip().lower().replace(' ', '_').replace('#', 'num')
                    val = parts[1].strip()
                    if val.replace('.', '').replace(',', '').isdigit():
                        stats['quast'][key] = val

    prokka_txt = prokka_dir / 'assembly.txt'
    if prokka_txt.exists():
        with open(prokka_txt) as f:
            for line in f:
                if ':' in line:
                    key, val = line.strip().split(':', 1)
                    stats['annotation'][key.strip()] = val.strip()

    return stats


# ─── Outbreak Investigation Pipeline ─────────────────────────────
def run_outbreak(job_id, sample_files, sample_metadata, threads=8, mode='diarrhea'):
    """Run pathogen identification pipeline.
    
    mode: 'kraken2' for comprehensive identification, 'diarrhea' for fast diarrhea screen
    """
    try:
        update_job(job_id, status='running', step='Preparing samples')

        job_dir = RESULTS_DIR / job_id
        job_dir.mkdir(parents=True, exist_ok=True)

        samples_dir = WORK_DIR / job_id / 'samples'
        samples_dir.mkdir(parents=True, exist_ok=True)

        sample_paths = []
        for i, (fname, fpath) in enumerate(sample_files):
            dest = samples_dir / fname
            shutil.copy2(fpath, dest)
            sample_paths.append((fname, dest, sample_metadata.get(i, 'unknown')))

        # Determine method
        use_kraken2 = False
        if mode == 'kraken2' and is_kraken2_ready():
            use_kraken2 = True
            panel = PATHOGEN_PANEL
            update_job(job_id, step='Running Kraken2 classification')
        elif mode == 'diarrhea':
            panel = DIARRHEA_PANEL
            update_job(job_id, step='Fast diarrhea screen (minimap2)')
        else:
            panel = PATHOGEN_PANEL
            update_job(job_id, step='Running minimap2 panel')

        results = {}
        total = len(sample_paths)

        for idx, (name, path, stype) in enumerate(sample_paths):
            progress = 10 + int(80 * idx / total)
            update_job(job_id,
                       step=f'Analyzing sample {idx+1}/{total}: {name}',
                       progress=progress)

            sample_result = {'name': name, 'type': stype, 'pathogens': {}}

            if use_kraken2:
                sample_result['pathogens'] = run_kraken2(path, job_dir, name, threads)
                sample_result['method'] = 'kraken2'
            else:
                sample_result['pathogens'] = run_minimap2_panel(
                    path, job_dir, name, threads, panel
                )
                sample_result['method'] = 'minimap2'

            results[name] = sample_result

        # Generate summary
        update_job(job_id, step='Generating outbreak report', progress=95)
        summary = generate_outbreak_summary(results)

        with open(job_dir / 'results.json', 'w') as f:
            json.dump({'samples': results, 'summary': summary}, f, indent=2)

        update_job(job_id, status='completed', progress=100,
                   result={'samples': results, 'summary': summary},
                   step='Complete')

    except Exception as e:
        update_job(job_id, status='error', error=str(e))


def get_panel_accession(name):
    """Map a taxon name to the NCBI accession of our reference panel."""
    name_lower = name.lower()
    if 'adenovirus' in name_lower or 'mastadenovirus' in name_lower:
        return 'NC_001454.1'
    elif 'campylobacter' in name_lower:
        return 'NC_002163.1'
    elif 'escherichia' in name_lower:
        return 'NC_002655.2'
    elif 'salmonella' in name_lower:
        return 'NC_003197.2'
    elif 'shigella' in name_lower:
        return 'NC_008258.1'
    elif 'vibrio' in name_lower:
        return 'NC_002505.1'
    elif 'clostridioides' in name_lower or 'difficile' in name_lower:
        return 'NC_009089.1'
    elif 'listeria' in name_lower:
        return 'NC_003210.1'
    elif 'yersinia' in name_lower:
        return 'NC_008800.1'
    elif 'staphylococcus' in name_lower:
        return 'NC_007795.1'
    elif 'bacillus' in name_lower:
        return 'NC_004722.1'
    elif 'aeromonas' in name_lower:
        return 'NC_008570.1'
    elif 'norovirus' in name_lower:
        return 'NC_008228.1'
    elif 'rotavirus' in name_lower:
        return 'NC_004526.1'
    elif 'astrovirus' in name_lower:
        return 'NC_001943.1'
    elif 'sapovirus' in name_lower:
        return 'NC_006269.1'
    return None


def get_alignment_coverage(reads_path, accession, threads=4):
    """Run a quick minimap2 alignment to compute coverage and depth for validation."""
    ref_file = REFS_DIR / f'{accession}.fasta'
    if not ref_file.exists():
        return 0.0, 0.0
        
    outprefix = WORK_DIR / f'tmp_{uuid.uuid4().hex[:6]}'
    sam = f'{outprefix}.sam'
    bam = f'{outprefix}.bam'
    
    # Map reads
    run_cmd(
        f'minimap2 -ax map-ont -t {threads} {ref_file} {reads_path} > {sam}',
        timeout=120
    )
    
    # Check if we got alignments
    sam_check, _ = run_cmd(f'wc -l < {sam}')
    try: sam_lines = int(sam_check.strip())
    except: sam_lines = 0
    
    if sam_lines < 10:
        run_cmd(
            f'minimap2 -ax sr -t {threads} {ref_file} {reads_path} > {sam}',
            timeout=120
        )
        
    run_cmd(f'samtools sort -@ 2 -o {bam} {sam}', timeout=60)
    
    depth_out, _ = run_cmd(f'samtools depth {bam}', timeout=60)
    
    depths = []
    for line in depth_out.strip().split('\n')[:100000]:
        parts = line.split('\t')
        if len(parts) >= 3:
            try: depths.append(int(parts[2]))
            except: pass
            
    avg_depth = sum(depths) / len(depths) if depths else 0
    
    ref_len_out, _ = run_cmd(f'grep -v "^>" {ref_file} | tr -d "\\n" | wc -c')
    try: ref_len = int(ref_len_out.strip())
    except: ref_len = 0
    
    covered = sum(1 for d in depths if d > 0)
    genome_cov = (covered / ref_len * 100) if ref_len > 0 else 0
    
    # Cleanup
    for ext in ['.sam', '.bam', '.bam.bai']:
        f = Path(f'{outprefix}{ext}')
        if f.exists():
            f.unlink()
            
    return round(genome_cov, 2), round(avg_depth, 1)


def run_kraken2(reads_path, job_dir, sample_name, threads):
    """Run Kraken2 classification."""
    results = {}
    report_file = job_dir / f'{sample_name}_kraken2_report.txt'
    output_file = job_dir / f'{sample_name}_kraken2_output.txt'

    cmd = (
        f'kraken2 --db {KRAKEN2_DB} '
        f'--threads {threads} '
        f'--report {report_file} '
        f'--output {output_file} '
        f'{reads_path} 2>&1'
    )
    out, rc = run_cmd(cmd, timeout=600)

    if report_file.exists():
        with open(report_file) as f:
            report_text = f.read()
        results = parse_kraken2_report(report_text)
        
        # Run hybrid validation (minimap2 mapping) on the key candidates
        for name, pdata in list(results.items()):
            pct = pdata.get('percentage', 0)
            if pct > 0.5:
                accession = get_panel_accession(name)
                if accession:
                    print(f"  Validating {name} ({accession}) with minimap2 alignment...")
                    cov, depth = get_alignment_coverage(reads_path, accession, threads)
                    pdata['genome_coverage'] = cov
                    pdata['avg_depth'] = depth
                    print(f"    => Coverage: {cov}%, Depth: {depth}x")

    return results


def parse_kraken2_report(report_text):
    """Parse Kraken2 report — only species/genus/strain level hits."""
    results = {}
    for line in report_text.strip().split('\n'):
        parts = line.strip().split('\t')
        if len(parts) >= 6:
            pct = float(parts[0])
            reads = int(parts[1])
            rank = parts[3].strip()
            taxid = parts[4]
            name = parts[5].strip()
            
            # Filter out host/human, synthetic constructs, and common non-pathogen contaminants
            name_lower = name.lower()
            is_host_or_contaminant = (
                'homo' in name_lower or
                'sapiens' in name_lower or
                'human' in name_lower or
                'synthetic' in name_lower or
                'artificial' in name_lower or
                'vector' in name_lower or
                taxid in ('9606', '9605', '9600')
            )
            if is_host_or_contaminant:
                continue
                
            if rank in ('S', 'S1', 'G') and pct > 0.1 and taxid != '0':
                results[name] = {
                    'percentage': pct,
                    'reads': reads,
                    'taxid': taxid,
                    'rank': rank,
                }
    return results


def run_minimap2_panel(reads_path, job_dir, sample_name, threads, panel=None):
    """Map reads against pathogen panel using minimap2."""
    if panel is None:
        panel = PATHOGEN_PANEL
    results = {}

    reads_path = Path(reads_path)
    if not reads_path.exists():
        print(f'[ERROR] Reads file not found: {reads_path}')
        return results

    for pathogen, accession in panel.items():
        ref_file = REFS_DIR / f'{accession}.fasta'

        if not ref_file.exists() or ref_file.stat().st_size < 100:
            print(f'  Downloading reference for {pathogen}...')
            ref = download_reference(accession)
            if ref is None:
                results[pathogen] = {
                    'accession': accession, 'mapping_rate': 0, 'avg_depth': 0,
                    'genome_coverage': 0, 'mapped_reads': 0, 'total_reads': 0,
                    'error': 'Reference download failed',
                }
                continue

        outprefix = job_dir / f'{sample_name}_vs_{pathogen}'
        sam = f'{outprefix}.sam'
        bam = f'{outprefix}.bam'

        # Map reads
        map_out, map_rc = run_cmd(
            f'minimap2 -ax map-ont -t {threads} {ref_file} {reads_path} > {sam}',
            timeout=600
        )

        sam_check, _ = run_cmd(f'wc -l < {sam}')
        try: sam_lines = int(sam_check.strip())
        except: sam_lines = 0

        if sam_lines < 10:
            map_out, map_rc = run_cmd(
                f'minimap2 -ax sr -t {threads} {ref_file} {reads_path} > {sam}',
                timeout=600
            )

        run_cmd(f'samtools sort -@ 4 -o {bam} {sam}', timeout=300)
        run_cmd(f'samtools index {bam}')

        flagstat, _ = run_cmd(f'samtools flagstat {bam}')
        depth_out, _ = run_cmd(f'samtools depth {bam}', timeout=300)

        mapped = 0
        total = 0
        for line in flagstat.split('\n'):
            if 'mapped (' in line and 'primary' not in line and 'singleton' not in line:
                try: mapped = int(line.split('+')[0].strip())
                except: pass
            if '+ 0 in total' in line:
                try: total = int(line.split('+')[0].strip())
                except: pass

        depths = []
        for line in depth_out.strip().split('\n')[:100000]:
            parts = line.split('\t')
            if len(parts) >= 3:
                try: depths.append(int(parts[2]))
                except: pass

        avg_depth = sum(depths) / len(depths) if depths else 0
        mapping_rate = (mapped / total * 100) if total > 0 else 0

        ref_len_out, _ = run_cmd(f'grep -v "^>" {ref_file} | tr -d "\\n" | wc -c')
        try: ref_len = int(ref_len_out.strip())
        except: ref_len = 0

        covered = sum(1 for d in depths if d > 0)
        genome_cov = (covered / ref_len * 100) if ref_len > 0 else 0

        results[pathogen] = {
            'accession': accession,
            'mapping_rate': round(mapping_rate, 2),
            'avg_depth': round(avg_depth, 1),
            'genome_coverage': round(genome_cov, 2),
            'mapped_reads': mapped,
            'total_reads': total,
        }

        if mapped > 0:
            print(f'  [HIT] {pathogen}: {mapped}/{total} mapped ({mapping_rate:.1f}%), '
                  f'depth={avg_depth:.1f}x, coverage={genome_cov:.1f}%')

        for ext in ['.sam', '.bam', '.bam.bai']:
            f = Path(f'{outprefix}{ext}')
            if f.exists():
                f.unlink()

    return results


def is_known_pathogen(name):
    """Check if the organism name matches a known diarrheal pathogen."""
    pathogen_keywords = [
        'salmonella', 'shigella', 'campylobacter', 'vibrio', 'cholerae', 
        'clostridioides', 'difficile', 'listeria', 'yersinia', 'staphylococcus', 
        'aureus', 'bacillus', 'cereus', 'aeromonas', 'norovirus', 'rotavirus', 
        'adenovirus', 'mastadenovirus', 'astrovirus', 'sapovirus', 'escherichia'
    ]
    name_lower = name.lower()
    return any(keyword in name_lower for keyword in pathogen_keywords)


def generate_outbreak_summary(results):
    """Generate a summary of the outbreak investigation."""
    summary = {
        'total_samples': len(results),
        'suspects': [],
        'top_hit_per_sample': {},
    }

    pathogen_hits = {}

    for sname, sdata in results.items():
        best = None
        best_effective_score = 0

        for pname, pdata in sdata['pathogens'].items():
            if sdata.get('method') == 'kraken2':
                score = pdata.get('percentage', 0)
                is_positive = score > 0.5
                
                # Use genome coverage if available to adjust the effective score
                cov = pdata.get('genome_coverage', 0)
                if cov > 0:
                    coverage_factor = 1.0 if cov > 10.0 else 0.01
                else:
                    coverage_factor = 1.0
                
                # Boost effective score for known pathogens to select them as best hit
                effective_score = score * 100 * coverage_factor if is_known_pathogen(pname) else score * coverage_factor
            else:
                rate = pdata.get('mapping_rate', 0)
                depth = pdata.get('avg_depth', 0)
                cov = pdata.get('genome_coverage', 0)
                score = depth
                is_positive = (rate > 1.0 and depth > 3 and cov > 10)
                effective_score = score

            if effective_score > best_effective_score:
                best_effective_score = effective_score
                best = pname

            if is_positive:
                if pname not in pathogen_hits:
                    pathogen_hits[pname] = []
                pathogen_hits[pname].append({
                    'sample': sname,
                    'mapping_rate': pdata.get('mapping_rate', pdata.get('percentage', 0)),
                    'avg_depth': pdata.get('avg_depth', pdata.get('percentage', 0)),
                    'genome_coverage': pdata.get('genome_coverage', 0),
                })

        if best:
            summary['top_hit_per_sample'][sname] = best

    for pname, hits in sorted(pathogen_hits.items(),
                               key=lambda x: (
                                   not is_known_pathogen(x[0]),  # Known pathogens first (False = 0)
                                   -len(x[1]),                    # Most positive samples first
                                   -sum(h['avg_depth'] for h in x[1]) / len(x[1])  # Highest abundance first
                               )):
        avg_rate = sum(h['mapping_rate'] for h in hits) / len(hits)
        avg_depth = sum(h['avg_depth'] for h in hits) / len(hits)
        avg_cov = sum(h['genome_coverage'] for h in hits) / len(hits)
        summary['suspects'].append({
            'pathogen': pname,
            'positive_samples': len(hits),
            'samples': [h['sample'] for h in hits],
            'avg_mapping_rate': round(avg_rate, 2),
            'avg_depth': round(avg_depth, 1),
            'avg_genome_coverage': round(avg_cov, 1),
        })

    return summary


# ─── Flask Routes ────────────────────────────────────────────────
@app.route('/')
def index():
    return render_template('index.html')


@app.route('/assembly')
def assembly_page():
    return render_template('assembly.html')


@app.route('/outbreak')
def outbreak_page():
    mode = request.args.get('mode', 'diarrhea')
    if mode == 'diarrhea':
        panel = DIARRHEA_PANEL
    else:
        panel = PATHOGEN_PANEL
    return render_template('outbreak.html',
                         kraken2_ready=is_kraken2_ready(),
                         pathogen_panel=panel, mode=mode)


@app.route('/api/assembly/start', methods=['POST'])
def start_assembly():
    """Start a genome assembly job."""
    files = request.files.getlist('files')
    if len(files) < 2:
        return jsonify({'error': 'Please upload 2 files (R1 and R2)'}), 400

    job_id = f"asm_{uuid.uuid4().hex[:8]}"
    job_dir = UPLOAD_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    paths = []
    for f in files:
        fpath = job_dir / f.filename
        f.save(str(fpath))
        paths.append(fpath)

    paths.sort(key=lambda p: p.name)

    threads = int(request.form.get('threads', 16))
    memory = int(request.form.get('memory', 32))

    jobs[job_id] = {
        'type': 'assembly', 'status': 'queued', 'progress': 0,
        'step': 'Queued', 'files': [p.name for p in paths],
        'created': datetime.now().isoformat(),
    }

    t = threading.Thread(
        target=run_assembly,
        args=(job_id, str(paths[0]), str(paths[1]), threads, memory),
        daemon=True
    )
    t.start()

    return jsonify({'job_id': job_id, 'status': 'queued'})


@app.route('/api/outbreak/start', methods=['POST'])
def start_outbreak():
    """Start an outbreak investigation job."""
    files = request.files.getlist('files')
    if not files:
        return jsonify({'error': 'Please upload at least one sample'}), 400

    job_id = f"out_{uuid.uuid4().hex[:8]}"
    job_dir = UPLOAD_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    sample_files = []
    sample_metadata = {}
    for i, f in enumerate(files):
        fpath = job_dir / f.filename
        f.save(str(fpath))
        sample_files.append((f.filename, str(fpath)))
        stype = request.form.get(f'type_{i}', 'unknown')
        sample_metadata[i] = stype

    threads = int(request.form.get('threads', 8))
    mode = request.form.get('mode', 'diarrhea')

    jobs[job_id] = {
        'type': 'outbreak', 'mode': mode, 'status': 'queued',
        'progress': 0, 'step': 'Queued',
        'files': [f[0] for f in sample_files],
        'created': datetime.now().isoformat(),
    }

    t = threading.Thread(
        target=run_outbreak,
        args=(job_id, sample_files, sample_metadata, threads, mode),
        daemon=True
    )
    t.start()

    return jsonify({'job_id': job_id, 'status': 'queued'})


@app.route('/api/job/<job_id>')
def get_job(job_id):
    if job_id not in jobs:
        return jsonify({'error': 'Job not found'}), 404
    return jsonify(jobs[job_id])


@app.route('/api/jobs')
def list_jobs():
    return jsonify(jobs)


@app.route('/results/<job_id>')
def results_page(job_id):
    if job_id not in jobs:
        return redirect(url_for('index'))
    return render_template('results.html', job_id=job_id)


@app.route('/api/results/<job_id>/json')
def results_json(job_id):
    rfile = RESULTS_DIR / job_id / 'results.json'
    if rfile.exists():
        with open(rfile) as f:
            return jsonify(json.load(f))
    return jsonify({'error': 'No results yet'}), 404


@app.route('/api/results/<job_id>/download/<filename>')
def download_file(job_id, filename):
    fpath = RESULTS_DIR / job_id / filename
    if fpath.exists():
        return send_file(str(fpath), as_attachment=True)
    return jsonify({'error': 'File not found'}), 404


@app.route('/api/kraken2/status')
def kraken2_status():
    """Check if Kraken2 database is ready."""
    ready = is_kraken2_ready()
    db_path = str(KRAKEN2_DB)
    size = 0
    if ready:
        try:
            size = sum(f.stat().st_size for f in KRAKEN2_DB.iterdir())
        except:
            pass
    return jsonify({
        'ready': ready,
        'database_path': db_path,
        'database_size_gb': round(size / 1024**3, 2),
    })


if __name__ == '__main__':
    print("=" * 60)
    print("  BioLab Service — Docker Edition")
    print("=" * 60)
    print(f"  Data dir:    {DATA_DIR}")
    print(f"  Refs dir:    {REFS_DIR}")
    print(f"  Kraken2 DB:  {KRAKEN2_DB} ({'✅ Ready' if is_kraken2_ready() else '⚠️  Not found — mount to /data/kraken2'})")
    print(f"  Results dir: {RESULTS_DIR}")
    print(f"  Access:      http://0.0.0.0:5050")
    print("=" * 60)

    app.run(host='0.0.0.0', port=5050, debug=False, threaded=True)
