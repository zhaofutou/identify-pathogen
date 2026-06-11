#!/usr/bin/env python3
"""
Outbreak Investigation Pathogen Identification Pipeline
=========================================================
Targets common diarrheal pathogens using reference-based mapping.
This is the standard approach for real outbreak investigations.
"""

import subprocess
import os
import json
import sys
from collections import defaultdict

# Common diarrheal pathogens with their reference genomes (NCBI accession)
PATHOGEN_PANEL = {
    # Bacteria
    'Salmonella_enterica': 'NC_003197.2',
    'Escherichia_coli_O157_H7': 'NC_002655.2',
    'Shigella_flexneri_2a': 'NC_008258.1',
    'Campylobacter_jejuni': 'NC_002163.1',
    'Vibrio_cholerae_O1': 'NC_002505.1',
    'Clostridioides_difficile': 'NC_009089.1',
    'Listeria_monocytogenes': 'NC_003210.1',
    'Aeromonas_hydrophila': 'NC_008570.1',
    # Viruses
    'Norovirus_GII': 'NC_008228.1',
    'Rotavirus_A': 'NC_004526.1',
    'Adenovirus_40': 'NC_001454.1',
    'Astrovirus': 'NC_001943.1',
    # Parasites
    'Cryptosporidium_parvum': 'NC_006982.1',
}

WORK_DIR = '/root/outbreak_investigation/analysis'
DATA_DIR = '/root/outbreak_investigation/fastq_data'


def download_reference(accession, refs_dir):
    """Download reference genome from NCBI"""
    ref_file = os.path.join(refs_dir, f'{accession}.fasta')
    if os.path.exists(ref_file) and os.path.getsize(ref_file) > 100:
        return ref_file

    os.makedirs(refs_dir, exist_ok=True)
    url = (f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
           f'db=nucleotide&id={accession}&rettype=fasta&retmode=text')
    result = subprocess.run(
        ['curl', '-sL', '-o', ref_file, url, '--connect-timeout', '30', '--max-time', '120'],
        capture_output=True, text=True, timeout=130
    )
    if os.path.exists(ref_file) and os.path.getsize(ref_file) > 100:
        return ref_file
    return None


def map_reads_minimap2(reads_file, ref_file, outprefix):
    """Map Nanopore reads to reference with minimap2"""
    sam_file = f'{outprefix}.sam'
    bam_file = f'{outprefix}.bam'

    # Map
    cmd = f'minimap2 -ax map-ont -t 8 --eqx {ref_file} {reads_file} 2>/dev/null'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
    with open(sam_file, 'w') as f:
        f.write(result.stdout)

    # Sort and index
    subprocess.run(f'samtools sort -@ 4 -o {bam_file} {sam_file}', shell=True,
                   capture_output=True, timeout=120)
    subprocess.run(f'samtools index {bam_file}', shell=True, capture_output=True)

    # Flagstat
    result = subprocess.run(f'samtools flagstat {bam_file}', shell=True,
                            capture_output=True, text=True)
    flagstat = result.stdout

    # Depth
    result = subprocess.run(f'samtools depth {bam_file}', shell=True,
                            capture_output=True, text=True, timeout=120)
    depths = []
    for line in result.stdout.strip().split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 3:
                depths.append(int(parts[2]))

    avg_depth = sum(depths) / len(depths) if depths else 0
    covered_positions = sum(1 for d in depths if d > 0)

    # Parse flagstat for mapping rate
    mapped_reads = 0
    total_reads = 0
    for line in flagstat.split('\n'):
        if 'mapped (' in line and 'primary' not in line:
            try:
                mapped_reads = int(line.split('+')[0].strip())
            except:
                pass
        if '+ 0 in total' in line:
            try:
                total_reads = int(line.split('+')[0].strip())
            except:
                pass

    mapping_rate = (mapped_reads / total_reads * 100) if total_reads > 0 else 0

    # Get ref length
    ref_len = 0
    with open(ref_file) as f:
        for line in f:
            if not line.startswith('>'):
                ref_len += len(line.strip())

    genome_coverage = (covered_positions / ref_len * 100) if ref_len > 0 else 0

    # Cleanup
    os.remove(sam_file)

    return {
        'mapped_reads': mapped_reads,
        'total_reads': total_reads,
        'mapping_rate': round(mapping_rate, 2),
        'avg_depth': round(avg_depth, 1),
        'genome_coverage_pct': round(genome_coverage, 2),
        'ref_length': ref_len,
    }


def main():
    os.makedirs(WORK_DIR, exist_ok=True)
    refs_dir = os.path.join(WORK_DIR, 'refs')
    os.makedirs(refs_dir, exist_ok=True)

    # Find available samples
    samples = []
    for f in sorted(os.listdir(DATA_DIR)):
        if f.endswith('.fastq') and not f.endswith('.gz'):
            samples.append((f.replace('.fastq', ''), os.path.join(DATA_DIR, f)))

    if not samples:
        print('ERROR: No .fastq files found in', DATA_DIR)
        sys.exit(1)

    print(f'Found {len(samples)} samples to analyze')
    print(f'Pathogen panel: {len(PATHOGEN_PANEL)} targets')
    print()

    # Sample metadata (simulating outbreak scenario)
    sample_metadata = {}
    for i, (name, _) in enumerate(samples):
        if i < len(samples) // 2:
            sample_metadata[name] = 'clinical_stool'
        else:
            sample_metadata[name] = 'food_environmental'

    # Download all references first
    print('Downloading reference genomes...')
    for pathogen, accession in PATHOGEN_PANEL.items():
        ref = download_reference(accession, refs_dir)
        status = 'OK' if ref else 'FAILED'
        print(f'  {pathogen}: {status}')

    # Analyze each sample
    all_results = {}
    for sample_name, reads_file in samples:
        print(f'\n{"="*60}')
        print(f'Analyzing: {sample_name} ({sample_metadata.get(sample_name, "unknown")})')
        print(f'{"="*60}')

        sample_results = {}
        for pathogen, accession in PATHOGEN_PANEL.items():
            ref_file = os.path.join(refs_dir, f'{accession}.fasta')
            if not os.path.exists(ref_file):
                continue

            print(f'  Mapping to {pathogen}...', end=' ', flush=True)
            outprefix = os.path.join(WORK_DIR, f'{sample_name}_vs_{pathogen}')

            try:
                result = map_reads_minimap2(reads_file, ref_file, outprefix)
                sample_results[pathogen] = result

                if result['mapping_rate'] > 1:
                    print(f'MAPPED! {result["mapping_rate"]}% rate, {result["avg_depth"]}x depth')
                else:
                    print(f'negative ({result["mapping_rate"]}%)')
            except Exception as e:
                print(f'ERROR: {e}')

            # Cleanup BAM files
            for ext in ['.bam', '.bam.bai']:
                bam = f'{outprefix}{ext}'
                if os.path.exists(bam):
                    os.remove(bam)

        all_results[sample_name] = {
            'metadata': sample_metadata.get(sample_name, 'unknown'),
            'results': sample_results,
        }

    # ========== GENERATE REPORT ==========
    print('\n\n' + '='*70)
    print('  OUTBREAK INVESTIGATION - PATHOGEN IDENTIFICATION REPORT')
    print('='*70)
    print(f'  Samples analyzed: {len(all_results)}')
    print(f'  Pathogen panel: {len(PATHOGEN_PANEL)} targets')
    print(f'  Platform: Oxford Nanopore (MinION)')
    print()

    # Identify positive hits
    positive_hits = defaultdict(list)
    for sample_name, data in all_results.items():
        for pathogen, result in data['results'].items():
            if result['mapping_rate'] > 5 and result['avg_depth'] > 3:
                positive_hits[pathogen].append({
                    'sample': sample_name,
                    'type': data['metadata'],
                    'mapping_rate': result['mapping_rate'],
                    'avg_depth': result['avg_depth'],
                    'genome_coverage': result['genome_coverage_pct'],
                })

    if positive_hits:
        print('  SUSPECT PATHOGENS DETECTED:')
        print('  ' + '-'*65)
        for pathogen, hits in sorted(positive_hits.items(), key=lambda x: -len(x[1])):
            n_clinical = sum(1 for h in hits if h['type'] == 'clinical_stool')
            n_food = sum(1 for h in hits if h['type'] == 'food_environmental')
            avg_rate = sum(h['mapping_rate'] for h in hits) / len(hits)
            avg_depth = sum(h['avg_depth'] for h in hits) / len(hits)
            avg_cov = sum(h['genome_coverage'] for h in hits) / len(hits)

            print(f'\n  >>> {pathogen}')
            print(f'      Positive samples: {len(hits)}/{len(all_results)} '
                  f'(clinical: {n_clinical}, food/env: {n_food})')
            print(f'      Avg mapping rate: {avg_rate:.1f}%')
            print(f'      Avg depth: {avg_depth:.1f}x')
            print(f'      Avg genome coverage: {avg_cov:.1f}%')

            for h in hits:
                print(f'        - {h["sample"]} ({h["type"]}): '
                      f'{h["mapping_rate"]}% rate, {h["avg_depth"]}x depth, '
                      f'{h["genome_coverage"]}% coverage')
    else:
        print('  No significant pathogen hits detected.')
        print('  Showing top hits per sample:')
        for sample_name, data in all_results.items():
            if data['results']:
                top = max(data['results'].items(), key=lambda x: x[1]['mapping_rate'])
                print(f'    {sample_name}: best = {top[0]} ({top[1]["mapping_rate"]}%)')

    # Conclusion
    print('\n' + '='*70)
    print('  CONCLUSION')
    print('='*70)
    if positive_hits:
        # Find the pathogen present in the most samples
        top_pathogen = max(positive_hits.items(), key=lambda x: len(x[1]))
        in_clinical = sum(1 for h in top_pathogen[1] if h['type'] == 'clinical_stool')
        in_food = sum(1 for h in top_pathogen[1] if h['type'] == 'food_environmental')

        print(f'  The most likely causative agent is: {top_pathogen[0]}')
        print(f'  Detected in {len(top_pathogen[1])} of {len(all_results)} samples')
        if in_food > 0:
            print(f'  Found in BOTH clinical ({in_clinical}) AND food/environmental ({in_food}) samples')
            print(f'  => Foodborne transmission CONFIRMED')
        else:
            print(f'  Found only in clinical samples ({in_clinical})')
            print(f'  => Food/environmental source not confirmed in this dataset')
    else:
        print('  No clear pathogen identified. Consider:')
        print('  - Expanding the pathogen panel')
        print('  - Metagenomic shotgun analysis')
        print('  - Toxin testing')

    # Save results
    report_file = os.path.join(WORK_DIR, 'outbreak_report.json')
    with open(report_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f'\n  Full results saved to: {report_file}')
    print('='*70)


if __name__ == '__main__':
    main()
