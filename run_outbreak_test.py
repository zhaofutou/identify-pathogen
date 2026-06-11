#!/usr/bin/env python3
"""
Quick outbreak analysis — runs the same pipeline as the web app
but directly in WSL2 for faster results.
"""
import subprocess
import os
import json
import sys
from collections import defaultdict
from datetime import datetime

PANEL = {
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
    'Norovirus_GII': 'NC_008228.1',
    'Rotavirus_A': 'NC_004526.1',
    'Adenovirus_40': 'NC_001454.1',
    'Astrovirus': 'NC_001943.1',
    'Sapovirus': 'NC_006269.1',
}

SAMPLES = [
    ('SRR35602535', '/mnt/c/Users/rdpuser/bacterial_assemble/outbreak_data/SRR35602535.fastq.gz', 'clinical_stool'),
    ('SRR35602538', '/mnt/c/Users/rdpuser/bacterial_assemble/outbreak_data/SRR35602538.fastq.gz', 'clinical_stool'),
    ('SRR35602539', '/mnt/c/Users/rdpuser/bacterial_assemble/outbreak_data/SRR35602539.fastq.gz', 'food_environmental'),
    ('SRR35602540', '/mnt/c/Users/rdpuser/bacterial_assemble/outbreak_data/SRR35602540.fastq.gz', 'food_environmental'),
]

REFS_DIR = '/root/biolab/refs'
WORK_DIR = '/root/biolab/work'
THREADS = 8


def download_ref(accession):
    ref_file = f'{REFS_DIR}/{accession}.fasta'
    if os.path.exists(ref_file) and os.path.getsize(ref_file) > 100:
        return ref_file
    os.makedirs(REFS_DIR, exist_ok=True)
    url = (f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
           f'db=nucleotide&id={accession}&rettype=fasta&retmode=text')
    subprocess.run(['curl', '-sL', '-o', ref_file, url,
                    '--connect-timeout', '30', '--max-time', '120'],
                   capture_output=True, timeout=130)
    if os.path.exists(ref_file) and os.path.getsize(ref_file) > 100:
        return ref_file
    return None


def map_sample_to_pathogen(reads_file, ref_file, outprefix):
    sam_file = f'{outprefix}.sam'
    bam_file = f'{outprefix}.bam'

    with open(sam_file, 'w') as f:
        proc = subprocess.run(
            ['minimap2', '-ax', 'map-ont', '-t', str(THREADS), ref_file, reads_file],
            stdout=f, stderr=subprocess.PIPE, timeout=600
        )

    subprocess.run(['samtools', 'sort', '-@', '4', '-o', bam_file, sam_file],
                   capture_output=True, timeout=300)
    subprocess.run(['samtools', 'index', bam_file], capture_output=True)

    # Flagstat
    result = subprocess.run(['samtools', 'flagstat', bam_file],
                           capture_output=True, text=True)
    mapped = 0
    total = 0
    for line in result.stdout.split('\n'):
        if 'mapped (' in line and 'primary' not in line and 'singleton' not in line:
            try: mapped = int(line.split('+')[0].strip())
            except: pass
        if '+ 0 in total' in line:
            try: total = int(line.split('+')[0].strip())
            except: pass

    # Depth
    result = subprocess.run(['samtools', 'depth', bam_file],
                           capture_output=True, text=True, timeout=300)
    depths = []
    for line in result.stdout.strip().split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 3:
                try: depths.append(int(parts[2]))
                except: pass

    avg_depth = sum(depths) / len(depths) if depths else 0
    mapping_rate = (mapped / total * 100) if total > 0 else 0

    # Ref length
    ref_len = 0
    with open(ref_file) as rf:
        for line in rf:
            if not line.startswith('>'):
                ref_len += len(line.strip())

    covered = sum(1 for d in depths if d > 0)
    genome_cov = (covered / ref_len * 100) if ref_len > 0 else 0

    # Cleanup
    for f in [sam_file, bam_file, f'{bam_file}.bai']:
        if os.path.exists(f):
            os.remove(f)

    return {
        'mapped_reads': mapped,
        'total_reads': total,
        'mapping_rate': round(mapping_rate, 2),
        'avg_depth': round(avg_depth, 1),
        'genome_coverage': round(genome_cov, 2),
        'ref_length': ref_len,
    }


def main():
    os.makedirs(REFS_DIR, exist_ok=True)
    os.makedirs(WORK_DIR, exist_ok=True)

    print(f'{"="*70}')
    print(f'  OUTBREAK INVESTIGATION — BioLab Pipeline Test')
    print(f'  {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    print(f'  Samples: {len(SAMPLES)}  |  Panel: {len(PANEL)} targets')
    print(f'{"="*70}')

    # Download references
    print('\n[1/3] Downloading references...')
    for name, acc in PANEL.items():
        ref = download_ref(acc)
        status = '✓' if ref else '✗'
        size = os.path.getsize(ref) if ref else 0
        print(f'  {status} {name} ({acc}) — {size:,} bp')

    # Analyze samples
    print(f'\n[2/3] Analyzing samples...')
    all_results = {}

    for sample_name, reads_file, metadata in SAMPLES:
        print(f'\n  --- {sample_name} ({metadata}) ---')
        sample_results = {}

        for pathogen, accession in PANEL.items():
            ref_file = f'{REFS_DIR}/{accession}.fasta'
            if not os.path.exists(ref_file):
                continue

            outprefix = f'{WORK_DIR}/{sample_name}_vs_{pathogen}'
            try:
                result = map_sample_to_pathogen(reads_file, ref_file, outprefix)
                sample_results[pathogen] = result

                if result['mapped_reads'] > 100:
                    print(f'    🔴 {pathogen}: {result["mapping_rate"]}% rate, '
                          f'{result["avg_depth"]}x depth, '
                          f'{result["genome_coverage"]}% coverage')
                elif result['mapped_reads'] > 0:
                    print(f'    🟡 {pathogen}: {result["mapped_reads"]} reads '
                          f'({result["mapping_rate"]}%)')
            except Exception as e:
                print(f'    ✗ {pathogen}: ERROR — {e}')

        all_results[sample_name] = {
            'metadata': metadata,
            'results': sample_results,
        }

    # Generate report
    print(f'\n[3/3] Generating report...\n')
    print(f'{"="*70}')
    print(f'  PATHOGEN IDENTIFICATION REPORT')
    print(f'{"="*70}')

    positive_hits = defaultdict(list)
    for sample_name, data in all_results.items():
        for pathogen, result in data['results'].items():
            if result['mapping_rate'] > 1 and result['avg_depth'] > 3 and result['genome_coverage'] > 10:
                positive_hits[pathogen].append({
                    'sample': sample_name,
                    'type': data['metadata'],
                    'mapping_rate': result['mapping_rate'],
                    'avg_depth': result['avg_depth'],
                    'genome_coverage': result['genome_coverage'],
                })

    if positive_hits:
        print(f'\n  SUSPECT PATHOGENS DETECTED:')
        print(f'  {"-"*65}')
        for pathogen, hits in sorted(positive_hits.items(),
                                      key=lambda x: -max(h['avg_depth'] for h in x[1])):
            n_clinical = sum(1 for h in hits if h['type'] == 'clinical_stool')
            n_food = sum(1 for h in hits if h['type'] == 'food_environmental')
            avg_depth = sum(h['avg_depth'] for h in hits) / len(hits)
            avg_cov = sum(h['genome_coverage'] for h in hits) / len(hits)

            print(f'\n  >>> {pathogen}')
            print(f'      Positive: {len(hits)}/{len(all_results)} samples '
                  f'(clinical: {n_clinical}, food/env: {n_food})')
            print(f'      Avg depth: {avg_depth:.1f}x  |  Avg coverage: {avg_cov:.1f}%')
            for h in hits:
                print(f'        - {h["sample"]} ({h["type"]}): '
                      f'{h["mapping_rate"]}%, {h["avg_depth"]}x, {h["genome_coverage"]}%')

        # Conclusion
        top = max(positive_hits.items(), key=lambda x: len(x[1]))
        print(f'\n{"="*70}')
        print(f'  CONCLUSION')
        print(f'{"="*70}')
        print(f'  Most likely causative agent: {top[0]}')
        print(f'  Detected in {len(top[1])} of {len(all_results)} samples')
        in_food = sum(1 for h in top[1] if h['type'] == 'food_environmental')
        if in_food > 0:
            print(f'  Found in BOTH clinical AND food/environmental samples')
            print(f'  => Foodborne transmission CONFIRMED')
    else:
        print('  No significant pathogen hits detected.')

    # Save results
    report_file = '/root/biolab/work/outbreak_report.json'
    with open(report_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f'\n  Full results saved to: {report_file}')
    print(f'{"="*70}')


if __name__ == '__main__':
    main()
