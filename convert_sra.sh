#!/bin/bash
set -e
cd /tmp/outbreak_fastq
for acc in SRR35602535 SRR35602538 SRR35602539 SRR35602540; do
    sra_path="/mnt/c/Users/rdpuser/bacterial_assemble/outbreak_data/${acc}.sra"
    if [ ! -f "$sra_path" ]; then
        echo "MISSING: $sra_path"
        continue
    fi
    echo "=== Converting $acc ==="
    fasterq-dump --split-3 --threads 8 --bufsize 1G --curcache 2G -O . "$sra_path" 2>&1 | tail -5
    gzip -f ${acc}*.fastq 2>/dev/null || true
    ls -lh ${acc}*.fastq.gz 2>/dev/null
    echo "Done: $acc"
    echo ""
done
echo "=== ALL DONE ==="
