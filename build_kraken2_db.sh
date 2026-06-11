#!/bin/bash
set -e
DBDIR=/root/outbreak_investigation/kraken2_db
mkdir -p $DBDIR/taxonomy
cd $DBDIR/taxonomy

echo "=== Step 1: Download NCBI taxonomy ==="
# These files are small and essential
for f in names.dmp nodes.dmp merged.dmp delnodes.dmp; do
  if [ ! -f "$f" ]; then
    echo "Downloading $f..."
    curl -sL -o "${f}.gz" "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/${f}.gz" --connect-timeout 30 --max-time 120
    gunzip -f "${f}.gz" 2>/dev/null || true
  fi
done

# Get accession2taxid map
echo "Downloading nucl.accession2taxid..."
if [ ! -f "nucl.accession2taxid" ]; then
  curl -sL -o "nucl.accession2taxid.gz" "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/nucl.accession2taxid.gz" --connect-timeout 30 --max-time 600
  gunzip -f "nucl.accession2taxid.gz" 2>/dev/null || true
fi

echo "=== Taxonomy files ==="
ls -lh $DBDIR/taxonomy/

echo ""
echo "=== Step 2: Download bacterial sequences ==="
# Download a pre-built library if possible, or build from individual genomes
LIBDIR=$DBDIR/library
mkdir -p $LIBDIR/bacteria

# Try to download pre-built bacterial library from an alternative source
echo "Trying to get pre-built bacterial sequences..."
cd $LIBDIR/bacteria

# Download 16S rRNA sequences for quick pathogen ID (small, fast)
if [ ! -f "16S_rRNA.fasta" ]; then
  echo "Downloading 16S_rRNA references..."
  curl -sL -o "16S_rRNA.fasta.gz" \
    "https://ftp.ncbi.nlm.nih.gov/blast/db/16S_rRNA.tar.gz" \
    --connect-timeout 30 --max-time 600 2>/dev/null || true
fi

echo "=== Current status ==="
du -sh $DBDIR/
ls -la $DBDIR/taxonomy/
