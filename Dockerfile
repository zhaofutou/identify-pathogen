FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Install all system dependencies in one layer
RUN apt-get update -qq && apt-get install -y -qq --no-install-recommends \
    # Bioinformatics tools
    minimap2 \
    samtools \
    bedtools \
    kraken2 \
    fastqc \
    spades \
    prokka \
    # Python
    python3 \
    python3-pip \
    # Network tools for reference downloads
    curl \
    wget \
    procps \
    # Cleanup
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip3 install --no-cache-dir flask werkzeug quast

# Create directory structure
# /data/refs     — reference genomes (persistent volume)
# /data/kraken2  — Kraken2 database (mount from host, ~8GB)
# /data/results  — analysis results (persistent volume)
# /data/uploads  — uploaded files (persistent volume)
# /data/work     — temp work files
RUN mkdir -p /app /data/uploads /data/refs /data/results /data/work /data/kraken2

# Copy application code
COPY webapp/app_docker.py /app/webapp/app.py
COPY webapp/templates/ /app/webapp/templates/
COPY pathogen_panel.py /app/
COPY entrypoint.sh /app/
RUN chmod +x /app/entrypoint.sh

# Set working directory
WORKDIR /app

# Expose web port
EXPOSE 5050

# Health check
HEALTHCHECK --interval=30s --timeout=5s --retries=3 \
    CMD curl -f http://localhost:5050/ || exit 1

# Default: run the web app
ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["web"]
