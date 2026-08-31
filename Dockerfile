# ==============================================================================
# Multi-stage build. pip install and the large reference-data downloads run
# in isolated stages so their memory usage doesn't stack cumulatively within
# one continuous build process.
#
# pyensembl's genome build (GTF parsing/indexing) is deliberately NOT done
# here at build time - it reliably triggered OOMKilled (exit 137) on Hugging
# Face's build environment, which runs with less memory than the runtime
# hardware tier you select. It's built once at container startup instead,
# in a background thread in app.py, using the real runtime hardware. See
# README for details on this tradeoff.
# ==============================================================================

# ---------- Stage: deps (install Python packages) ---------------------------
FROM python:3.11-slim AS deps

RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN useradd -m -u 1000 appuser
USER appuser
ENV HOME=/home/appuser \
    PATH=/home/appuser/.local/bin:$PATH

WORKDIR /home/appuser/app
COPY --chown=appuser:appuser requirements.txt .
RUN pip install --no-cache-dir --user -r requirements.txt

# ---------- Stage: refdata (download large reference data files) -----------
FROM python:3.11-slim AS refdata

RUN apt-get update && apt-get install -y --no-install-recommends curl \
    && rm -rf /var/lib/apt/lists/*
RUN useradd -m -u 1000 appuser
USER appuser
WORKDIR /home/appuser/app

# ClinVar variant summary (~50MB from NCBI's FTP)
RUN curl -fsSL -o clinvar_variant_summary.txt.gz \
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

# UniProt's human ID-mapping reference file (~142MB uncompressed - too large
# for git), decompressed to the exact filename core.py expects.
RUN curl -fsSL -o HUMAN_9606_idmapping.dat.gz \
    https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz \
    && gunzip HUMAN_9606_idmapping.dat.gz

# ---------- Final stage: assemble the runtime image -------------------------
FROM python:3.11-slim

RUN useradd -m -u 1000 appuser
USER appuser
ENV HOME=/home/appuser \
    PATH=/home/appuser/.local/bin:$PATH \
    PYENSEMBL_CACHE_DIR=/home/appuser/.cache/pyensembl

WORKDIR /home/appuser/app

COPY --from=deps --chown=appuser:appuser /home/appuser/.local /home/appuser/.local
COPY --from=refdata --chown=appuser:appuser /home/appuser/app/clinvar_variant_summary.txt.gz .
COPY --from=refdata --chown=appuser:appuser /home/appuser/app/HUMAN_9606_idmapping.dat .

# App code + your data files (mane_transcripts.csv,
# postsynaptic_genes_aiupred_score_list.json, and any pre-fetched gene .pkl
# files you want to ship with the image - see README)
COPY --chown=appuser:appuser . .

EXPOSE 7860
# Bind to $PORT if the platform provides one (Render, etc.), otherwise fall
# back to 7860 (matches local `docker run -p 7860:7860` and Hugging Face
# Spaces' fixed port convention). Shell form (not exec array form) so the
# ${PORT:-7860} expansion actually happens.
CMD gunicorn --bind 0.0.0.0:${PORT:-7860} --workers 2 --timeout 300 app:app