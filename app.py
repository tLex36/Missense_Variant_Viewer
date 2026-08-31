"""
app.py — Flask web app for the Roche Lab Variant Viewer.

Replaces the original Streamlit UI with:
  - a landing page: gene name box -> transcript dropdown -> Generate button
  - a tabbed results view (Missense / Truncations), each with a Plotly chart,
    a data table, and the 3D structure viewer

All actual analysis logic lives in core.py / gene_pipeline.py / structure_viewer.py,
extracted unchanged from the original Streamlit script.
"""
import logging
import os
import subprocess
import threading
import traceback

from flask import Flask, render_template, request, jsonify

import core
from gene_pipeline import get_or_build_gene_dict
from structure_viewer import build_structure_html

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("variant_viewer")

app = Flask(__name__)

# Fail fast and clearly if the templates/static folders aren't where Flask expects
# them, rather than surfacing a mid-request TemplateNotFound stack trace.
_here = os.path.dirname(os.path.abspath(__file__))
_missing = [
    p for p in ("templates/index.html", "static/style.css", "static/app.js")
    if not os.path.exists(os.path.join(_here, p))
]
if _missing:
    raise RuntimeError(
        "Missing expected file(s): " + ", ".join(_missing) + ". "
        f"These must live in templates/ and static/ subfolders next to app.py "
        f"(app.py is at {_here}). Nothing in this app writes these files - "
        f"they need to already be in place before startup."
    )

# ---------------------------------------------------------------------------
# pyensembl's reference genome (GTF parsing + indexing) is too memory-heavy
# to build reliably inside Hugging Face's Docker build environment (it was
# previously baked in at build time, but that step reliably hit an OOMKilled
# build error - see README). It's built here instead, once, at container
# startup, on the real runtime hardware. This runs in a background thread so
# the app starts serving HTTP immediately rather than blocking container
# startup (avoiding any platform health-check timeout); gene-lookup routes
# return 503 with a clear message until it's done. pyensembl's own install
# command is idempotent, so this is a fast no-op on any restart where the
# data already exists (e.g. persistent storage is attached).
# ---------------------------------------------------------------------------
_pyensembl_ready = threading.Event()


def _prepare_pyensembl():
    try:
        logger.info(
            "Ensuring pyensembl reference data is installed - this can take "
            "several minutes on a fresh container with no persistent storage."
        )
        subprocess.run(
            ["pyensembl", "install", "--release", "112", "--species", "human"],
            check=True,
        )
        logger.info("pyensembl reference data ready.")
    except Exception:
        logger.exception("Failed to prepare pyensembl reference data")
    finally:
        _pyensembl_ready.set()


threading.Thread(target=_prepare_pyensembl, daemon=True).start()


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/transcripts")
def api_transcripts():
    if not _pyensembl_ready.is_set():
        return jsonify({
            "error": "The server is still preparing reference genome data "
                     "(first startup can take several minutes) - please try again shortly."
        }), 503

    gene = request.args.get("gene", "").strip().upper()
    if not gene:
        return jsonify({"error": "No gene provided"}), 400

    try:
        transcript_options = core.get_transcripts(gene)
    except Exception as e:
        logger.exception("Failed to fetch transcripts for %s", gene)
        return jsonify({"error": f"Couldn't find gene '{gene}': {e}"}), 404

    if not transcript_options:
        return jsonify({"error": f"No protein-coding transcripts found for '{gene}'"}), 404

    try:
        mane_id = core.get_mane_transcript(gene)
    except Exception:
        mane_id = None
    try:
        canonical_id = core.get_uniprot_canonical_transcript(gene)
    except Exception:
        canonical_id = None

    transcript_options = sorted(
        transcript_options,
        key=lambda t: (0 if t[0] == mane_id else 1 if t[0] == canonical_id else 2)
    )

    results = []
    for tid, aa_len in transcript_options:
        if tid == mane_id and tid == canonical_id:
            label = f"{tid}  ({aa_len} aa)  — MANE/Canonical Transcript"
        elif tid == mane_id:
            label = f"{tid}  ({aa_len} aa)  — MANE Transcript"
        elif tid == canonical_id:
            label = f"{tid}  ({aa_len} aa)  — Canonical Transcript"
        else:
            label = f"{tid}  ({aa_len} aa)"
        results.append({"id": tid, "label": label})

    return jsonify({"gene": gene, "transcripts": results})


@app.route("/api/generate", methods=["POST"])
def api_generate():
    if not _pyensembl_ready.is_set():
        return jsonify({
            "error": "The server is still preparing reference genome data "
                     "(first startup can take several minutes) - please try again shortly."
        }), 503

    data = request.get_json(force=True)
    gene = (data.get("gene") or "").strip().upper()
    transcript_id = data.get("transcript_id")
    window_size = int(data.get("window_size", 5))
    normalize_missense = bool(data.get("normalize_missense", False))
    flip_gnomad = bool(data.get("flip_gnomad", False))
    intron_size = int(data.get("intron_size", 200))

    if not gene or not transcript_id:
        return jsonify({"error": "gene and transcript_id are required"}), 400

    try:
        gene_dict = get_or_build_gene_dict(gene, transcript_id)

        df_trunc = core.create_truncation_df(
            gene_dict['gnomad_truncation_variants'], gene_dict['clinvar_truncation_variants'],
            gene_dict['gnomad_genomic_variants'], transcript_id,
            gene_dict['exon_intron_map'], gene_dict['cds_regions']
        )
        df_missense = core.create_missense_df(
            gene_dict['clinvar_missense_VEP_annotations'], gene_dict['clinvar_missense_variants'],
            gene_dict['gnomad_missense_VEP_annotations'], gene_dict['gnomad_missense_variants'],
            gene_dict['gnomad_genomic_variants'], transcript_id
        )

        last_exon = core.get_last_coding_exon_number(
            transcript_id, gene_dict['exon_intron_map'], gene_dict['cds_regions']
        )

        def highlight_last_exon_cell(val):
            return "background-color: yellow" if val == last_exon else ""

        fig_missense = core.plot_missense_cluster_chart(
            gene_dict['clinvar_missense_variants'], gene_dict['gnomad_missense_variants'],
            df_missense, gene_dict['aiupred_binding'], gene_dict['aiupred_disorder'],
            gene_dict['sequence'], window_size, normalize_missense, flip_gnomad
        )
        fig_trunc = core.plot_exon_map_plotly(
            gene_dict['exon_intron_map'], gene_dict['clinvar_truncation_variants'],
            gene_dict['gnomad_truncation_variants'], df_trunc, cds_map=gene_dict['cds_regions'],
            intron_size=intron_size, gene_name=gene, transcript_id=transcript_id
        )

        structure_html = build_structure_html(gene_dict, df_missense, df_trunc, transcript_id)

        missense_table_html = df_missense.style.applymap(
            core.highlight_by_database, subset=['database']
        ).to_html()
        trunc_table_html = df_trunc.style.applymap(
            highlight_last_exon_cell, subset=['exon']
        ).applymap(core.highlight_by_database, subset=['database']).to_html()

        return jsonify({
            "gene": gene,
            "transcript_id": transcript_id,
            "length": gene_dict.get("length"),
            "missense_chart_html": fig_missense.to_html(
                full_html=False, include_plotlyjs=False, config={"scrollZoom": True}
            ),
            "missense_table_html": missense_table_html,
            "trunc_chart_html": fig_trunc.to_html(
                full_html=False, include_plotlyjs=False, config={"scrollZoom": True}
            ),
            "trunc_table_html": trunc_table_html,
            "structure_html": structure_html,
        })

    except Exception as e:
        logger.exception("Failed to generate analysis for %s / %s", gene, transcript_id)
        return jsonify({"error": str(e), "trace": traceback.format_exc()}), 500


if __name__ == "__main__":
    # Hugging Face Docker Spaces routes traffic to port 7860 by default
    app.run(host="0.0.0.0", port=7860, debug=False)