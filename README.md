---
title: Missense Variant Viewer
emoji: 🧬
colorFrom: blue
colorTo: indigo
sdk: docker
app_port: 7860
pinned: false
---

# Roche Lab Variant Viewer — Flask version

**[Live demo on Hugging Face Spaces](https://huggingface.co/spaces/YOUR-USERNAME/YOUR-SPACE-NAME)**
*(update this link once your Space is live — see the deployment steps at the bottom of this file)*

This is a refactor of `web_app_06_18_26.py` (the original Streamlit app) into a
Flask backend + plain HTML/CSS/JS frontend, structured for deployment as a
Hugging Face Docker Space.

## What changed vs. the original

All analysis logic is unchanged — extracted, not rewritten, from your original
script:

| File | What it contains |
|---|---|
| `core.py` | Every pure data-fetching / consequence-calculation / dataframe / Plotly chart function (lines 1–1642 of the original) |
| `gene_pipeline.py` | The 3-tier caching logic (`get_or_build_gene_dict`) that decides whether to load from a `.pkl`, partially reuse one, or fetch everything fresh |
| `structure_viewer.py` | The full 3Dmol.js structure-viewer HTML/JS builder (variant density, binding/disorder, hydrophobicity, charge coloring, surface mesh controls — all of it) |
| `app.py` | Flask routes wiring it together: `/`, `/api/transcripts`, `/api/generate` |
| `templates/index.html` + `static/` | The new frontend: gene box → transcript dropdown → Generate → tabbed results (Missense / Truncations), each with its Plotly chart + table, plus a shared structure viewer below |

One design note: the original app rendered **one** structure viewer below both
tabs (its coloring options pull from missense *and* truncation data together),
so I kept that — the structure panel stays visible under whichever tab is
active, rather than duplicating it per tab.

**Security fix applied:** the hardcoded NCBI API key from the original script
is gone — it wasn't actually referenced anywhere in the code (dead constant),
so `core.py` now reads it from an `NCBI_API_KEY` environment variable if you
ever need it.

## Files you still need to supply

I don't have your actual data files, so these aren't included — copy them
into the project root (next to `app.py`) before testing:

- `HUMAN_9606_idmapping.dat` — **don't copy this one manually**; it's 141MB
  uncompressed (too large for git) and downloaded fresh instead. Get it with:
  ```bash
  curl -fsSL -o HUMAN_9606_idmapping.dat.gz \
      https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
  gunzip HUMAN_9606_idmapping.dat.gz
  ```
  (The Dockerfile does this same download automatically at build time, so you
  only need this command for local testing.)
- `mane_transcripts.csv`
- `postsynaptic_genes_aiupred_score_list.json` — if this file is over 10MB and
  you're pushing to Hugging Face (their per-file limit is 10MiB, stricter
  than GitHub's 100MB), track it with Git LFS rather than committing it
  directly:
  ```bash
  git lfs install
  git lfs track "postsynaptic_genes_aiupred_score_list.json"
  git add .gitattributes postsynaptic_genes_aiupred_score_list.json
  ```
- Any existing `{GENE}-{TRANSCRIPT}_sheet.pkl` cache files you want to reuse
  locally. These are excluded from git by default (`.gitignore`) since
  they're regenerated automatically and can individually exceed 10MB; if you
  want specific genes preloaded in a deployment, `git lfs track` those
  particular `.pkl` files the same way as above instead of removing the
  `.gitignore` rule wholesale.

## Testing locally

**1. Install dependencies** (this pulls in pyensembl, biopython, torch, etc. —
expect this to take a few minutes and a few GB of disk space):

```bash
cd webapp
python3 -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

**2. Build the pyensembl reference database** (one-time, several GB, several
minutes — this is the step the Dockerfile bakes into the image at build time):

```bash
pyensembl install --release 112 --species human
```

Note: `pyensembl` keeps an internal table mapping each Ensembl release number
to its reference assembly, and that table only covers releases known as of
whatever version you have installed. If you ever see `ValueError: No genome
for homo_sapiens in Ensembl release ###`, it means your installed pyensembl
predates that release — run `pip install --upgrade pyensembl` and retry.

Similarly, `varcode` (used for variant-effect prediction) needs to keep pace
with how newer Ensembl releases represent things — an old varcode can throw
errors like `module 'varcode.effects.effect_classes' has no attribute
'...'` against release 112+ data. Same fix: `pip install --upgrade varcode`.

**After any one-off `pip install --upgrade <package>` command**, re-run
`pip install -r requirements.txt` afterward to make sure the rest of the
pinned stack (numpy in particular) hasn't drifted out of sync — an isolated
`--upgrade` only resolves the package you named, not everything already
installed, which is exactly what caused the numpy/varcode conflict earlier
in this same setup process.

(Note: unlike many CLIs, pyensembl doesn't support `python -m pyensembl` — it
only exposes a console-script entry point, so it must be run as the bare
`pyensembl` command. If that gives "command not found", your venv's `bin/`
isn't on PATH, or the earlier `pip install` didn't actually succeed — check
with `pip show pyensembl`.)

**3. Copy your data files** (`mane_transcripts.csv`,
`postsynaptic_genes_aiupred_score_list.json`, any `.pkl` caches) into the
`webapp/` folder, next to `app.py`. Download `HUMAN_9606_idmapping.dat`
separately using the command in "Files you still need to supply" above.

**4. Run it:**

```bash
python app.py
```

Then open `http://localhost:7860` in your browser. Type a gene name, tab or
click out of the field to trigger the transcript lookup, pick a transcript,
and click Generate. First run for a never-seen gene will take a while (same
1–5 minute caveat as the original Streamlit app); it'll be fast on repeat runs
since the `.pkl` cache is reused.

## Testing with Docker (matches the actual Hugging Face deployment)

```bash
docker build -t variant-viewer .
docker run -p 7860:7860 variant-viewer
```

Note: the Dockerfile's `COPY --chown=appuser:appuser . .` step picks up
whatever's in the `webapp/` folder at build time — so put your data files
there *before* building the image if you want them baked in.

## Deploying to Hugging Face Spaces

1. Create a new Space → SDK: **Docker** → visibility: your choice.
2. Push this folder's contents (including your data files, or a `.gitignore`
   for `*.pkl` if you'd rather they get generated fresh at runtime — see the
   caveat below) to the Space's git repo.
3. Hugging Face builds the Dockerfile automatically. Build time will be
   longer than usual because of the pyensembl reference-data step — that's
   expected.

**Important caveat about caching on Spaces:** the gene `.pkl` cache files are
written to the container's working directory, which is **ephemeral** on the
free tier — anything not baked into the image is wiped whenever the Space
restarts or redeploys (including its sleep/wake cycle after 48h of
inactivity). Two ways to handle this:

- Bake a handful of your most-used gene `.pkl` files directly into the image
  (just have them in the folder when you build/push) — good fit for the
  "3–5 preloaded gene" demo idea you mentioned earlier.
- Enable Hugging Face's **Persistent Storage** add-on (paid) and point the
  app's working directory at the mounted persistent volume, so every gene
  ever fetched stays cached across restarts.

**Important caveat about startup time:** pyensembl's reference-genome build
(parsing the full human GTF annotation) turned out to reliably exceed the
memory available in Hugging Face's *build* environment, which runs with a
tighter memory allowance than the runtime hardware you select — this
produced `OOMKilled` build failures even on CPU Basic's 16GB runtime tier.
Because of this, `app.py` builds it once at **container startup** instead
(in a background thread, so the app still responds to HTTP immediately) —
gene-lookup requests return a `503` with a clear "still preparing reference
data" message until that finishes, which can take several minutes on a
fresh container. Without persistent storage, this setup work repeats on
every restart, not just the very first deploy. If this becomes a real
annoyance, the Persistent Storage add-on mentioned above would let this
step be skipped on subsequent restarts too, the same way it helps with the
gene cache.

One more thing already handled: `torch` has been removed from both
`core.py`'s imports and `requirements.txt` — it was dead weight (only used
inside a commented-out AIUPred scoring class), and it was also the source of
a platform-specific install error on macOS (`torch==2.4.1+cpu` is a
Linux/Windows-only wheel). If you ever add live AIUPred scoring back in,
re-add plain `torch==2.4.1` (no `+cpu` suffix) for local/macOS development,
and keep the `+cpu` pinned version with the extra index URL for the Docker
image, since that one runs on Linux.

## Publishing this: GitHub + Hugging Face Spaces together

One codebase, pushed to two git remotes — GitHub for the source/portfolio,
Hugging Face for the running app.

**1. Before pushing anywhere**, check:
- `mane_transcripts.csv` and `postsynaptic_genes_aiupred_score_list.json` are
  present, and any gene `.pkl` caches you want to ship preloaded.
  `HUMAN_9606_idmapping.dat` should **not** be committed — it's too large for
  git (141MB) and is downloaded fresh by the Dockerfile at build time; your
  `.gitignore` already excludes it.
- Nothing sensitive is committed — `git status` and eyeball the diff before
  the first commit. `NCBI_API_KEY` already reads from an environment
  variable rather than being hardcoded (see `core.py`).
- Decide public vs. private for both the GitHub repo and the Space, if any
  of the underlying data or results aren't ready to be public.

**2. Push to GitHub:**
```bash
git init                      # skip if already a repo
git add .
git commit -m "Initial commit"
git remote add origin https://github.com/<your-username>/<repo-name>.git
git branch -M main
git push -u origin main
```

**3. Create the Space:** on huggingface.co, New Space → SDK: **Docker** →
name it → choose public/private. Hugging Face gives you a git URL like
`https://huggingface.co/spaces/<your-username>/<space-name>`.

**4. Push the same code to the Space**, as a second remote on the same repo:
```bash
git remote add space https://huggingface.co/spaces/<your-username>/<space-name>
git push space main
```
You'll be prompted for credentials — use your Hugging Face username and an
**access token** (Settings → Access Tokens, "write" scope) as the password,
not your account password.

**5. Watch it build:** open the Space's "Logs" tab. Given the pyensembl
reference-data and ClinVar-cache steps baked into the Dockerfile, expect the
first build to take a while (potentially 15+ minutes) — this is normal, not
a hang.

**6. Once it's live**, copy the Space's URL and update the link at the very
top of this README (and optionally add a badge):
```markdown
[![Hugging Face Spaces](https://img.shields.io/badge/%F0%9F%A4%97%20Hugging%20Face-Spaces-blue)](https://huggingface.co/spaces/<your-username>/<space-name>)
```
Commit and push that README update to `origin` (GitHub) — the Space itself
doesn't need the badge/link change, since it doesn't link to itself.

**Keeping both in sync later:** after any further changes, `git push origin
main` and `git push space main` separately (or `git push origin main &&
git push space main` in one line) to update both.