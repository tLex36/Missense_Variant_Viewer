"""
core.py — refactored from the original Streamlit app (web_app_06_18_26.py).
Contains every pure data-fetching, consequence-calculation, dataframe-building,
and Plotly chart-building function, UNCHANGED in logic from the original script.
Only Streamlit-specific decorators/calls have been swapped for framework-agnostic
equivalents (functools.lru_cache instead of @st.cache_data, logging instead of st.write).
"""
import pickle
import requests
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import random
import xml.etree.ElementTree as ET
import json
import time
import zipfile
import io
import openpyxl
import pandas as pd
import re
import os
from scipy.stats import binom
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from Bio.PDB import MMCIFParser
from functools import lru_cache
import logging

logger = logging.getLogger("variant_viewer")

import sqlite3
_original_connect = sqlite3.connect
def _patched_connect(*args, **kwargs):
    kwargs['check_same_thread'] = False
    return _original_connect(*args, **kwargs)
sqlite3.connect = _patched_connect


from varcode import Variant
from pyensembl import EnsemblRelease


# SECURITY NOTE: the original script had a real NCBI API key hardcoded here.
# It isn't actually referenced anywhere else in the app (dead constant), so it's
# replaced with an env var read. If you do need it, set NCBI_API_KEY as an
# environment variable / HF Spaces secret - never commit it to the repo.
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
ENSEMBL_VEP_PAYLOAD_SIZE = 200

# @st.cache_resource
# def get_genome():
#     return EnsemblRelease(112)

import threading
_thread_local = threading.local()

def get_genome():
    if not hasattr(_thread_local, 'genome'):
        _thread_local.genome = EnsemblRelease(112)
    return _thread_local.genome

aa_map = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    # Add more if needed
}

@lru_cache(maxsize=1)
def load_transcript_to_uniprot():
    df = pd.read_csv("HUMAN_9606_idmapping.dat", sep="\t", header=None, names=["uniprot", "type", "value"])
    return dict(
        zip(df[df["type"] == "Ensembl_TRS"]["value"].str.split(".").str[0],
            df[df["type"] == "Ensembl_TRS"]["uniprot"].str.split("-").str[0])
    )

TRANSCRIPT_TO_UNIPROT = load_transcript_to_uniprot()

@lru_cache(maxsize=256)
def get_alphafold_cif(uniprot_id):
        # fetch metadata to get the correct URL
        api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        response = requests.get(api_url)
        response.raise_for_status()
        
        data = response.json()[0]
        cif_url = data["cifUrl"]  # correct URL from the API
        
        cif_response = requests.get(cif_url)
        cif_response.raise_for_status()
        return cif_response.text

#class AIUPredMPS(AIUPred):
#    def _setup_device(self, force_cpu, gpu_num):
#        if torch.backends.mps.is_available() and not force_cpu:
#            self.device = torch.device("mps")
#        elif torch.cuda.is_available() and not force_cpu:
#            self.device = torch.device(f"cuda:{gpu_num}")
#        else:
#            self.device = torch.device("cpu")
def get_mane_transcript(gene_symbol):
    mane_df = pd.read_csv("mane_transcripts.csv")
    row = mane_df[mane_df["symbol"] == gene_symbol]
    return row.iloc[0]["Ensembl_nuc"].split(".")[0]#[["Ensembl_nuc"]].to_dict("records")
def get_ensemblVEP_aa_consequences(variants):

    payload = create_ensemblVEP_payload(variants)
    results = []
    for batch in payload:
        batch_result = get_amino_acid_consequences_batch(batch)
        if batch_result is None:
            return None
        results.extend(batch_result)
    return results
    #Now each element of the list is 1 variant, with it's ensemblVEP data packed into it. ie. v[0] = all the VEP info
def create_ensemblVEP_payload(variants):
    payload_size = ENSEMBL_VEP_PAYLOAD_SIZE
    #Pass in the list of grch38_pos variants
    #Returns the list of payloads to iterate through in 200 variant chunks
    #convert the variant strings from #-######-ref-alt format to #:#####-#####/alt
    variants = [v.split("-")[0]+":"+v.split("-")[1]+"-"+v.split("-")[1]+"/"+v.split("-")[3] for v in variants]
    #split the variants into as many 200/payload_size chunked lists so that I can pass these onto the get_amino_acid_consequences_batch() function
    payload = [variants[i:i+payload_size]for i in range(0,len(variants),payload_size)]
    return payload

def get_amino_acid_consequences_batch(variants): #pass in a list of up to 200 variants and get that list back
    url = "https://rest.ensembl.org/vep/human/region"
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    payload = {"variants": variants, "LoF":1}  # list of "chrom:pos-pos/alt" strings
    response = requests.post(url, headers=headers, json=payload)
    try:
        response.raise_for_status()
        return response.json()
    except requests.HTTPError:
        print(f"EnsemblVEP HTTP error: {response.status_code}")
        return None
    
def parse_protein_consequence(name):
    aa_map = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    }
    match = re.search(r"p\.([A-Za-z]+)(\d+)([A-Za-z]+)", name)
    if match:
        ref = aa_map.get(match.group(1), match.group(1))
        loc = int(match.group(2))
        alt = aa_map.get(match.group(3), match.group(3))
        return ref, loc, alt
    return "N/A", "N/A", "N/A"

def my_parse_variants(variation):
    v = {}
    v["clinvar_id"]   = variation.get("VariationID", "N/A")
    v["name"]         = variation.get("VariationName", "N/A")
    v["accession"]    = variation.get("Accession", "N/A")

    interp = variation.find(".//Interpretation")
    if interp is not None:
        desc = interp.find("Description")
        v["clinical_significance"] = desc.text if desc is not None else "N/A"
    else:
        v["clinical_significance"] = "N/A"
    consequences = list({
        mc.get("Type", "")
        for mc in variation.findall(".//MolecularConsequence")
        if mc.get("Type")
    })
    v["molecular_consequences"] = consequences

    location = variation.find(
        ".//SimpleAllele/Location/SequenceLocation[@Assembly='GRCh38']"
    )
    if location is None:
        for loc in variation.findall(".//SequenceLocation[@Assembly='GRCh38']"):
            if loc.get("referenceAlleleVCF") is not None:
                location = loc
                break

    if location is not None:
        v["chromosome"]     = location.get("Chr", "N/A")
        v["position_start"] = location.get("start", "N/A")
        v["position_stop"]  = location.get("stop", "N/A")
        v["ref_allele"]     = location.get("referenceAlleleVCF", "N/A")
        v["alt_allele"]     = location.get("alternateAlleleVCF", "N/A")
        v['variant_id']     = "-".join([v['chromosome'],v['position_start'],v['ref_allele'],v['alt_allele']])
    else:
        v["chromosome"] = v["position_start"] = v["position_stop"] = "N/A"
        v["ref_allele"] = v["alt_allele"] = "N/A"

    #I want to extract from the name field the protein consequence of the missense variant and store that as a protein_location, protein_ref, and protein_alt
    #The name field is usually in the format "NM_000123.4:c.123A>G (p.Lys41Arg)" or similar. I can use a regex to extract the protein change information from the parentheses. The protein change is usually in the format "p.AAA123BBB" where AAA is the reference amino acid, 123 is the position, and BBB is the alternate amino acid. I can use a regex like r"\(p\.([A-Za-z]+)(\d+)([A-Za-z]+)\)" to extract these components.
    name = v["name"] or ""

    v["protein_ref"], v["protein_location"], v["protein_alt"] = parse_protein_consequence(name)
    review = variation.find(".//ReviewStatus")
    v["review_status"] = review.text if review is not None else "N/A"

    xref = variation.find(".//XRef[@DB='dbSNP']")
    v["rsid"] = f"rs{xref.get('ID')}" if xref is not None else "N/A"

    return v

def get_sequence_from_transcript(transcript_id):
    genome = EnsemblRelease(112)
    transcript = genome.transcript_by_id(transcript_id)
    protein_sequence = transcript.protein_sequence
    return protein_sequence

CLINVAR_SUMMARY_CACHE = "clinvar_variant_summary.txt.gz"
CLINVAR_SUMMARY_URL   = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"



def fetch_all_clinvar(gene):
    """
    Fetches all ClinVar variants for a gene from NCBI's pre-built FTP variant
    summary file. Downloads and caches the file on first run (~50 MB).
    Returns a list of dicts; each dict with GRCh38 VCF coordinates includes
    a 'variant_id' key in chr-pos-ref-alt format.
    """
    start = time.time()
    
    if not os.path.exists(CLINVAR_SUMMARY_CACHE):
        print("  Downloading ClinVar variant summary (one-time, ~50 MB)...")
        resp = requests.get(CLINVAR_SUMMARY_URL, stream=True, timeout=600)
        resp.raise_for_status()
        with open(CLINVAR_SUMMARY_CACHE, "wb") as f:
            for chunk in resp.iter_content(chunk_size=65536):
                f.write(chunk)
        print("  Download complete.")

    CLINVAR_COLS = [
        "GeneSymbol", "Assembly", "Chromosome", "PositionVCF",
        "ReferenceAlleleVCF", "AlternateAlleleVCF", "VariationID",
        "Name", "Start", "Stop", "ClinicalSignificance", "ReviewStatus"
    ]

    matching = []
    for chunk in pd.read_csv(
        CLINVAR_SUMMARY_CACHE,
        sep="\t",
        compression="gzip",
        chunksize=50000,
        usecols=lambda c: c in CLINVAR_COLS or "dbSNP" in c,
        low_memory=False,
        dtype=str,
    ):
        rows = chunk[(chunk["GeneSymbol"] == gene) & (chunk["Assembly"] == "GRCh38")]
        if len(rows):
            matching.append(rows)

    if not matching:
        print(f"  No GRCh38 ClinVar records found for {gene}")
        return []

    df = pd.concat(matching, ignore_index=True)

    has_vcf = (
        df["ReferenceAlleleVCF"].notna() &
        df["AlternateAlleleVCF"].notna() &
        ~df["ReferenceAlleleVCF"].isin(["na", "N/A", "-", "nan"]) &
        ~df["AlternateAlleleVCF"].isin(["na", "N/A", "-", "nan"])
    )
    df = df[has_vcf].copy()

    df["variant_id"] = (
        df["Chromosome"] + "-" +
        df["PositionVCF"]  + "-" +
        df["ReferenceAlleleVCF"] + "-" +
        df["AlternateAlleleVCF"]
    )

    rs_col = next((c for c in df.columns if "dbSNP" in c), None)

    variants = []
    for _, row in df.iterrows():
        rsid_raw = row[rs_col] if rs_col else ""
        variants.append({
            "variant_id":            row["variant_id"],
            "clinvar_id":            row.get("VariationID", "N/A"),
            "name":                  row.get("Name", "N/A"),
            "chromosome":            row["Chromosome"],
            "position_start":        row.get("Start", "N/A"),
            "position_stop":         row.get("Stop", "N/A"),
            "ref_allele":            row["ReferenceAlleleVCF"],
            "alt_allele":            row["AlternateAlleleVCF"],
            "clinical_significance": row.get("ClinicalSignificance", "N/A"),
            "review_status":         row.get("ReviewStatus", "N/A"),
            "rsid":                  f"rs{rsid_raw}" if rsid_raw not in ["", "nan", "-1"] else "N/A",
            "molecular_consequences": [],
        })

    print(f"  {gene}: {len(variants)} ClinVar variants with GRCh38 VCF coordinates")
    print(f"fetch_all_clinvar took {round(time.time() - start, 2)}s")
    return variants

def fetch_all_clinvar_for_gene(gene, api_key=None):
    """
    Fetch all ClinVar records for a gene. Paginates esearch to collect all IDs,
    then batch-fetches VCV records by ID directly (avoids WebEnv pagination bugs).
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    esearch_batch = 500
    fetch_batch   = 200
    all_ids = []
    total   = None
    retstart = 0

    # Step 1: paginate esearch to collect all matching IDs
    while True:
        search_params = {
            "db":      "clinvar",
            "term":    f"{gene}[gene] NOT \"copy number gain\"[Type] NOT \"copy number loss\"[Type]",
            "retmax":  esearch_batch,
            "retstart": retstart,
            "retmode": "json",
        }
        if api_key:
            search_params["api_key"] = api_key
        try:
            resp = requests.get(f"{base_url}/esearch.fcgi", params=search_params, timeout=30)
            resp.raise_for_status()
        except Exception as e:
            print(f"    Warning: ClinVar esearch failed for {gene}: {e}")
            return None, 0

        esearch = resp.json().get("esearchresult", {})
        if total is None:
            total = int(esearch.get("count", 0))
            print(f"    esearch total for {gene}: {total}")
        if total == 0:
            return None, 0

        ids = esearch.get("idlist", [])
        all_ids.extend(ids)
        if not ids or len(all_ids) >= total:
            break
        retstart += esearch_batch
        time.sleep(0.34 if not api_key else 0.11)

    print(f"    collected {len(all_ids)} IDs for {gene}")

    # Step 2: batch-fetch VCV records by direct ID
    all_variations = []
    for i in range(0, len(all_ids), fetch_batch):
        batch = all_ids[i:i + fetch_batch]
        time.sleep(0.34 if not api_key else 0.11)
        fetch_params = {
            "db":      "clinvar",
            "id":      ",".join(batch),
            "rettype": "vcv",
            "retmode": "xml",
        }
        if api_key:
            fetch_params["api_key"] = api_key
        try:
            fetch_resp = requests.get(f"{base_url}/efetch.fcgi", params=fetch_params, timeout=60)
            fetch_resp.raise_for_status()
            batch_root = ET.fromstring(fetch_resp.content)
            batch_vars = batch_root.findall(".//VariationArchive")
            all_variations.extend(batch_vars)
            print(f"    efetch IDs {i}–{i+len(batch)}: {len(batch_vars)} records")
        except Exception as e:
            print(f"    Warning: efetch failed for batch {i}–{i+len(batch)}: {e}")

    combined = ET.Element("ClinVarResult")
    combined.extend(all_variations)
    return combined, total

def fetch_all_gnomad(gene):
 
    #This function will query gnomAD for missense variants in the given gene symbol and return a list of missense variants with their positions and allele frequencies
    #I can use the gnomAD API to fetch this information. The endpoint for fetching variants by gene is "https://gnomad.broadinstitute.org/api/variants/search?query=gene:{gene_symbol}&variantType=missense"
    GNOMAD_API = "https://gnomad.broadinstitute.org/api"
    query = """
    query ($geneSymbol: String!) {
    gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
            variants(dataset: gnomad_r4) {
                variant_id
                pos
                ref
                alt
                consequence
                hgvsc
                hgvsp
                genome {
                    ac
                    an
                    homozygote_count
                    hemizygote_count
                }
                exome {
                    ac
                    an
                    homozygote_count
                    hemizygote_count
                }
            }
        }
    }
    """
    try:
        response = requests.post(GNOMAD_API, json={"query": query,"variables": {"geneSymbol": gene}})
        data = response.json()
        variants = data["data"]["gene"]["variants"]
        #gnomad_missense_list = [v for v in variants if "missense_variant" in v["consequence"]]
        return variants
    except requests.RequestException as e:
        print(f"Failed to fetch gnomAD missense variants for {gene}: {e}")
        return []
    
def calculate_consequences(variants_list, transcript_ids):
    genome = get_genome()
    gene_effects_dict = {}
    print(f"  calculate_consequences: processing {len(variants_list)} variants")
    for i in variants_list:
        parts = i.split("-")
        if len(parts) != 4:
            continue
        chrom, pos, ref, alt = parts
        if not pos.isdigit():
            continue
        if not all(set(s.upper()).issubset({'A','T','G','C'}) for s in (ref, alt)):
            continue
        try:
            v = Variant(contig=chrom, start=int(pos), ref=ref, alt=alt, genome=genome)

            effects_all = list(v.effects(raise_on_error=True))
            effects = [k for k in effects_all
                    if k.transcript_id is not None and k.transcript_id.split('.')[0] in transcript_ids]
            gene_effects_dict[i] = effects[0] if effects else None
            if i == variants_list[0]:  # just inspect the first variant
                print(f" First variant: {i}")
                print(f"  All effects: {[(type(k).__name__, k.transcript_id) for k in effects_all]}")
                print(f"  Matching effects: {effects}")
        except Exception as e:
            print(f"  Error on {i}: {e}")
            continue

    print(f"  calculate_consequences: done — {len(gene_effects_dict)} results")
    return gene_effects_dict

def filter_consequences_by_type(consequences_dict, types):
    types = [s.lower() for s in types]
    if "missense" in types:
        print("found missense")
        types.append("substitution")
    filtered_list = []
    for i in consequences_dict:
        if (type(consequences_dict[i]).__name__).lower() in types:
            filtered_list.append({i:consequences_dict[i]})
    return filtered_list

def fetch_gnomad_population_data(variants_list):
    #How to fetch the XX/XY data for each individual variant, by batching the query requests for the individual variant data
    #I pass in a list of variant_ids and it returns population breakdowns, including XX/XY
    results = {}
    batch_size = 20
    for i in range(0,len(variants_list),batch_size):
        batch = variants_list[i:i+batch_size] #Extract the variant ids from the list dictionaries of effects ie["X-123445-A-T": "effect"]
        data = fetch_gnomad_population_data_batch(batch)
        if data:
            results.update(data)
        else:   
            print("Empty response — likely rate limited")
        time.sleep(2)

    return results

def fetch_gnomad_population_data_batch(variant_ids, retries=3):
    aliases = ""
    for i, vid in enumerate(variant_ids):
        aliases += """
        v%d: variant(variantId: "%s", dataset: gnomad_r4) {
            genome { populations { id ac an } }
            exome { populations { id ac an } }
        }
        """ % (i, vid)
    
    query = "{ %s }" % aliases
    #print(query)
    for attempt in range(retries):
        try:
            response = requests.post(GNOMAD_API, json={"query": query}, timeout=30)
            response.raise_for_status()
            #print(response.status_code)
            #print(response.text)  # see what actually came back
            data = response.json()['data']
            remapped = {}
            for i, vid in enumerate(variant_ids):
                remapped[vid] = data.get(f"v{i}")
            return remapped
           
        except Exception as e:
            wait = 2 ** attempt  # 1s, 2s, 4s
            #print(f"Attempt {attempt+1} failed: {e}. Waiting {wait}s...")
            time.sleep(wait)
    
    return None
def get_exon_intron_coor(transcript_id):
    genome = get_genome()
    transcript = genome.transcript_by_id(transcript_id)
    exons = sorted(transcript.exons, key=lambda e: e.start)

    exon_map = [{"type": "exon", "start": e.start, "end": e.end} for e in exons]
    intron_map = [
        {"type": "intron", "start": exons[i].end + 1, "end": exons[i+1].start - 1}
        for i in range(len(exons) - 1)
    ]
    return sorted(exon_map + intron_map, key=lambda x: x["start"])
def get_cds_map(transcript_id):
    genome = get_genome()
    transcript = genome.transcript_by_id(transcript_id)
    cds_regions = [{"start": s, "end": e} for s, e in transcript.coding_sequence_position_ranges]
    return sorted(cds_regions, key=lambda x: x["start"])

def adjust_display_map(display_map,intron_size):
    """
    Pass in a exon/intron map, with a desired intron_size to get back an adjusted exon/intron map with shrunken introns
    display_map: exon/intron map formated as list of [{'start':int, 'end':int},...]
    returns the display_map with shrunken introns
    """
    adj_display_map = []
    offset = display_map[0]['start']
    for region in display_map:
        if region['type'] == "exon":
            exon_width = region['end']-region['start']
            adj_display_map.append({'type':'exon','start':offset,'end':offset+exon_width})
            #print(offset)
            offset += exon_width
        if region['type'] == "intron":
            adj_display_map.append({'type':'intron','start':offset,'end':offset+intron_size})
            offset += intron_size
    return adj_display_map

def adjust_trunc_variant_positions(x_coor, display_map, intron_size):
    """
    Pass in the truncation variants list, display_map, and intron_size.
    This will return the truncation variants list with an adjusted_start and adjusted_end value
    """
    adj_display_map = adjust_display_map(display_map,intron_size)
    adj_x_coor = []
    for x in x_coor:
        #print(f"X = {x}")
        i = max([int(j) for j in range(len(display_map)) if int(display_map[j]['start'])<=x])
        adj_start =adj_display_map[i]['start']+(x - display_map[i]['start'])
        adj_x_coor.append(adj_start)
    return adj_x_coor

def adjust_cds_map(cds_map, display_map, intron_size):

    adj_display_map = adjust_display_map(display_map,intron_size)
    offset = min([i['start'] for i in cds_map])
    for r in cds_map:
        start = r['start']
        end = r['end']
        i = max([j for j in range(len(display_map)) if display_map[j]['start']<=start])
        #print(f"determined index: {i} brings cds_map[start]:{start} closest to display_map[start]{display_map[i]['start']}")
        adj_start = adj_display_map[i]['start']+(start - display_map[i]['start'])
        adj_end = adj_display_map[i]['start']+(start - display_map[i]['start'])+(end-start)

        r['adj_start'] = adj_start
        r['adj_end'] = adj_end
    return cds_map

def find_exon_number(pos, region_map):
    exons = [r for r in region_map if r["type"] == "exon"]
    for i, exon in enumerate(exons, start=1):
        if int(exon["start"]) <= int(pos) <= int(exon["end"]):
            return i
    return None  # in intron or not found

def make_truncation_trace(variants, x_coor, y_coor, name, color='blue', visibility_arg = False):
    return go.Scatter(
        x=x_coor, y=y_coor,
        yaxis="y2",
        mode="markers",
        marker=dict(color=color, size=8),
        customdata=
           list(zip(variants['id'],variants['protein_consequence'],variants['consequence_terms'])),
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>"
            "consequence: %{customdata[1]}<br>"
            "protein_loc: %{customdata[2]}<br>"
            "<extra></extra>"
        ),
        name=name,
        visible = visibility_arg
    )

def make_missense_trace(variants, x_coor, y_coor, name, color='blue', visibility_arg = False, list_name = None):
    return go.Scatter(
        x=x_coor, y=y_coor,
        yaxis="y2",
        mode="markers",
        marker=dict(color=color, size=8),
        customdata=
           list(zip(variants['id'],variants['protein_consequence'],[list_name]*len(variants['database']),variants['polyphen prediction'],variants['sift prediction'])),
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>"
            "protein_loc: %{customdata[1]}<br>"
            "list_name: %{customdata[2]}<br>"
            "polyphen: %{customdata[3]}<br>"
            "sift: %{customdata[4]}<br>"
            "<extra></extra>"
        ),
        name=name,
        visible = visibility_arg
    )

def plot_exon_map_plotly(display_map, clinvar_truncating_variants = None, gnomad_truncating_variants = None, df_trunc = None, cds_map = None, intron_size = None, box_height=1,gene_name = None, transcript_id = None):
    """
    display_map: output of build_display_map()
    truncating_variants: list of dicts with keys 'display_x' and any variant data to show on hover
    """
    fig = go.Figure()

    # Draw exons and introns
    x_exon, y_exon = [], []

    if intron_size is None:
        x_exon, y_exon = [], []
        x_intron, y_intron = [], []

        for region in display_map:
            if region["type"] == "exon":
                x_exon += [region["start"], region["start"], region["end"], region["end"], None]
                y_exon += [0, box_height, box_height, 0, None]
        # All introns as one line trace
            if region["type"] == "intron":
                mid = box_height / 2
                x_intron += [region["start"], region["end"], None]
                y_intron += [mid, mid, None]

    if not intron_size is None:
        adj_display_map = adjust_display_map(display_map,intron_size)
        x_exon, y_exon = [], []
        x_intron, y_intron = [], []

        for region in adj_display_map:
            if region['type'] == 'exon':
                x_exon += [region["start"], region["start"], region["end"], region["end"], None]
                y_exon += [0, box_height, box_height, 0, None]
            if region["type"] == "intron":
                mid = box_height / 2
                x_intron += [region["start"], region["end"], None]
                y_intron += [mid, mid, None]
    
    x_cds, y_cds = [], []
    if not cds_map is None:
        if intron_size is None:
            for c in cds_map:
                x_cds += [c['start'], c['start'],c['end'],c['end'],None]
                y_cds += [0, box_height, box_height, 0, None]
        if not intron_size is None:
            adj_cds_map = adjust_cds_map(cds_map,display_map,intron_size)
            for c in adj_cds_map:
                x_cds += [c['adj_start'], c['adj_start'],c['adj_end'],c['adj_end'],None]
                y_cds += [0, box_height, box_height, 0, None]
    fig.add_trace(go.Scatter(
        x=x_exon, y=y_exon,
        fill="toself", #fillcolor="steelblue",
        fillcolor = "purple",
        line=dict(color="purple"),
        mode="lines", showlegend=False, hoverinfo="skip"
    ))
    fig.add_trace(go.Scatter(
        x=x_cds, y=y_cds,
        fill="toself", #fillcolor="steelblue",
        fillcolor = "steelblue",
        line=dict(color="steelblue"),
        mode="lines", showlegend=False, hoverinfo="skip"
    ))
    fig.add_trace(go.Scatter(
        x=x_intron, y=y_intron,
        mode="lines", line=dict(color="black", width=1),
        showlegend=False, hoverinfo="skip"
    ))
    is_sex_chrom = df_trunc['id'].iloc[0][0] in ['X','Y']
    clinvar_df = df_trunc[df_trunc['database'] == "ClinVar"]
    clinvar_x_coor = adjust_trunc_variant_positions(clinvar_df['location'],display_map,intron_size)
    clinvar_y_coor = [box_height/2 +random.uniform(0.3,0.8) for _ in range(len(clinvar_x_coor))]
    gnomad_df = df_trunc[df_trunc['database']=='GnomAD']
    gnomad_x_coor = adjust_trunc_variant_positions(gnomad_df['location'],display_map,intron_size)
    gnomad_y_coor = [box_height/2 +random.uniform(-0.2,0.3) for _ in range(len(gnomad_x_coor))]
    counted_gnomad_df = df_trunc[(df_trunc['database']=='GnomAD') & (df_trunc['total_count'] > 0)]
    counted_gnomad_x_coor = adjust_trunc_variant_positions(counted_gnomad_df['location'],display_map,intron_size)
    counted_gnomad_y_coor = [box_height/2 +random.uniform(-0.2,0.3) for _ in range(len(counted_gnomad_x_coor))]
    homo_df = df_trunc[(df_trunc['database']=='GnomAD') & (df_trunc['homozygote_count'] > 0)]
    homo_x_coor = adjust_trunc_variant_positions(homo_df['location'],display_map,intron_size)
    homo_y_coor = [box_height/2 +random.uniform(-0.2,0.3) for _ in range(len(homo_x_coor))]
    het_df = df_trunc[(df_trunc['database']=='GnomAD') & (df_trunc['heterozygote count'] > 0)]
    het_x_coor = adjust_trunc_variant_positions(het_df['location'],display_map,intron_size)
    het_y_coor = [box_height/2 +random.uniform(-0.2,0.3) for _ in range(len(het_x_coor))]

    fig.add_trace(make_truncation_trace(clinvar_df,clinvar_x_coor,clinvar_y_coor,"ClinVar truncating variants",color='red',visibility_arg=True))
    fig.add_trace(make_truncation_trace(gnomad_df, gnomad_x_coor, gnomad_y_coor, "All GnomAD truncating variants", color='blue',visibility_arg = True))
    fig.add_trace(make_truncation_trace(counted_gnomad_df, counted_gnomad_x_coor, counted_gnomad_y_coor, "Counted GnomAD truncating variants", color='blue'))
    fig.add_trace(make_truncation_trace(homo_df, homo_x_coor, homo_y_coor, "Homozygote GnomAD truncating variants", color='blue'))
    fig.add_trace(make_truncation_trace(het_df, het_x_coor, het_y_coor, "Heterozygote GnomAD truncating variants", color='blue'))
    
    if is_sex_chrom:
        #this would be XY count since only males can be hemizygous 
        xy_hemi_df = df_trunc[(df_trunc['database']=='GnomAD') & (df_trunc['hemizygote count'] > 0)]
        xy_hemi_x_coor = adjust_trunc_variant_positions(xy_hemi_df['location'],display_map,intron_size)
        xy_hemi_y_coor = [box_height/2 +random.uniform(-0.2,0.3) for _ in range(len(xy_hemi_x_coor))]
        fig.add_trace(make_truncation_trace(xy_hemi_df, xy_hemi_x_coor, xy_hemi_y_coor, "XY - Hemizygote truncating variants", color='blue'))

        #This would be XX count since only females can be homo or het for an X allele
        xx_homo_df = df_trunc[(df_trunc['database']=='GnomAD') & ((df_trunc['heterozygote count'] > 0) | (df_trunc['homozygote_count'] > 0))]
        xx_homo_het_x_coor = adjust_trunc_variant_positions(xx_homo_df['location'],display_map,intron_size)
        xx_homo_het_y_coor = [box_height/2 +random.uniform(-0.2,0.3) for _ in range(len(xx_homo_het_x_coor))]
        fig.add_trace(make_truncation_trace(df_trunc, xx_homo_het_x_coor, xx_homo_het_y_coor, "XX - Homo& Het truncating variants", color='blue'))

    if is_sex_chrom:
        n_existing = len(fig.data) - 7
    else:
        n_existing = len(fig.data) - 5
        
    x_min = (min(x_cds[::5]))-400
    x_max = (max(x_cds[2::5]))+400
    #strand = list(set([v['strand'] for v in gnomad_truncating_variants]))[0]

    ########Layout configuration!!!
    fig.update_layout(
    
        xaxis=dict(title="Position (bp from start of chromosome)", showgrid=False, range=[x_min, x_max]),
        yaxis=dict(visible=False, range=[-0.5, box_height + 0.5]),
        yaxis2=dict(
            overlaying="y",
            visible=False,
            range=[-0.5, box_height + 0.5]  # match primary y range
        ),
        height=400,
        #width=1400,
        plot_bgcolor="white"
    )
    if not is_sex_chrom:
        fig.update_layout(
        updatemenus=[dict(
            type="dropdown",
            x=0.0, y=1.15,
            xanchor="left",
            yanchor="top",
            active = 1,
            buttons=[
                dict(label="All GnomAD truncating variants",
                    method="update",
                    args=[{"visible": [True]*n_existing + [True,  True, False, False, False]}]),
                dict(label="Counted GnomAD truncating variants",
                    method="update",
                    args=[{"visible": [True]*n_existing + [True, False,  True, False, False]}]),
                dict(label="Homozygote GnomAD truncating variants",
                    method="update",
                    args=[{"visible": [True]*n_existing + [True, False, False, True, False]}]),
                dict(label="Heterozygote GnomAD truncating variants",
                    method="update",
                    args=[{"visible": [True]*n_existing + [True, False, False, False, True ]}]),
            ]
        )],margin=dict(t=80)  # add top margin to make room for buttons
    )
    if is_sex_chrom:
            fig.update_layout(
            updatemenus=[dict(
                type="dropdown",
                x=0.0, y=1.15,
                xanchor="left",
                yanchor="top",
                active = 1,
                buttons=[
                    dict(label="All GnomAD truncating variants",
                        method="update",
                        args=[{"visible": [True]*n_existing + [True,  True, False, False, False, False, False]}]),
                    dict(label="Counted GnomAD truncating variants",
                        method="update",
                        args=[{"visible": [True]*n_existing + [True, False,  True, False, False, False, False]}]),
                    dict(label="Homozygote GnomAD truncating variants",
                        method="update",
                        args=[{"visible": [True]*n_existing + [True, False, False, True, False, False, False]}]),
                    dict(label="Heterozygote GnomAD truncating variants",
                        method="update",
                        args=[{"visible": [True]*n_existing + [True, False, False, False, True, False, False]}]),
                    dict(label="XY - Hemizygote truncating variants",
                        method="update",
                        args=[{"visible": [True]*n_existing + [True, False, False, False, False, True, False]}]),
                    dict(label="XX - Homo& Het truncating variants",
                        method="update",
                        args=[{"visible": [True]*n_existing + [True, False, False, False, False, False, True ]}]),
                ]
            )],margin=dict(t=80)  # add top margin to make room for buttons
        )
    fig.update_layout(
        xaxis=dict(
            rangeslider=dict(visible=True),   # scrollable minimap below chart
            type="linear"
        )
    )
    fig.update_layout(
        xaxis=dict(showline=True, linewidth=2, linecolor="black", mirror=True),
        yaxis=dict(showline=True, linewidth=2, linecolor="black", mirror=True)
    )
    fig.update_layout(yaxis=dict(range=[-0.5, box_height + 3]))
         

    fig.update_layout(xaxis=dict(title="Position (bp from start of chromosome)", showgrid=False))

    fig.update_layout(dragmode="pan")
    fig.update_xaxes(
        range=[x_min, x_max],
        minallowed=x_min,
        maxallowed=x_max,
        autorange=False
    )
    title=dict(text=f"{gene_name}: Early truncations: {transcript_id}",
            font=dict(size=24),
            x=0.5,          # center horizontally
            xanchor="center")

    return fig
def is_missense(variation):
    for mc in variation.findall(".//MolecularConsequence"):
        so_term  = mc.get("Type", "").lower()
        function = mc.get("Function", "").lower()
        if "missense" in so_term or "missense" in function:
            return True
    # HGVS protein fallback
    for hgvs in variation.findall(".//HGVSExpression[@Type='hgvs, protein']"):
        text = hgvs.text or ""
        if text and not any(x in text for x in ["Ter", "*", "fs", "del", "ins", "="]):
            return True
    return False

def get_cleaned_missense_variants_data(VEP_consequences,transcript_id):
    results = []
    for v in VEP_consequences: 
        tc = v['transcript_consequences']
        tc = [i for i in tc if i['transcript_id'] == transcript_id][0]
        if any("missense" in s for s in tc['consequence_terms']):
            results.append(tc)
    return results

def get_uniprot_canonical_transcript(gene_symbol):
    # Search for the reviewed (Swiss-Prot) human entry
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "format": "json",
        "fields": "xref_ensembl,xref_refseq,cc_alternative_products"
    }
    response = requests.get(url, params=params)
    response.raise_for_status()
    results = response.json()["results"]
    
    if not results:
        return None

    entry = results[0]
    transcripts = []

    canonical_isoform = [isoform['isoformIds'][0] for isoform in entry['comments'][0]['isoforms'] if isoform['isoformSequenceStatus'] == "Displayed"][0] #Get the protein isoform ID from Uniprot that has a Displayed status
    canonical_transcript = [t['id'] for t in entry["uniProtKBCrossReferences"] if 'isoformId' in t and t['isoformId'] == canonical_isoform]#using the ID of the canonical isoform, I check which transcripts transcribe that ID and return it as a list
    canonical_transcript = [i.split(".")[0] for i in canonical_transcript if i.startswith("ENST")]#find the transcripts that start with ENST and return them dropping anything after the "." (which indicates the version number)
    
    return canonical_transcript[0]

def get_protein_sequences(gene_name, *transcript_ids):
    """
    Get amino acid sequences for a gene across multiple transcripts.

    Args:
        gene_name (str): Gene symbol for reference.
        *transcript_ids (str): Any number of Ensembl transcript IDs (ENST...).

    Returns:
        dict: Maps each transcript_id to its amino acid sequence string.
    """
    url = "https://rest.ensembl.org/sequence/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    results = {}
    for transcript_id in transcript_ids:
        response = requests.get(
            f"{url}/{transcript_id}",
            headers=headers,
            params={"type": "protein"}
        )
        response.raise_for_status()
        results[transcript_id] = response.json()["seq"]
    
    return results
def get_aiupred_scores(sequence):
    #This function will take a protein sequence as input and return the aiupred binding score and disorder score for each residue in the sequence
    predictor = AIUPredMPS()
    binding_scores = predictor.predict_binding(sequence).tolist()
    disorder_scores = predictor.predict_disorder(sequence).tolist()
    return binding_scores, disorder_scores
def fetch_snvs_for_gene(gene, api_key=None):
    """
    Fetch all SNV records from ClinVar for a given gene symbol.
    Returns parsed XML root, or None on failure.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    search_params = {
        "db": "clinvar",
        "term": f"{gene}[gene] AND \"single nucleotide variant\"[Type]",
        "retmax": 5000,
        "retmode": "json",
        "usehistory": "y",
    }
    if api_key:
        search_params["api_key"] = api_key

    try:
        search_resp = requests.get(
            f"{base_url}/esearch.fcgi",
            params=search_params,
            timeout=30
        )
        search_resp.raise_for_status()
        esearch = search_resp.json().get("esearchresult", {})

        total = int(esearch.get("count", 0))
        if total == 0:
            return None, 0

        time.sleep(0.34 if not api_key else 0.11)

        fetch_params = {
            "db": "clinvar",
            "query_key": esearch.get("querykey"),
            "WebEnv": esearch.get("webenv"),
            "rettype": "vcv",
            "retmode": "xml",
            "retmax": 5000,
        }
        if api_key:
            fetch_params["api_key"] = api_key

        fetch_resp = requests.get(
            f"{base_url}/efetch.fcgi",
            params=fetch_params,
            timeout=60
        )
        fetch_resp.raise_for_status()
        return ET.fromstring(fetch_resp.content), total

    except Exception as e:
        print(f"    Warning: ClinVar query failed for {gene}: {e}")
        return None, 0
def is_truncation(variation):
    truncating = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "nonsense"}

    # How to determine if an early truncation
    # HGVS protein fallback
    for hgvs in variation.findall(".//HGVSExpression[@Type='hgvs, protein']"):
        text = hgvs.text or ""
        if text and any(x in text for x in ["Ter", "*", "fs", "del", "ins"]):
            return True
    return False
def fetch_gnomad_missense_variants(gene_symbol):
    #This function will query gnomAD for missense variants in the given gene symbol and return a list of missense variants with their positions and allele frequencies
    #I can use the gnomAD API to fetch this information. The endpoint for fetching variants by gene is "https://gnomad.broadinstitute.org/api/variants/search?query=gene:{gene_symbol}&variantType=missense"
    GNOMAD_API = "https://gnomad.broadinstitute.org/api"
    query = """
    query ($geneSymbol: String!) {
    gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
            variants(dataset: gnomad_r4) {
                variant_id
                pos
                ref
                alt
                consequence
                hgvsc
                hgvsp
            }
        }
    }
    """
    try:
        response = requests.post(GNOMAD_API, json={"query": query,"variables": {"geneSymbol": gene_symbol}})
        data = response.json()
        variants = data["data"]["gene"]["variants"]
        gnomad_missense_list = [v for v in variants if "missense_variant" in v["consequence"]]
        return gnomad_missense_list
    except requests.RequestException as e:
        print(f"Failed to fetch gnomAD missense variants for {gene_symbol}: {e}")
        return []
    
def assemble_ensemblVEP_aa_consequences(total_ensemblVEP_consequences,*transcript_ids):
    results = []
    #I want a list of key:value pairs. The key is the transcript_ID and the pair is the missense variants
    for v in total_ensemblVEP_consequences: 
        tc = v['transcript_consequences']
        tc = [i for i in tc if i['transcript_id'] in transcript_ids]
        for t in tc:
            aa = t.get("amino_acids","N/A")
            if aa != "N/A":
                #print("processing amino acid: "+ str(aa))
                results.append({t['transcript_id']:str(aa[0])+str(t['protein_start'])+str(aa[-1])})
    
    tmp_results = {}
    for t in transcript_ids:
        tmp_results[t]  = ([r[t] for r in results if list(r.keys())[0] == t])
    results = tmp_results
    
    #print("successfully fetched all protein consequences for the transcripts: "+str([i for i in transcript_ids]))
    return results

def assemble_indexed_residue_list(residues,prot_len):
    residue_nums = []
    indexed_residues = [0]*prot_len
    for i in residues:
        if isinstance(i, int):
            residue_nums.append(i)
        else:
            match = int(re.findall(r'\d+',i)[0])
            residue_nums.append(match)   
    for i in range(prot_len):
        if i in residue_nums:
            indexed_residues[i] = 1
    return indexed_residues 

def convolve_v3(rsd_index_list, protein_len, window_size=5):
    n = window_size*2
    p = 0.5
    r_values = list(range(n+1))
    dist = [binom.pmf(r,n,p) for r in r_values ]
    # plt.bar(r_values, dist)
    # plt.show()
    norm = [float(i)/max(dist) for i in dist]
    
    convolved_list = []
    for i in range(protein_len):
        tmp = 0
        
        offset_numerator = max([(0 if (i-protein_len+int(len(dist)/2)) <0 else (i-protein_len+int(len(dist)/2))), 
                             (0 if -1*(i-int(len(dist)/2)) < 0 else -1*(i-int(len(dist)/2)))])
        #print(offset_numerator)
        
        for j in range(len(dist)):
            #if i - int(len(dist)/2) >0:
            #checks if the j values are within the range of the i values I guess *** go back and check this
            if i-int(len(dist)/2)+j >= 0 and i-int(len(dist)/2)+j < protein_len:
                tmp += rsd_index_list[i-int(len(dist)/2)+j]*dist[j]
            
        convolved_list.append(tmp*(1+(offset_numerator/(len(dist)/2))))
    return convolved_list

def plot_missense_cluster_chart(missense_clinvar, missense_gnomad, df_missense, binding_scores, disorder_scores, sequence, window_size = 5, normalize_missense = True, flip_gnomad = False):

    #Features I want to add to my chart:
    #   - rolling window slider
    #   - amino acid substitution matrix (checkbox)
    #   - adjust rolling window according to disorder value (checkbox)
    #Steps: Do I need to pre-generate all possible combinations of values (i.e. slider rolling window size, aa sub matrix, disorder adjusted rolling window size)
    DISORDER_THRESHOLD = 0.5
    length = len(sequence)
    box_height=1
    is_sex_chrom = df_missense['id'][0].split('-')[0] in ['X','Y']
    x = list(range(1, length + 1))

    clinvar_df = df_missense[df_missense['database'] == "ClinVar"]
    clinvar_x_coor = list(clinvar_df['position'])
    clinvar_y_coor = [-1*((i % 8) * (1 / 20) + 0.3) for i in range(len(clinvar_x_coor))]
    clinvar_scatter = make_missense_trace(clinvar_df, clinvar_x_coor,clinvar_y_coor, "ClinVar", color='red',visibility_arg=True,list_name = "clinvar")
    clinvar_density = convolve_v3(assemble_indexed_residue_list(clinvar_x_coor,length),length, window_size)

    gnomad_df = df_missense[df_missense['database'] == "GnomAD"]
    gnomad_x_coor = list(gnomad_df['position'])
    print(f"len gnomad_x_coor {len(gnomad_x_coor)}")
    gnomad_y_coor = [(-1 if flip_gnomad else 1)*((i % 8) * (1 / 20) + 0.3) for i in range(len(gnomad_x_coor))]
    gnomad_scatter = make_missense_trace(gnomad_df, gnomad_x_coor,gnomad_y_coor, "GnomAD", color='blue',visibility_arg='legendonly',list_name = "gnomad")
    gnomad_density = convolve_v3(assemble_indexed_residue_list(gnomad_x_coor,length),length, window_size)
   
    counted_gnomad_df = df_missense[(df_missense['database']=='GnomAD') & (df_missense['total_count'] > 0)]
    counted_gnomad_x_coor = list(counted_gnomad_df['position'])
    counted_gnomad_y_coor = [(-1 if flip_gnomad else 1)*((i % 8) * (1 / 20) + 0.3) for i in range(len(counted_gnomad_x_coor))]
    counted_gnomad_scatter = make_missense_trace(counted_gnomad_df, counted_gnomad_x_coor,counted_gnomad_y_coor, "Counted GnomAD", color='blue',list_name = "count_gnomad")
    counted_gnomad_density = convolve_v3(assemble_indexed_residue_list(counted_gnomad_x_coor,length),length, window_size)

    homo_df = df_missense[(df_missense['database']=='GnomAD') & (df_missense['homozygote_count'] > 0)]
    homo_x_coor = list(homo_df['position'])
    homo_y_coor = [(-1 if flip_gnomad else 1)*((i % 8) * (1 / 20) + 0.3) for i in range(len(homo_x_coor))]
    homo_scatter = make_missense_trace(homo_df, homo_x_coor,homo_y_coor, "Homozygote GnomAD", color='blue',list_name = "homo")
    homo_density = convolve_v3(assemble_indexed_residue_list(homo_x_coor,length),length, window_size)

    het_df = df_missense[(df_missense['database']=='GnomAD') & (df_missense['heterozygote count'] > 0)]
    het_x_coor = list(het_df['position'])
    #print(f"het_x_coor: {het_x_coor}")
    het_y_coor = [(-1 if flip_gnomad else 1)*((i % 8) * (1 / 20) + 0.3) for i in range(len(het_x_coor))]
    het_scatter = make_missense_trace(het_df, het_x_coor,het_y_coor, "Heterozygote GnomAD", color='blue',list_name = "het")
    het_density = convolve_v3(assemble_indexed_residue_list(het_x_coor,length),length,window_size)
    #print(f"het_density: {het_density}")

    if is_sex_chrom:
        #this would be XY count since only males can be hemizygous 
        xy_hemi_df = df_missense[(df_missense['database']=='GnomAD') & (df_missense['hemizygote count'] > 0)]
        xy_hemi_x_coor = list(xy_hemi_df['position'])
        xy_hemi_y_coor =  [(-1 if flip_gnomad else 1)*((i % 8) * (1 / 20) + 0.3) for i in range(len(xy_hemi_x_coor))]
        xy_hemi_scatter = make_missense_trace(xy_hemi_df, xy_hemi_x_coor, xy_hemi_y_coor, "XY - Hemizygote", color='blue',list_name = "hemi")
        xy_hemi_density = convolve_v3(assemble_indexed_residue_list(xy_hemi_x_coor,length),length,window_size)
        

    if normalize_missense:
        clinvar_density = [float(i)/max(clinvar_density) for i in clinvar_density]
        gnomad_density = [float(i)/max(gnomad_density) for i in gnomad_density]
        counted_gnomad_density = [float(i)/max(counted_gnomad_density) for i in counted_gnomad_density]
        homo_density = [float(i)/max(homo_density) for i in homo_density]
        het_density = [float(i)/max(het_density) for i in het_density]
        if is_sex_chrom:
            xy_hemi_density =[float(i)/max(xy_hemi_density) for i in xy_hemi_density]

    x = list(range(1, length + 1))
    y_diff_arr = np.array(np.subtract(gnomad_density, clinvar_density))
    counted_y_diff_arr =  np.array(np.subtract(counted_gnomad_density, clinvar_density))
    homo_y_diff_arr = np.array(np.subtract(homo_density, clinvar_density))
    het_y_diff_arr = np.array(np.subtract(het_density, clinvar_density))
    #print(f"het_y_diff_arr: {het_y_diff_arr}")
    if is_sex_chrom:
        xy_hemi_y_diff_arr = np.array(np.subtract(xy_hemi_density, clinvar_density))

    fill_clinvar_density_trace =go.Scatter(x=x, y=clinvar_density,fill="tozeroy", name = 'ClinVar densities', fillcolor="rgba(255,0,0,0.2)",
        line=dict(color="red"),showlegend=True, hoverinfo="skip", legendgroup='ClinVar densities', visible = 'legendonly')
    fill_gnomad_density_trace =go.Scatter(x=x, y=gnomad_density,fill="tozeroy", name = 'GnomAD densities', fillcolor="rgba(0,0,255,0.2)",
        line=dict(color="blue"),showlegend=True, hoverinfo="skip", legendgroup='GnomAD densities', visible = 'legendonly')
    fill_counted_gnomad_density_trace =go.Scatter(x=x, y=counted_gnomad_density,fill="tozeroy", name = 'Counted_GnomAD densities', fillcolor="rgba(0,0,255,0.2)",
        line=dict(color="blue"),showlegend=True, hoverinfo="skip", legendgroup='Counted GnomAD densities', visible = False)
    fill_homo_density_trace =go.Scatter(x=x, y=homo_density,fill="tozeroy", name = 'Homo_GnomAD densities', fillcolor="rgba(0,0,255,0.2)",
        line=dict(color="blue"),showlegend=True, hoverinfo="skip", legendgroup='Homo GnomAD densities', visible = False)
    fill_het_density_trace =go.Scatter(x=x, y=het_density,fill="tozeroy", name = 'Het_GnomAD densities', fillcolor="rgba(0,0,255,0.2)",
        line=dict(color="blue"),showlegend=True, hoverinfo="skip", legendgroup='Het GnomAD densities', visible = False)
    if is_sex_chrom:
        fill_xy_hemi_density_trace =go.Scatter(x=x, y=xy_hemi_density,fill="tozeroy", name = 'XY_Hemi_GnomAD densities', fillcolor="rgba(0,0,255,0.2)",
            line=dict(color="blue"),showlegend=True, hoverinfo="skip", legendgroup='Hemi GnomAD densities', visible = False)
        
    y_diff_trace = go.Scatter(x=x, y=list(y_diff_arr),mode="lines", name="tolerance",line=dict(color="purple"), hoverinfo="skip", legendgroup='tolerance')
    fill_blue_trace = go.Scatter(x=x, y=[v if v > 0 else 0 for v in y_diff_arr],fill="tozeroy", fillcolor="rgba(0,0,255,0.4)",line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip",legendgroup='tolerance')
    fill_red_trace =go.Scatter(x=x, y=[v if v < 0 else 0 for v in y_diff_arr],fill="tozeroy", fillcolor="rgba(255,0,0,0.4)",
        line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip", legendgroup='tolerance')
    
    counted_y_diff_trace = go.Scatter(x=x, y=list(counted_y_diff_arr),mode="lines", name="tolerance_counted",line=dict(color="purple"), hoverinfo="skip", legendgroup='tolerance_counted',visible=False)
    counted_fill_blue_trace = go.Scatter(x=x, y=[v if v > 0 else 0 for v in counted_y_diff_arr],fill="tozeroy", fillcolor="rgba(0,0,255,0.4)",line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip",legendgroup='tolerance_counted',visible=False)
    counted_fill_red_trace =go.Scatter(x=x, y=[v if v < 0 else 0 for v in counted_y_diff_arr],fill="tozeroy", fillcolor="rgba(255,0,0,0.4)",
        line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip", legendgroup='tolerance_counted',visible=False)   
    
    homo_y_diff_trace = go.Scatter(x=x, y=list(homo_y_diff_arr),mode="lines", name="tolerance_homo",line=dict(color="purple"), hoverinfo="skip", legendgroup='tolerance_homo',visible=False)
    homo_fill_blue_trace = go.Scatter(x=x, y=[v if v > 0 else 0 for v in homo_y_diff_arr],fill="tozeroy", fillcolor="rgba(0,0,255,0.4)",line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip",legendgroup='tolerance_homo',visible=False)
    homo_fill_red_trace =go.Scatter(x=x, y=[v if v < 0 else 0 for v in homo_y_diff_arr],fill="tozeroy", fillcolor="rgba(255,0,0,0.4)",
        line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip", legendgroup='tolerance_homo',visible=False)   
    
    het_y_diff_trace = go.Scatter(x=x, y=list(het_y_diff_arr),mode="lines", name="tolerance_het",line=dict(color="purple"), hoverinfo="skip", legendgroup='tolerance_het',visible=False)
    het_fill_blue_trace = go.Scatter(x=x, y=[v if v > 0 else 0 for v in het_y_diff_arr],fill="tozeroy", fillcolor="rgba(0,0,255,0.4)",line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip",legendgroup='tolerance_het',visible=False)
    het_fill_red_trace =go.Scatter(x=x, y=[v if v < 0 else 0 for v in het_y_diff_arr],fill="tozeroy", fillcolor="rgba(255,0,0,0.4)",
        line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip", legendgroup='tolerance_het',visible=False)
    
    if is_sex_chrom:
        xy_hemi_y_diff_trace = go.Scatter(x=x, y=list(xy_hemi_y_diff_arr),mode="lines", name="tolerance_hemi",line=dict(color="purple"), hoverinfo="skip", legendgroup='tolerance_hemi',visible=False)
        xy_hemi_fill_blue_trace = go.Scatter(x=x, y=[v if v > 0 else 0 for v in xy_hemi_y_diff_arr],fill="tozeroy", fillcolor="rgba(0,0,255,0.4)",line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip",legendgroup='tolerance_hemi',visible=False)
        xy_hemi_fill_red_trace =go.Scatter(x=x, y=[v if v < 0 else 0 for v in xy_hemi_y_diff_arr],fill="tozeroy", fillcolor="rgba(255,0,0,0.4)",
            line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip", legendgroup='tolerance_hemi',visible=False) 

#=============================================================================================================================================================================================================================================
    if binding_scores is None:
        binding_scores = [0]*len(x)
    if disorder_scores is None:
        disorder_scores = [0] * len(x)
    binding_score_trace = go.Scatter(x=x, y=list(binding_scores),mode="lines", name="AIUPred Binding",
        line=dict(color="green"), hoverinfo="skip", legendgroup='aiupred_binding', visible = 'legendonly')
    fill_binding_pink_trace = go.Scatter(x=x, y=[v if v >= DISORDER_THRESHOLD else 0 for v in binding_scores],fill="tozeroy", 
        fillcolor="rgba(255, 192, 203, 0.3)",line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip",legendgroup='aiupred_binding',visible = 'legendonly')
    
    disorder_score_trace = go.Scatter(x=x, y=list(disorder_scores),mode="lines", name="AIUPred Disorder",
        line=dict(color="yellow"), hoverinfo="skip", legendgroup='aiupred_disorder', visible = 'legendonly')
    fill_disorder_green_trace = go.Scatter(x=x, y=[v if v >= DISORDER_THRESHOLD else 0 for v in disorder_scores],fill="tozeroy", 
        fillcolor="rgba(0,128,0,0.3)",line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip",legendgroup='aiupred_disorder',visible = 'legendonly')

        #This would be XX count since only females can be homo or het for an X allele
        # xx_homo_df = df_trunc[(df_trunc['database']=='GnomAD') & ((df_trunc['heterozygote count'] > 0) | (df_trunc['homozygote_count'] > 0))]
        # xx_homo_het_x_coor = xx_homo_df['position']
        # xx_homo_het_y_coor = [-1*((i % 8) * (1 / 20) + 0.3) for i in range(len(het_x_coor))]
        #xx_homo_het_trace = make_truncation_trace(df_trunc, xx_homo_het_x_coor, xx_homo_het_y_coor, "XX - Homo & Het", color='blue')


    # ── Amino acid sequence heatmap track ─────────────────────────────────────────
    AA_GROUP = {
        'A':0,'V':0,'I':0,'L':0,'M':0,'F':0,'W':0,'P':0,
        'S':1,'T':1,'N':1,'Q':1,
        'R':2,'K':2,'H':2,
        'D':3,'E':3,
        'C':4,'G':4,'Y':4,
    }
    AA_COLORSCALE = [
        [0.0,   '#FFD700'], [0.125, '#FFD700'],
        [0.125, '#90EE90'], [0.375, '#90EE90'],
        [0.375, '#6495ED'], [0.625, '#6495ED'],
        [0.625, '#FF6347'], [0.875, '#FF6347'],
        [0.875, '#C0C0C0'], [1.0,   '#C0C0C0'],
    ]
    gnomad_by_pos  = gnomad_df.groupby('position')['alternate'].apply(list).to_dict()
    clinvar_by_pos = clinvar_df.groupby('position')['alternate'].apply(list).to_dict()
    max_gnomad  = max((len(v) for v in gnomad_by_pos.values()), default=0)
    max_clinvar = max((len(v) for v in clinvar_by_pos.values()), default=0)
    ref_row_idx = max_clinvar
    n_rows      = max_clinvar + 1 + max_gnomad
    z_matrix    = [[None] * length for _ in range(n_rows)]
    text_matrix = [[''  ] * length for _ in range(n_rows)]
    y_labels    = (
        [f'ClinVar {i+1}' for i in range(max_clinvar - 1, -1, -1)] +
        ['Ref'] +
        [f'GnomAD {i+1}' for i in range(max_gnomad)]
    )
    for col, aa in enumerate(sequence):
        z_matrix[ref_row_idx][col]    = AA_GROUP.get(aa, 4)
        text_matrix[ref_row_idx][col] = aa
    for pos, alts in gnomad_by_pos.items():
        col = int(pos) - 1
        for stack_i, aa in enumerate(alts):
            row = ref_row_idx + 1 + stack_i
            if 0 <= col < length and row < n_rows:
                z_matrix[row][col]    = AA_GROUP.get(aa, 4)
                text_matrix[row][col] = aa
    for pos, alts in clinvar_by_pos.items():
        col = int(pos) - 1
        for stack_i, aa in enumerate(alts):
            row = ref_row_idx - 1 - stack_i
            if 0 <= col < length and row >= 0:
                z_matrix[row][col]    = AA_GROUP.get(aa, 4)
                text_matrix[row][col] = aa
    seq_heatmap = go.Heatmap(
        z=z_matrix,
        x=list(range(1, length + 1)),
        y=y_labels,
        text=text_matrix,
        texttemplate="%{text}",
        textfont=dict(size=10),
        colorscale=AA_COLORSCALE,
        zmin=0, zmax=4,
        showscale=False,
        xgap=1, ygap=1,
        xaxis='x2',
        yaxis='y3',
        name='Residue Track',
        showlegend=True,
        visible='legendonly',
        hovertemplate='Pos: %{x}<br>AA: %{text}<extra></extra>',
    )

    fig = go.Figure()

    fig.add_trace(y_diff_trace) #true
    fig.add_trace(fill_blue_trace) #truee
    fig.add_trace(fill_red_trace) #true

    fig.add_trace(counted_y_diff_trace) #false
    fig.add_trace(counted_fill_blue_trace) #false
    fig.add_trace(counted_fill_red_trace) #false

    fig.add_trace(homo_y_diff_trace) #false
    fig.add_trace(homo_fill_blue_trace) #false
    fig.add_trace(homo_fill_red_trace) #false

    fig.add_trace(het_y_diff_trace) #false
    fig.add_trace(het_fill_blue_trace) #false
    fig.add_trace(het_fill_red_trace) #false
    
    if is_sex_chrom:
        fig.add_trace(xy_hemi_y_diff_trace) #false
        fig.add_trace(xy_hemi_fill_blue_trace) #false
        fig.add_trace(xy_hemi_fill_red_trace) #false

    # fig.add_trace(clinvar_missense_scatter)
    # fig.add_trace(gnomad_missense_scatter)
    fig.add_trace(fill_clinvar_density_trace) #legendonly
    fig.add_trace(fill_gnomad_density_trace) #legendonly -> false
    fig.add_trace(fill_counted_gnomad_density_trace) #False ->Legendonly
    fig.add_trace(fill_homo_density_trace)  #False ->Legendonly
    fig.add_trace(fill_het_density_trace)  #False ->Legendonly
    if is_sex_chrom:
        fig.add_trace(fill_xy_hemi_density_trace)  #False ->Legendonly

    fig.add_trace(binding_score_trace) #"legendonly"
    fig.add_trace(fill_binding_pink_trace) #"legendonly"
    fig.add_trace(disorder_score_trace) #"legendonly"
    fig.add_trace(fill_disorder_green_trace) #"legendonly"

    fig.add_trace(clinvar_scatter) #True * len(clinvar_x_coor)
    fig.add_trace(gnomad_scatter) #True -> False * len(gnomad_x_coor)
    fig.add_trace(counted_gnomad_scatter) #False -> True * len(counted_gnomad_x_coor)
    fig.add_trace(homo_scatter) #False -> True * len(homo_x_coor)
    fig.add_trace(het_scatter) #False -> True * len(het_x_coor)
    if is_sex_chrom:
        fig.add_trace(xy_hemi_scatter) #False -> True * len(xy_hemi_x_coor)
    fig.add_trace(seq_heatmap)
    button_layer_1_height = 1.08

    true3 = [True]*3
    false3 = [False]*3
    per_layout = ["legendonly"]*4
    if not is_sex_chrom:
        fig.update_layout(
            updatemenus=[
                dict(
                    active=1,
                    buttons=list([
                        dict(label="Total",
                            method="update",
                            args=[{"visible": (true3+ (false3 * 3) + ["legendonly"] + ["legendonly", False, False, False] + per_layout 
                                    + [True] + [True] +[False]
                                    + [False] + [False])},
                                {"title": "GnomAD",
                                    "annotations": []}]),
                        dict(label="Counted",
                            method="update",
                            args=[{"visible": (false3 + true3 + (false3 * 2) + ["legendonly"] + [False, "legendonly", False, False] + per_layout 
                                    + [True] + [False] +[True]
                                    + [False]+ [False])
                                    },
                                {"title": "Counted",
                                    "annotations": []}]),
                        dict(label="Homo",
                            method="update",
                            args=[{"visible": ((false3 * 2)+ true3 + (false3 * 1) + ["legendonly"] + [False, False, "legendonly", False] + per_layout 
                                    + [True] + [False]+[False] 
                                    + [True] + [False])
                                    },
                                {"title": "Homozygote",
                                    "annotations": []}]),
                        dict(label="Het",
                            method="update",
                            args=[{"visible": ((false3 * 3)+ true3 + (false3*0) + ["legendonly"] + [False, False, False , "legendonly"] + per_layout 
                                    + [True] + [False] +[False]
                                    + [False] + [True])
                                    },
                                {"title": "Het",
                                    "annotations": []}]),
                    ]),
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0,
                xanchor="left",
                y=button_layer_1_height,
                yanchor="top"
                )
            ])
    if is_sex_chrom:
        fig.update_layout(
            updatemenus=[
                dict(
                    active=1,
                    buttons=list([
                        dict(label="Total",
                            method="update",
                            args=[{"visible": (true3+ (false3 * 4) + ["legendonly"] + ["legendonly", False, False, False, False] + per_layout 
                                    + [True] + [False] +[False]
                                    + [False] + [False]+[False])},
                                {"title": "GnomAD",
                                    "annotations": []}]),
                        dict(label="Counted",
                            method="update",
                            args=[{"visible": (false3 + true3 + (false3 * 3) + ["legendonly"] + [False, "legendonly", False, False, False] + per_layout 
                                    + [True] + [False] +[True]
                                    + [False]+ [False]+[False])
                                    },
                                {"title": "Counted",
                                    "annotations": []}]),
                        dict(label="Homo",
                            method="update",
                            args=[{"visible": ((false3 * 2)+ true3 + (false3 * 2) + ["legendonly"] + [False, False, "legendonly", False, False] + per_layout 
                                    + [True] + [False]+[False] 
                                    + [True] + [False]+[False])
                                    },
                                {"title": "Homozygote",
                                    "annotations": []}]),
                        dict(label="Het",
                            method="update",
                            args=[{"visible": ((false3 * 3)+ true3 + false3 + ["legendonly"] + [False, False, False , "legendonly", False] + per_layout 
                                    + [True] + [False] +[False]
                                    + [False] + [True]+[False])
                                    },
                                {"title": "Het",
                                    "annotations": []}]),
                        dict(label="Hemi(XY)",
                            method="update",
                            args=[{"visible": ((false3 * 4)+ true3  + ["legendonly"] + [False, False, False , False, "legendonly"] + per_layout 
                                    + [True] + [False]+[False]
                                    + [False] + [False]+[True])
                                    },
                                {"title": "Hemi",
                                    "annotations": []}]),
                    ]),
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0,
                xanchor="left",
                y=button_layer_1_height,
                yanchor="top"
                )
            ])

    for menu in fig.layout.updatemenus:
        for btn in menu.buttons:
            btn.args[0]['visible'] = list(btn.args[0]['visible']) + ['legendonly']

    fig.update_layout(
        yaxis=dict(range=[-1.1, 1.1], fixedrange=True, domain=[0.30, 1.0]),
        xaxis2=dict(
            matches='x',
            overlaying='x',
            visible=False,
        ),
        yaxis2=dict(overlaying="y", visible=False, fixedrange=True, range=[-1.1, 1.1]),
        yaxis3=dict(
            domain=[0.23, 0.41],
            fixedrange=True,
            showticklabels=False,
            ticks="",
            showgrid=False,
            zeroline=False,
            side='right',
        ),
        height=600,
        margin=dict(b=50),
    )
    fig.update_layout(dragmode="pan")
    fig.update_xaxes(
        range=[0, length],
        minallowed=0,
        maxallowed=length,
        autorange=False
    )

    fig.update_layout(
        xaxis=dict(
            autorange=False,
            rangeslider=dict(
                visible=True,
                range=[0, len(binding_scores)],
                thickness=0.08,
            )
        )
    )
    return fig

def get_last_coding_exon_number(transcript_id,region_map = None, cds_map = None):
    if region_map is None:
        region_map = get_exon_intron_coor(transcript_id)
    if cds_map is None:
        cds_regions = get_cds_map(transcript_id)

    exons = [r for r in region_map if r["type"] == "exon"]

    last_coding_exon = None
    for i, exon in enumerate(exons, start=1):
        for cds in cds_map:
            if cds["start"] <= exon["end"] and cds["end"] >= exon["start"]:
                last_coding_exon = i
                break
    return last_coding_exon

def highlight_by_database(val):
    if val == 'ClinVar':
        return "background-color: rgb(254, 169, 164)"
    if val == 'GnomAD':
        return "background-color: rgb(169, 153, 253)"

def get_uniprot_id(transcript_id):
    return TRANSCRIPT_TO_UNIPROT[transcript_id]

def create_truncation_df(trunc_gnomad, trunc_clinvar , gnomad_raw_data, transcript_id, region_map = None, cds_map = None):
    
    if region_map == None:
        region_map = get_exon_intron_coor(transcript_id)
    if cds_map == None:
        cds_map = get_cds_map(transcript_id)

    clinvar_df = []
    for i in trunc_clinvar:
        id = list(i.keys())[0]
        exon_num = find_exon_number(id.split('-')[1], region_map)
        tmp = {
            'id':id,
            'location' : int(id.split('-')[1]),
            'codons': f"{id.split('-')[2]}/{id.split('-')[3]}", 
            'exon':exon_num,
            'protein start':int(re.search(r'\d+', i[id].short_description).group()),
            'protein_consequence': i[id].short_description,
            'database':'ClinVar',
            'consequence_terms':type(i[id]).__name__,
            }
        clinvar_df.append(tmp)
    clinvar_df = pd.DataFrame(clinvar_df)

   
#=======================================================================================================
#   Getting the population data to add to my gnomad dictionary
#=======================================================================================================
    gnomad_df = []
    gnomad_raw_df = pd.DataFrame(gnomad_raw_data)
    is_sex_chrom = gnomad_raw_df['variant_id'].iloc[0][0] in ['X','Y']
    for i in trunc_gnomad:
        id = list(i.keys())[0]
        pop_dict = {}

#=======================================================================================================
#       Getting the population data to add to my gnomad dictionary
        g_sex = [0,0,0,0]
        ex_sex = [0,0,0,0]
        sum_sex = []
        exome = gnomad_raw_df[gnomad_raw_df['variant_id'] == id]['exome'].iloc[0]
        genome = gnomad_raw_df[gnomad_raw_df['variant_id'] == id]['genome'].iloc[0]
        if genome is not None:
            if is_sex_chrom:
                #allele_count, homo, hemi, hetero
                g_sex = [genome['ac'],genome['homozygote_count'],genome['hemizygote_count'],(genome['ac']-((2*genome['homozygote_count'])-genome['hemizygote_count']))]
            else:
                #total, homo, het
                g_sex = [genome['ac'],genome['homozygote_count'],(genome['ac']-genome['homozygote_count'])]
        if exome is not None:
            if is_sex_chrom:
                ex_sex = [exome['ac'],exome['homozygote_count'],exome['hemizygote_count'],(exome['ac']-((2*exome['homozygote_count'])+exome['hemizygote_count']))]
            else:
                ex_sex = [exome['ac'],exome['homozygote_count'],(exome['ac']-(2*exome['homozygote_count']))]
        sum_sex = [a + b for a, b in zip(ex_sex,g_sex)]
#=======================================================================================================        
        exon_num = find_exon_number(id.split('-')[1], region_map)

        tmp = {
            'id':id,
            'location' : int(id.split('-')[1]),
            'codons': f"{id.split('-')[2]}/{id.split('-')[3]}", 
            'exon':exon_num,
            'protein start':int(re.search(r'\d+', i[id].short_description).group()),
            'protein_consequence': i[id].short_description,
            'database':'GnomAD',
            'consequence_terms':type(i[id]).__name__,
            'total_count' : sum_sex[0],
            'homozygote_count' : sum_sex[1],
            'heterozygote count' : sum_sex[2],
            }
        if is_sex_chrom:
            tmp['hemizygote count'] = sum_sex[3]
        gnomad_df.append(tmp)
    gnomad_df = pd.DataFrame(gnomad_df)
    #print(gnomad_df)
    df = pd.concat([clinvar_df, gnomad_df], ignore_index = True)
    return df



def create_missense_df(clinvar_VEP, clinvar_missense_consq, gnomad_VEP, gnomad_missense_consq, gnomad_raw_data, transcript_id):

    df = []
    for i in clinvar_missense_consq:
        id = list(i.keys())[0]
        vep = None
        if clinvar_VEP is not None:
            vep_matches = [x for x in clinvar_VEP if x['id'].split(':')[1].split('-')[0] == id.split('-')[1] and x['id'].split('/')[1] == id.split('-')[3]]
            if vep_matches:
                vep_set = vep_matches[0]
                transcript_matches = [x for x in vep_set['transcript_consequences'] if x['transcript_id'] == transcript_id]
                if transcript_matches:
                    vep = transcript_matches[0]

        match = re.search(r'p\.([A-Z])(\d+)([A-Z])', i[id].short_description)
        df.append({
            'id' : id,
            'protein_consequence' : i[id].short_description,
            'position': int(re.search(r'\d+', i[id].short_description).group()),
            'variation': match.group(1)+"/"+match.group(3),
            'reference' : match.group(1),
            'alternate': match.group(3),
            'cDNA Position': int(id.split('-')[1]),
            'codons': f"{id.split('-')[2]}/{id.split('-')[3]}", 
            'database' : "ClinVar",
            'consequence_terms': type(i[id]).__name__,
            'polyphen prediction': vep.get('polyphen_prediction',None) if vep else None,
            'polyphen score': vep.get('polyphen_score',None) if vep else None,
            'sift prediction': vep.get('sift_prediction',None) if vep else None,
            'sift score': vep.get('sift_score',None) if vep else None,
            'impact':vep.get('impact',None) if vep else None
        })
    gnomad_raw_df = pd.DataFrame(gnomad_raw_data)
    is_sex_chrom = gnomad_raw_df['variant_id'].iloc[0][0] in ['X','Y']
    for i in gnomad_missense_consq:

        id = (list(i.keys())[0])
        g_sex = [0,0,0,0]
        ex_sex = [0,0,0,0]
        sum_sex = []
        #print(f"MISSENSE CHART DEBUGGING: {id}")
        exome = gnomad_raw_df[gnomad_raw_df['variant_id'] == id]['exome'].iloc[0]
        genome = gnomad_raw_df[gnomad_raw_df['variant_id'] == id]['genome'].iloc[0]
        if genome is not None:
            if is_sex_chrom:
                g_sex = [genome['ac'],genome['homozygote_count'],genome['hemizygote_count'],(genome['ac']-((2*genome['homozygote_count'])-genome['hemizygote_count']))]
            else:
                #total, homo, het
                g_sex = [genome['ac'],genome['homozygote_count'],(genome['ac']-genome['homozygote_count'])]
        if exome is not None:
            if is_sex_chrom:
                ex_sex = [exome['ac'],exome['homozygote_count'],exome['hemizygote_count'],(exome['ac']-((2*exome['homozygote_count'])+exome['hemizygote_count']))]
            else:
                ex_sex = [exome['ac'],exome['homozygote_count'],(exome['ac']-(2*exome['homozygote_count']))]
        sum_sex = [a + b for a, b in zip(ex_sex,g_sex)]
        if gnomad_VEP is not None:
            vep_matches = [x for x in gnomad_VEP if x['id'].split(':')[1].split('-')[0] == id.split('-')[1] and x['id'].split('/')[1] == id.split('-')[3]]
            if vep_matches:
                vep_set = vep_matches[0]
                transcript_matches = [x for x in vep_set['transcript_consequences'] if x['transcript_id'] == transcript_id]
                if transcript_matches:
                    vep = transcript_matches[0]

        match = re.search(r'p\.([A-Z])(\d+)([A-Z])', i[id].short_description)
        to_append = {
            'id' : id,
            'protein_consequence' : i[id].short_description,
            'position': int(re.search(r'\d+', i[id].short_description).group()),
            'variation': match.group(1)+"/"+match.group(3),
            'reference' : match.group(1),
            'alternate': match.group(3),
            'cDNA Position': int(id.split('-')[1]),
            'codons': f"{id.split('-')[2]}/{id.split('-')[3]}", 
            'database' : "GnomAD",
            'consequence_terms': type(i[id]).__name__,
            'polyphen prediction': vep.get('polyphen_prediction',None) if vep else None,
            'polyphen score': vep.get('polyphen_score',None) if vep else None,
            'sift prediction': vep.get('sift_prediction',None) if vep else None,
            'sift score': vep.get('sift_score',None) if vep else None,
            'impact':vep.get('impact',[]) if vep else None,
            'total_count' : sum_sex[0],
            'homozygote_count' : sum_sex[1],
            'heterozygote count' : sum_sex[2],
            }
        if is_sex_chrom:
            to_append['hemizygote count'] = sum_sex[3]
        df.append(to_append)

    df = pd.DataFrame(df)
    return df


@lru_cache(maxsize=512)
def get_transcripts(gene_name):
    genome = get_genome()
    genes = genome.genes_by_name(gene_name)
    pc_transcripts = [
        (t.transcript_id, len((genome.transcript_by_id(t.transcript_id)).protein_sequence))
        for t in genes[0].transcripts if t.biotype == "protein_coding"
    ]
    return pc_transcripts
