"""
gene_pipeline.py — refactored from the original app's main data-fetching block.
Handles the three-tier cache logic exactly as the original script did:
  1) full local pickle cache hit -> load directly
  2) partial match (different transcript, same gene) -> reuse genomic variant calls, recompute the rest
  3) no cache -> fetch everything fresh from ClinVar / gnomAD / Ensembl
Returns the same gene_dict structure the original Streamlit app built.
"""
import os
import glob
import json
import time
import pickle
import logging

from core import (
    get_sequence_from_transcript,
    calculate_consequences,
    filter_consequences_by_type,
    get_ensemblVEP_aa_consequences,
    get_exon_intron_coor,
    get_cds_map,
    fetch_all_clinvar,
    fetch_all_gnomad,
)

logger = logging.getLogger("variant_viewer")

AIUPRED_SCORE_FILE = "postsynaptic_genes_aiupred_score_list.json"


def get_or_build_gene_dict(gene, transcript_id):
    """
    Fetches (or loads from cache) all ClinVar/gnomAD/AIUPred data needed to
    render the missense and truncation charts + structure viewer for one
    gene/transcript pair. Mirrors the original Streamlit app's caching logic.
    """
    start = time.time()
    matches = glob.glob(f"{gene}-*_sheet.pkl")
    if os.path.exists(f"{gene}-{transcript_id}_sheet.pkl"):
        with open(f"{gene}-{transcript_id}_sheet.pkl", "rb") as file:
            gene_dict = pickle.load(file)
        if gene_dict.get("clinvar_missense_VEP_annotations") is None:
            missense_clinvar = gene_dict["clinvar_missense_variants"]
            clinvar_missense_VEP_annotations = get_ensemblVEP_aa_consequences([list(i.keys())[0] for i in missense_clinvar])
            with open(f"{gene}-{transcript_id}_sheet.pkl", "wb") as file:
                pickle.dump(gene_dict, file)
        # if gene_dict.get("gnomad_missense_VEP_annotations") is None:
        #     missense_gnomad = gene_dict["gnomad_missense_variants"]
        #     gnomad_missense_VEP_annotations = get_ensemblVEP_aa_consequences([list(i.keys())[0] for i in missense_gnomad])
        #     with open(f"{gene}-{transcript_id}_sheet.pkl", "wb") as file:
        #         pickle.dump(gene_dict, file)

    elif matches:
        with open(matches[0],"rb") as file:
            old_gene_dict = pickle.load(file)

        sequence = get_sequence_from_transcript(transcript_id)
        length = len(sequence)
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully fetched protein sequence for {gene} locally... ")

        consq_clinvar = calculate_consequences([i['variant_id'] for i in old_gene_dict["clinvar_genomic_variants"] if 'variant_id' in i], transcript_id)
        consq_gnomad = calculate_consequences([i['variant_id'] for i in old_gene_dict["gnomad_genomic_variants"]], transcript_id)



        missense_clinvar = filter_consequences_by_type(consq_clinvar, ['missense'])
        missense_gnomad = filter_consequences_by_type(consq_gnomad, ['missense'])

        print(f"filtered for {len(missense_clinvar)} clinvar missense variants")
        print(f"filtered for {len(missense_gnomad)} gnomad missense variants")
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully calculated {gene} missense(clinvar: {len(missense_clinvar)} / gnomad: {len(missense_gnomad)}) from clinvar and gnomad variants... ==> (writing locally)")

        trunc_clinvar = filter_consequences_by_type(consq_clinvar, ['FrameShift','FrameShiftTruncation','PrematureStop','StartLoss','StopLoss'])
        trunc_gnomad = filter_consequences_by_type(consq_gnomad, ['FrameShift','FrameShiftTruncation','PrematureStop','StartLoss','StopLoss'])
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully calculated {gene} truncation (clinvar: {len(trunc_clinvar)} / gnomad: {len(trunc_gnomad)})from clinvar and gnomad variants... ==> (writing locally)")


        clinvar_missense_VEP_annotations = get_ensemblVEP_aa_consequences([list(i.keys())[0] for i in missense_clinvar])
        if clinvar_missense_VEP_annotations is not None:
            print(f"fetched for {len(clinvar_missense_VEP_annotations)} clinvar missense VEP annotations")
            logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully fetched {gene} ClinVar missense variants consequence predictions from EnsemblVEP... ==> (writing locally)")
        else:
            print("failed to fetch clinvar missense VEP annotations")
            logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Failed to fetch {gene} ClinVar missense variants consequence predictions from EnsemblVEP... Will reattempt on rerun")   
        region_map = get_exon_intron_coor(transcript_id)
        print(f"fetched intron exon map: {len(region_map)} ")
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully fetched exon/intron map...")
        cds_regions = get_cds_map(transcript_id)
        print(f"fetched coding region map: {len(cds_regions)} ")
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully fetched coding region coordinates...")

        with open("postsynaptic_genes_aiupred_score_list.json", "r") as f:
            score_dict = json.load(f)
        score_dict = {k: v for d in score_dict for k, v in d.items()}
        if score_dict.get(gene,{}).get(transcript_id) is None:
            binding_scores = None
            disorder_scores = None
            logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— No locally stored precalculated AIUPred3 Binding/Disorder Prediction Scores found for gene: {gene}; transcript_id {transcript_id}")

        else:
            binding_scores = score_dict[gene][transcript_id]['binding_scores']
            disorder_scores = score_dict[gene][transcript_id]['disorder_scores']
            logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Found locally stored precalculated AIUPred3 Binding/Disorder Prediction Scores for gene: {gene}; transcript_id {transcript_id}")

        with open("postsynaptic_genes_aiupred_score_list.json", "r") as f:
            score_dict = json.load(f)

        gene_dict = {
            "gene_name" :                               gene,
            "transcript_id" :                           transcript_id,
            "sequence":                                 sequence,
            "length":                                   length,
            "clinvar_genomic_variants" :                old_gene_dict["clinvar_genomic_variants"],
            "gnomad_genomic_variants" :                 old_gene_dict["gnomad_genomic_variants"],
            "clinvar_consequence" :                     consq_clinvar,
            "gnomad_consequence" :                      consq_gnomad,
            "clinvar_missense_variants":                missense_clinvar,
            "gnomad_missense_variants" :                missense_gnomad,
            "clinvar_truncation_variants" :             trunc_clinvar,
            "gnomad_truncation_variants" :              trunc_gnomad,
            "clinvar_missense_VEP_annotations" :        clinvar_missense_VEP_annotations,
            "gnomad_missense_VEP_annotations" :         None,
            "exon_intron_map" :                         region_map,
            "cds_regions" :                             cds_regions,
            "aiupred_binding" :                         binding_scores,
            "aiupred_disorder" :                        disorder_scores
        }
        with open (f"{gene}-{transcript_id}_sheet.pkl","wb") as file:
            pickle.dump(gene_dict,file)
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully stored variant information locally...")

    if not os.path.exists(f"{gene}-{transcript_id}_sheet.pkl"):
        sequence = get_sequence_from_transcript(transcript_id)
        length = len(sequence)
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully fetched protein sequence for {gene} locally... ")

        gene_all_clinvar = fetch_all_clinvar(gene)
        print(f"fetched {len(gene_all_clinvar)} clinvar variants")
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully fetched all {gene} variants from clinvar... ==> (writing locally)")

        gene_all_gnomad = fetch_all_gnomad(gene)
        print(f"fetched {len(gene_all_gnomad)} gnomad variants")
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully fetched all {gene} variants from gnomad... ==> (writing locally)")

        consq_clinvar = calculate_consequences([i['variant_id'] for i in gene_all_clinvar if 'variant_id' in i], transcript_id)
        print(f"calculated {len(consq_clinvar)} clinvar consequences")
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully calculated all {gene} variants from ClinVar... ==> (writing locally)")

        consq_gnomad = calculate_consequences([i['variant_id'] for i in gene_all_gnomad], transcript_id)
        print(f"calculated {len(consq_gnomad)} gnomad consequences")
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully calculated all {gene} variants from GnomAD... ==> (writing locally)")

        trunc_clinvar = filter_consequences_by_type(consq_clinvar, ['FrameShift','FrameShiftTruncation','PrematureStop','StartLoss','StopLoss'])
        trunc_gnomad = filter_consequences_by_type(consq_gnomad, ['FrameShift','FrameShiftTruncation','PrematureStop','StartLoss','StopLoss'])
        print(f"filtered for {len(trunc_clinvar)} clinvar truncating variants")
        print(f"filtered for {len(trunc_gnomad)} gnomad truncating variants")
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully extracted {gene} early truncation variants from clinvar and gnomad... ==> (writing locally)")

        missense_clinvar = filter_consequences_by_type(consq_clinvar, ['missense'])
        missense_gnomad = filter_consequences_by_type(consq_gnomad, ['missense'])
        print(f"filtered for {len(missense_clinvar)} clinvar missense variants")
        print(f"filtered for {len(missense_gnomad)} gnomad missense variants")
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully extracted {gene} missense variants from clinvar and gnomad... ==> (writing locally)")

        clinvar_missense_VEP_annotations = get_ensemblVEP_aa_consequences([list(i.keys())[0] for i in missense_clinvar])
        if clinvar_missense_VEP_annotations is not None:
            print(f"fetched for {len(clinvar_missense_VEP_annotations)} clinvar missense VEP annotations")
            logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully fetched {gene} ClinVar missense variants consequence predictions from EnsemblVEP... ==> (writing locally)")
        else:
            print("failed to fetch clinvar missense VEP annotations")
            logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Failed to fetch {gene} ClinVar missense variants consequence predictions from EnsemblVEP... Will reattempt on rerun")   

        region_map = get_exon_intron_coor(transcript_id)
        print(f"fetched intron exon map: {len(region_map)} ")
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully fetched exon/intron map...")
        cds_regions = get_cds_map(transcript_id)
        print(f"fetched coding region map: {len(cds_regions)} ")
        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully fetched coding region coordinates...")

        with open("postsynaptic_genes_aiupred_score_list.json", "r") as f:
            score_dict = json.load(f)
        score_dict = {k: v for d in score_dict for k, v in d.items()}
        if score_dict.get(gene,{}).get(transcript_id) is None:
            binding_scores = None
            disorder_scores = None
            logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— No locally stored precalculated AIUPred3 Binding/Disorder Prediction Scores found for gene: {gene}; transcript_id {transcript_id}")

        else:
            binding_scores = score_dict[gene][transcript_id]['binding_scores']
            disorder_scores = score_dict[gene][transcript_id]['disorder_scores']
            logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Found locally stored precalculated AIUPred3 Binding/Disorder Prediction Scores for gene: {gene}; transcript_id {transcript_id}")

        with open("postsynaptic_genes_aiupred_score_list.json", "r") as f:
            score_dict = json.load(f)

        gene_dict = {
            "gene_name" :                               gene,
            "transcript_id" :                           transcript_id,
            "sequence":                                 sequence,
            "length":                                   length,
            "clinvar_genomic_variants" :                gene_all_clinvar,
            "gnomad_genomic_variants" :                 gene_all_gnomad,
            "clinvar_consequence" :                     consq_clinvar,
            "gnomad_consequence" :                      consq_gnomad,
            "clinvar_missense_variants":                missense_clinvar,
            "gnomad_missense_variants" :                missense_gnomad,
            "clinvar_truncation_variants" :             trunc_clinvar,
            "gnomad_truncation_variants" :              trunc_gnomad,
            "clinvar_missense_VEP_annotations" :        clinvar_missense_VEP_annotations,
            "gnomad_missense_VEP_annotations" :         None,
            "exon_intron_map" :                         region_map,
            "cds_regions" :                             cds_regions,
            "aiupred_binding" :                         binding_scores,
            "aiupred_disorder" :                        disorder_scores
        }
        with open (f"{gene}-{transcript_id}_sheet.pkl","wb") as file:
            pickle.dump(gene_dict,file)

        logger.info(f"{str(round(time.time()-start, 1)) + 's':<15}— Successfully stored variant information locally...")

    return gene_dict
