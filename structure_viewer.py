"""
structure_viewer.py — refactored from the original app's "Adding the structure figure"
section. Builds the same self-contained HTML/JS 3Dmol.js viewer document (variant density
overlays, binding/disorder/hydrophobicity/charge coloring, surface mesh controls, etc.)
that the original app rendered via st.components.v1.html(). Logic is unchanged; this
version returns the HTML string instead of handing it to Streamlit, so it can be served
directly by a Flask route and embedded in an <iframe>.
"""
import uuid
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from core import assemble_indexed_residue_list, convolve_v3, get_alphafold_cif, get_uniprot_id


def build_structure_html(gene_dict, df_missense, df_trunc, transcript_id):
    """
    Returns a full standalone HTML document (string) for the 3D structure viewer,
    with all variant-density and physicochemical coloring options wired up client-side.
    """
    is_sex_chrom = df_trunc['id'].iloc[0][0] in ['X','Y']
    clinvar_indexed = [i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_missense[df_missense['database'] == "ClinVar"]['position'])],gene_dict['length'])]

    clinvar_homo_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_missense[df_missense['database'] == "ClinVar"]['position'])],gene_dict['length'])],
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['homozygote_count'] > 0)]['position'])],gene_dict['length']))
    )]

    clinvar_homo_het_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_missense[df_missense['database'] == "ClinVar"]['position'])],gene_dict['length'])],
    [1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['homozygote_count'] > 0)]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['heterozygote count'] > 0) ]['position'])],gene_dict['length']))
    ]))]

    clinvar_homo_het_hemi_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_missense[df_missense['database'] == "ClinVar"]['position'])],gene_dict['length'])],
    [1 if a or b or c else 0 for a, b, c in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['homozygote_count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['heterozygote count'] > 0)]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['hemizygote count'] > 0) ]['position'])],gene_dict['length']))
    ]))] if is_sex_chrom else []

    clinvar_homo_hemi_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_missense[df_missense['database'] == "ClinVar"]['position'])],gene_dict['length'])],
    [1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['homozygote_count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['hemizygote count'] > 0) ]['position'])],gene_dict['length']))
    ]))] if is_sex_chrom else []

    clinvar_het_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_missense[df_missense['database'] == "ClinVar"]['position'])],gene_dict['length'])],
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['heterozygote count'] > 0) ]['position'])],gene_dict['length']))
    )]

    clinvar_het_hemi_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_missense[df_missense['database'] == "ClinVar"]['position'])],gene_dict['length'])],
    [1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['heterozygote count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['hemizygote count'] > 0) ]['position'])],gene_dict['length']))
    ]))] if is_sex_chrom else []

    clinvar_hemi_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_missense[df_missense['database'] == "ClinVar"]['position'])],gene_dict['length'])],
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['hemizygote count'] > 0) ]['position'])],gene_dict['length'])
    ))] if is_sex_chrom else []

    homo_indexed = assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['homozygote_count'] > 0) ]['position'])],gene_dict['length'])

    homo_het_indexed = [1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['homozygote_count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['heterozygote count'] > 0) ]['position'])],gene_dict['length']))
    ] if is_sex_chrom else []
    homo_het_hemi_indexed = [1 if a or b or c else 0 for a, b, c in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['homozygote_count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['heterozygote count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['hemizygote count'] > 0) ]['position'])],gene_dict['length']))
    ] if is_sex_chrom else []
    homo_hemi_indexed = [1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['homozygote_count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['hemizygote count'] > 0) ]['position'])],gene_dict['length']))
    ] if is_sex_chrom else []
    het_indexed = assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['heterozygote count'] > 0) ]['position'])],gene_dict['length'])

    het_hemi_indexed = [1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['heterozygote count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['hemizygote count'] > 0) ]['position'])],gene_dict['length']))
    ] if is_sex_chrom else []
    hemi_indexed = assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['hemizygote count'] > 0) ]['position'])],gene_dict['length']) if is_sex_chrom else []

    ROLLING_WINDOW_SLIDER_VALUE  = 10

    def normalize(input_list, range_min =0 , range_max = 1):
        mn, mx = min(input_list), max(input_list)
        return [(i-mn)/(mx-mn) for i in input_list]

    clinvar_density = normalize(convolve_v3(assemble_indexed_residue_list([str(i) for i in list(df_missense[df_missense['database'] == "ClinVar"]['position'])],gene_dict['length']),gene_dict['length']),ROLLING_WINDOW_SLIDER_VALUE)
    homo_density = normalize(convolve_v3(assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['homozygote_count'] > 0)]['position'])],gene_dict['length']),gene_dict['length']),ROLLING_WINDOW_SLIDER_VALUE)
    homo_het_density = normalize(convolve_v3([1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['homozygote_count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['heterozygote count'] > 0) ]['position'])],gene_dict['length']))
    ],gene_dict['length']),ROLLING_WINDOW_SLIDER_VALUE)
    homo_het_hemi_density = normalize(convolve_v3([1 if a or b or c else 0 for a, b, c in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['homozygote_count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['heterozygote count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['hemizygote count'] > 0) ]['position'])],gene_dict['length']))
    ],gene_dict['length']),ROLLING_WINDOW_SLIDER_VALUE) if is_sex_chrom else []
    homo_hemi_density = normalize(convolve_v3([1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['homozygote_count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['hemizygote count'] > 0) ]['position'])],gene_dict['length']))
    ],gene_dict['length']),ROLLING_WINDOW_SLIDER_VALUE) if is_sex_chrom else []
    het_hemi_density = normalize(convolve_v3([1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['heterozygote count'] > 0) ]['position'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['hemizygote count'] > 0) ]['position'])],gene_dict['length']))
    ],gene_dict['length']),ROLLING_WINDOW_SLIDER_VALUE) if is_sex_chrom else []

    het_density = normalize(convolve_v3(assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['heterozygote count'] > 0)]['position'])],gene_dict['length']),gene_dict['length']),ROLLING_WINDOW_SLIDER_VALUE)
    hemi_density = normalize(convolve_v3(assemble_indexed_residue_list([str(i) for i in list(df_missense[(df_missense['database'] == "GnomAD") & (df_missense['hemizygote count'] > 0)]['position'])],gene_dict['length']),gene_dict['length']),ROLLING_WINDOW_SLIDER_VALUE) if is_sex_chrom else []

    clinvar_homo_density = np.subtract(homo_density,clinvar_density)
    clinvar_het_density = np.subtract(het_density,clinvar_density)
    clinvar_hemi_density = np.subtract(hemi_density,clinvar_density) if is_sex_chrom else []
    clinvar_homo_het_density = np.subtract(homo_het_density,clinvar_density)
    clinvar_het_hemi_density = np.subtract(het_hemi_density,clinvar_density) if is_sex_chrom else []
    clinvar_homo_hemi_density = np.subtract(homo_hemi_density,clinvar_density) if is_sex_chrom else []
    clinvar_homo_het_hemi_density = np.subtract(homo_het_hemi_density,clinvar_density) if is_sex_chrom else []

    clinvar_trunc_indexed = [j/2 for j in assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "ClinVar")]['protein start'])],gene_dict['length'])]
    homo_trunc_indexed = assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['homozygote_count'] > 0) ]['protein start'])],gene_dict['length'])
    het_trunc_indexed = assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['heterozygote count'] > 0) ]['protein start'])],gene_dict['length'])
    hemi_trunc_indexed = assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['hemizygote count'] > 0) ]['protein start'])],gene_dict['length']) if is_sex_chrom else []

    clinvar_homo_trunc_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_trunc[df_trunc['database'] == "ClinVar"]['protein start'])],gene_dict['length'])],
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['homozygote_count'] > 0)]['protein start'])],gene_dict['length']))
    )]
    clinvar_homo_het_trunc_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_trunc[df_trunc['database'] == "ClinVar"]['protein start'])],gene_dict['length'])],
    [1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['homozygote_count'] > 0) ]['protein start'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['heterozygote count'] > 0) ]['protein start'])],gene_dict['length']))
    ]))]
    clinvar_homo_het_hemi_trunc_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_trunc[df_trunc['database'] == "ClinVar"]['protein start'])],gene_dict['length'])],
    [1 if a or b or c else 0 for a, b, c in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['homozygote_count'] > 0) ]['protein start'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['heterozygote count'] > 0) ]['protein start'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['hemizygote count'] > 0) ]['protein start'])],gene_dict['length']))
    ]))] if is_sex_chrom else []
    clinvar_homo_hemi_trunc_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_trunc[df_trunc['database'] == "ClinVar"]['protein start'])],gene_dict['length'])],
    [1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['homozygote_count'] > 0) ]['protein start'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['hemizygote count'] > 0) ]['protein start'])],gene_dict['length']))
    ]))] if is_sex_chrom else []
    clinvar_het_hemi_trunc_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_trunc[df_trunc['database'] == "ClinVar"]['protein start'])],gene_dict['length'])],
    [1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['heterozygote count'] > 0) ]['protein start'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['hemizygote count'] > 0) ]['protein start'])],gene_dict['length']))
    ]))] if is_sex_chrom else []
    clinvar_het_trunc_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_trunc[df_trunc['database'] == "ClinVar"]['protein start'])],gene_dict['length'])],
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['heterozygote count'] > 0)]['protein start'])],gene_dict['length']))
    )]
    clinvar_hemi_trunc_indexed = [
        sum(i) for i in list(zip([i/2 for i in assemble_indexed_residue_list([str(i) for i in list(df_trunc[df_trunc['database'] == "ClinVar"]['protein start'])],gene_dict['length'])],
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['hemizygote count'] > 0)]['protein start'])],gene_dict['length']))
    )] if is_sex_chrom else []

    homo_het_trunc_indexed = [1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['homozygote_count'] > 0) ]['protein start'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['heterozygote count'] > 0) ]['protein start'])],gene_dict['length']))
    ]
    homo_hemi_trunc_indexed = [1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['homozygote_count'] > 0) ]['protein start'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['hemizygote count'] > 0) ]['protein start'])],gene_dict['length']))
    ] if is_sex_chrom else []
    het_hemi_trunc_indexed = [1 if a or b else 0 for a, b in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['heterozygote count'] > 0) ]['protein start'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['hemizygote count'] > 0) ]['protein start'])],gene_dict['length']))
    ] if is_sex_chrom else []
    homo_het_hemi_trunc_indexed = [1 if a or b or c else 0 for a, b, c in zip(
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['homozygote_count'] > 0) ]['protein start'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['heterozygote count'] > 0) ]['protein start'])],gene_dict['length']),
        assemble_indexed_residue_list([str(i) for i in list(df_trunc[(df_trunc['database'] == "GnomAD") & (df_trunc['hemizygote count'] > 0) ]['protein start'])],gene_dict['length']))
    ] if is_sex_chrom else []

    #from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.colors as mcolors

    instance_color_map = {0: 'gray', 0.5: 'red',1: 'blue', 1.5: 'green'}  # adjust colors as needed
    def value_to_color(val):
        color_map = {
            0:   '#808080',  # gray
            0.5: '#ff0000',  # red
            1.0: '#0000ff',  # blue
            1.5: '#00ff00',  # green
        }
        return color_map.get(val, '#808080')  # default to gray if not found
    cmap = plt.cm.RdBu_r
    seismic_cmap = plt.cm.seismic
    colors = [mcolors.to_hex(seismic_cmap(i / 99)) for i in range(100)]
    viridis_cmap = plt.cm.viridis
    viridis_colors = [mcolors.to_hex(viridis_cmap(i / 99)) for i in range(100)]

    def value_to_index(val):
        return int(float(val) * 99)

    residue_numbers = list(range(0,gene_dict['length']))
    binding_arr = gene_dict.get('aiupred_binding',{})
    disorder_arr = gene_dict.get('aiupred_disorder',{})

    js_binding = str({r: viridis_colors[value_to_index(v)] for r, v in zip(residue_numbers, binding_arr)}) if binding_arr else str({r: viridis_colors[0] for r in residue_numbers})
    js_disorder = str({r: viridis_colors[value_to_index(v)] for r, v in zip(residue_numbers, disorder_arr)}) if disorder_arr else str({r: viridis_colors[0] for r in residue_numbers})



    js_clinvar_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_indexed)})
    js_clinvar_homo_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_homo_indexed)})
    js_clinvar_homo_het_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_homo_het_indexed)})
    js_clinvar_homo_het_hemi_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_homo_het_hemi_indexed)}) if is_sex_chrom else []
    js_clinvar_homo_hemi_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_homo_hemi_indexed)}) if is_sex_chrom else []
    js_clinvar_het_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_het_indexed)})
    js_clinvar_het_hemi_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_het_hemi_indexed)}) if is_sex_chrom else []
    js_clinvar_hemi_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_hemi_indexed)}) if is_sex_chrom else []

    js_homo_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,homo_indexed)})
    js_homo_het_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,homo_het_indexed)})
    js_homo_het_hemi_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,homo_het_hemi_indexed)}) if is_sex_chrom else []
    js_homo_hemi_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,homo_hemi_indexed)}) if is_sex_chrom else []

    js_het_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,het_indexed)})
    js_het_hemi_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,het_hemi_indexed)}) if is_sex_chrom else []
    js_hemi_instance = str({r: value_to_color(v) for r, v in zip(residue_numbers,hemi_indexed)}) if is_sex_chrom else []

    js_clinvar_convolved = str({r: colors[int((v)*99)] for r, v in zip(residue_numbers,clinvar_density)})
    js_clinvar_homo_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,clinvar_homo_density)})
    js_clinvar_homo_het_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,clinvar_homo_het_density)})
    js_clinvar_homo_het_hemi_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,clinvar_homo_het_hemi_density)}) if is_sex_chrom else []
    js_clinvar_homo_hemi_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,clinvar_homo_hemi_density)}) if is_sex_chrom else []
    js_clinvar_het_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,clinvar_het_density)})
    js_clinvar_het_hemi_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,clinvar_hemi_density)}) if is_sex_chrom else []
    js_clinvar_hemi_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,clinvar_hemi_density)}) if is_sex_chrom else []

    js_homo_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,homo_density)})
    js_homo_het_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,homo_het_density)})
    js_homo_het_hemi_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,homo_het_hemi_density)}) if is_sex_chrom else []
    js_homo_hemi_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,homo_hemi_density)}) if is_sex_chrom else []

    js_het_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,het_density)})
    js_het_hemi_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,het_hemi_density)}) if is_sex_chrom else []
    js_hemi_convolved = str({r: colors[int(v*99)] for r, v in zip(residue_numbers,hemi_density)}) if is_sex_chrom else []

    js_clinvar_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_trunc_indexed)})
    js_clinvar_homo_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_homo_trunc_indexed)})
    js_clinvar_homo_het_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_homo_het_trunc_indexed)})
    js_clinvar_homo_het_hemi_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_homo_het_hemi_trunc_indexed)}) if is_sex_chrom else []
    js_clinvar_homo_hemi_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_homo_hemi_trunc_indexed)}) if is_sex_chrom else []
    js_clinvar_het_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_het_trunc_indexed)})
    js_clinvar_het_hemi_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_het_hemi_trunc_indexed)}) if is_sex_chrom else []
    js_clinvar_hemi_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,clinvar_hemi_trunc_indexed)}) if is_sex_chrom else []

    js_homo_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,homo_trunc_indexed)})
    js_homo_het_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,homo_het_trunc_indexed)})
    js_homo_het_hemi_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,homo_het_hemi_trunc_indexed)}) if is_sex_chrom else []
    js_homo_hemi_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,homo_hemi_trunc_indexed)}) if is_sex_chrom else []

    js_het_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,het_trunc_indexed)})
    js_het_hemi_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,het_hemi_trunc_indexed)}) if is_sex_chrom else []
    js_hemi_trunc = str({r: value_to_color(v) for r, v in zip(residue_numbers,hemi_trunc_indexed)}) if is_sex_chrom else []

    viewer_bg = "#0f1116"  # dark theme, matching the site's look
    cif_string_af = get_alphafold_cif(get_uniprot_id(transcript_id))
    cif_escaped = cif_string_af.replace('\\', '\\\\').replace('`', '\\`')

    #Dropdown works with structure
    import tempfile, webbrowser, uuid

    viewer_id = f"viewer_{uuid.uuid4().hex[:8]}"
    sex_fill = ['<label><input type="checkbox" value="d" onchange="onMissenseCheck()"> GnomAD Hemi</label>',
                '<label><input type="checkbox" value="h"  onchange="onTruncationCheck()"> GnomAD Hemi</label>',
                '<label><input type="checkbox" value="d1" onchange="onGroup3Check()"> GnomAD Hemi</label>',
                'checked.includes("d1");',
                'checked.includes("d");',
                'checked.includes("h");'
                ] if is_sex_chrom else ["","","","false","false","false"]

    html_content = f"""<!DOCTYPE html>
    <html>
    <head>
    <style>
        .dropdown-container {{ position: relative; display: inline-block; font-family: sans-serif; }}
        .dropdown-panel {{ display: none; position: absolute; top: 110%; left: 0; background: white;
            border: 1px solid #ccc; border-radius: 8px; min-width: 220px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.15); z-index: 999; padding: 8px 0; }}
        .toggle-track {{ position: relative; display: inline-block; width: 40px; height: 20px; }}
        .toggle-track input {{ display: none; }}
        .toggle-knob {{ position: absolute; cursor: pointer; top:0; left:0; right:0; bottom:0;
            background: #ccc; border-radius: 20px; transition: .3s; }}
        .toggle-knob:before {{ content: ""; position: absolute; width:16px; height:16px;
            left:2px; bottom:2px; background: red; border-radius: 50%; transition: .3s; }}
        input:checked + .toggle-knob {{ background: blue; }}
        input:checked + .toggle-knob:before {{ transform: translateX(20px); }}
        .category-header {{ padding: 8px 14px; font-weight: bold; cursor: pointer;
            font-size: 13px; color: #333; user-select: none; }}
        .category-header:hover {{ background: #f5f5f5; }}
        .category-items {{ display: none; padding: 0 14px 8px 24px; }}
        .category-items label {{ display: block; font-size: 12px; margin: 4px 0; cursor: pointer; }}
    </style>
    </head>
    <body>
    <div style="width:100%; display:flex; flex-direction:column; align-items:center; gap:10px; padding:16px;">

        <!-- Row 1: PDB input -->
        <div style="display:flex; gap:8px; align-items:center;">
            <input id="pdb-input" type="text" placeholder="Enter PDB ID (e.g. 4HHB)"
                style="padding:6px; font-size:13px; width:220px;">
            <button onclick="loadPDB()" style="padding:6px 12px; font-size:13px;">Load</button>
            <span id="pdb-status" style="font-size:12px; color:grey;"></span>
        </div>

        <!-- Row 2: Dropdowns -->
        <div style="display:flex; gap:12px; align-items:flex-start;">

            <div class="dropdown-container">
                <button onclick="togglePanel('panel', event)">Cartoon ▼</button>
                <div class="dropdown-panel" id="panel">



                    <div style = "border-bottom:1px solid #eee;">
                        <div class="category-header" onclick="toggleCategory('group1', this)">
                            <span class="arrow">▶</span> Variant Coloring
                        </div>
                        <div class="category-items" id="group1">

                        <div style="font-size:11px; font-weight:bold; color:#888; padding:0px 0 0px 0; text-transform:uppercase; letter-spacing:0.5px;">
                                Missense
                            </div>

                        <div style="padding:-20px -20px; display:flex; align-items:center; gap:8px;">
                            <span style="font-size:12px;">Instance</span>
                            <label class="toggle-track">
                                <input type="checkbox" id="cartoon-toggle" onchange="onToggle()">
                                <span class="toggle-knob"></span>
                            </label>
                            <span style="font-size:12px;">Convolved</span>
                        </div>

                            <label><input type="checkbox" value="a" onchange="onMissenseCheck()"> ClinVar</label>
                            <label><input type="checkbox" value="b" onchange="onMissenseCheck()"> GnomAD Homo</label>
                            <label><input type="checkbox" value="c" onchange="onMissenseCheck()"> GnomAD Het</label>
                            {sex_fill[0]}
                            <div style="font-size:11px; font-weight:bold; color:#888; padding:4px 0 2px 0; text-transform:uppercase; letter-spacing:0.5px;">
                                Truncation
                            </div>
                            <label><input type="checkbox" value="e" onchange="onTruncationCheck()"> ClinVar</label>
                            <label><input type="checkbox" value="f" onchange="onTruncationCheck()"> GnomAD Homo</label>
                            <label><input type="checkbox" value="g" onchange="onTruncationCheck()"> GnomAD Het</label>
                            {sex_fill[1]}
                        </div>
                    </div>
                    <div>
                        <div class="category-header" onclick="toggleCategory('group2', this)">
                            <span class="arrow">▶</span> By Type
                        </div>
                        <div class="category-items" id="group2">
                            <label><input type="radio" name="group2-radio" value="spectrum" onchange="onGroup2Radio()"> Spectrum</label>
                            <label><input type="radio" name="group2-radio" value="binding"  onchange="onGroup2Radio()"> Binding</label>
                            <label><input type="radio" name="group2-radio" value="disorder" onchange="onGroup2Radio()"> Disorder</label>
                            <label><input type="radio" name="group2-radio" value="chain"    onchange="onGroup2Radio()"> Chain</label>
                            <label><input type="radio" name="group2-radio" value="residue"  onchange="onGroup2Radio()"> Residue type</label>
                            <label><input type="radio" name="group2-radio" value="ss"       onchange="onGroup2Radio()"> 2° Structure</label>
                            <label><input type="radio" name="group2-radio" value="hydropho" onchange="onGroup2Radio()"> Hydrophobicity</label>
                            <label><input type="radio" name="group2-radio" value="charges"  onchange="onGroup2Radio()"> Partial Charges</label>

                        </div>
                    </div>
                </div>
            </div>

            <div class="dropdown-container">
                <button onclick="togglePanel('panel-2', event)">Surface Mesh ▼</button>
                <div class="dropdown-panel" id="panel-2">
                <div style="padding:8px 14px; border-bottom:1px solid #eee;">
                    <label style="font-size:12px; display:block; margin-bottom:4px;">
                        Opacity: <span id="opacity-val">30</span>%
                    </label>
                <input type="range" id="opacity-slider" min="0" max="100" value="30"
                    oninput="document.getElementById('opacity-val').textContent = this.value;"
                    onchange="updateSurfaceMesh();"
                    style="width:100%; cursor:pointer;">
                </div>

                    <div style = "border-bottom:1px solid #eee;">
                        <div class="category-header" onclick="toggleCategory('group3', this)">
                            <span class="arrow">▶</span> Variant Coloring
                        </div>
                        <div class="category-items" id="group3">
                        <div style="padding:10px 14px; border-bottom:1px solid #eee; display:flex; align-items:center; gap:10px;">
                            <span style="font-size:12px;">Instance</span>
                            <label class="toggle-track">
                                <input type="checkbox" id="surface-toggle" onchange="updateSurfaceMesh()">
                                <span class="toggle-knob"></span>
                            </label>
                            <span style="font-size:12px;">Convolved</span>
                        </div>
                        <div style="font-size:11px; font-weight:bold; color:#888; padding:0px 0 0px 0; text-transform:uppercase; letter-spacing:0.5px;">
                                Missense
                            </div>
                            <label><input type="checkbox" value="a1" onchange="onGroup3Check()"> ClinVar</label>
                            <label><input type="checkbox" value="b1" onchange="onGroup3Check()"> GnomAD Homo</label>
                            <label><input type="checkbox" value="c1" onchange="onGroup3Check()"> GnomAD Het</label>
                            {sex_fill[2]}
                        </div>
                    </div>
                    <div>
                        <div class="category-header" onclick="toggleCategory('group4', this)">
                            <span class="arrow">▶</span> By Type
                        </div>
                        <div class="category-items" id="group4">
                            <label><input type="radio" name="group4-radio" value="white"        onchange="onGroup4Radio()"> White</label>
                            <label><input type="radio" name="group4-radio" value="binding"      onchange="onGroup4Radio()"> Binding</label>
                            <label><input type="radio" name="group4-radio" value="disorder"     onchange="onGroup4Radio()"> Disorder</label>
                            <label><input type="radio" name="group4-radio" value="chain"        onchange="onGroup4Radio()"> Chain</label>
                            <label><input type="radio" name="group4-radio" value="residue"      onchange="onGroup4Radio()"> Residue type</label>
                            <label><input type="radio" name="group4-radio" value="ss"           onchange="onGroup4Radio()"> 2° Structure </label>
                            <label><input type="radio" name="group4-radio" value="hydropho"     onchange="onGroup4Radio()"> Hydrophobicity </label>
                            <label><input type="radio" name="group4-radio" value="charges"      onchange="onGroup4Radio()"> Partial Charges </label>
                        </div>
                    </div>
                </div>
            </div>

        </div>

        <!-- Row 3: 3Dmol viewer -->
        <div id="{viewer_id}" style="width:80%; height:800px; position:relative;"></div>

    </div>

    <script>
        // ── 3Dmol viewer ────────────────────────────────────────────
        var viewer = null;

        var clinvar_InstanceColors ={js_clinvar_instance};
        var clinvar_homo_InstanceColors = {js_clinvar_homo_instance};
        var clinvar_homo_het_InstanceColors = {js_clinvar_homo_het_instance};
        var clinvar_homo_het_hemi_InstanceColors = {js_clinvar_homo_het_hemi_instance};
        var clinvar_homo_hemi_InstanceColors = {js_clinvar_homo_hemi_instance};
        var clinvar_het_InstanceColors = {js_clinvar_het_instance};
        var clinvar_het_hemi_InstanceColors = {js_clinvar_het_hemi_instance};
        var clinvar_hemi_InstanceColors = {js_clinvar_hemi_instance};
        var homo_InstanceColors = {js_homo_instance};
        var homo_het_InstanceColors = {js_homo_het_instance};
        var homo_het_hemi_InstanceColors = {js_homo_het_hemi_instance};
        var homo_hemi_InstanceColors = {js_homo_hemi_instance};
        var het_InstanceColors = {js_het_instance};
        var het_hemi_InstanceColors = {js_het_hemi_instance};
        var hemi_InstanceColors = {js_hemi_instance};

        var clinvar_ConvolvedColors = {js_clinvar_convolved};
        var clinvar_homo_ConvolvedColors = {js_clinvar_homo_convolved};
        var clinvar_homo_het_ConvolvedColors = {js_clinvar_homo_het_convolved};
        var clinvar_homo_het_hemi_ConvolvedColors = {js_clinvar_homo_het_hemi_convolved};
        var clinvar_homo_hemi_ConvolvedColors = {js_clinvar_homo_hemi_convolved};
        var clinvar_het_ConvolvedColors = {js_clinvar_het_convolved};
        var clinvar_het_hemi_ConvolvedColors = {js_clinvar_het_hemi_convolved};
        var clinvar_hemi_ConvolvedColors = {js_clinvar_hemi_convolved};
        var homo_ConvolvedColors = {js_homo_convolved};
        var homo_het_ConvolvedColors = {js_homo_het_convolved};
        var homo_het_hemi_ConvolvedColors = {js_homo_het_hemi_convolved};
        var homo_hemi_ConvolvedColors = {js_homo_hemi_convolved};
        var het_ConvolvedColors = {js_het_convolved};
        var het_hemi_ConvolvedColors = {js_het_hemi_convolved};
        var hemi_ConvolvedColors = {js_hemi_convolved};

        var clinvar_TruncColors ={js_clinvar_trunc};
        var clinvar_homo_TruncColors = {js_clinvar_homo_trunc};
        var clinvar_homo_het_TruncColors = {js_clinvar_homo_het_trunc};
        var clinvar_homo_het_hemi_TruncColors = {js_clinvar_homo_het_hemi_trunc};
        var clinvar_homo_hemi_TruncColors = {js_clinvar_homo_hemi_trunc};
        var clinvar_het_TruncColors = {js_clinvar_het_trunc};
        var clinvar_het_hemi_TruncColors = {js_clinvar_het_hemi_trunc};
        var clinvar_hemi_TruncColors = {js_clinvar_hemi_trunc};
        var homo_TruncColors = {js_homo_trunc};
        var homo_het_TruncColors = {js_homo_het_trunc};
        var homo_het_hemi_TruncColors = {js_homo_het_hemi_trunc};
        var homo_hemi_TruncColors = {js_homo_hemi_trunc};
        var het_TruncColors = {js_het_trunc};
        var het_hemi_TruncColors = {js_het_hemi_trunc};
        var hemi_TruncColors = {js_hemi_trunc};

        var bindingColors = {js_binding}
        var disorderColors = {js_disorder}
        var toleranceValues = {{}};

        var script = document.createElement('script');
        script.src = 'https://3Dmol.org/build/3Dmol-min.js';

        script.onload = function() {{
            viewer = $3Dmol.createViewer(document.getElementById("{viewer_id}"), {{backgroundColor: "{viewer_bg}"}});
            viewer.addModel(`{cif_escaped}`, 'cif');
            viewer.setStyle({{}}, {{cartoon: {{color: 'grey'}}}});
            viewer.zoomTo();
            viewer.render();

            if (!pdbId || pdbId == "AF") {{
                url = 'https://alphafold.ebi.ac.uk/files/AF-' + {get_uniprot_id(transcript_id)} + '-F1-model_v4.cif';
            }}
            else {{
                url = 'https://files.rcsb.org/download/' + pdbId + '.cif';
             }}
        }};


        // ── Dropdown logic ──────────────────────────────────────────
        function togglePanel(id, event) {{
            event.stopPropagation();
            var panel = document.getElementById(id);
            var isOpen = panel.style.display === 'block';
            ['panel', 'panel-2'].forEach(function(pid) {{
                document.getElementById(pid).style.display = 'none';
            }});
            if (!isOpen) panel.style.display = 'block';
        }}

        document.addEventListener('click', function(e) {{
            var containers = document.querySelectorAll('.dropdown-container');
            var inside = false;
            containers.forEach(function(c) {{ if (c.contains(e.target)) inside = true; }});
            if (!inside) {{
                ['panel', 'panel-2'].forEach(function(pid) {{
                    document.getElementById(pid).style.display = 'none';
                }});
            }}
        }});

        // ── Surface Mesh logic ──────────────────────────────────────────
        var currentSurface = null;

        function updateSurfaceMesh() {{
            if (!viewer) return;

            var checked = Array.from(document.querySelectorAll('#group3 input:checked'))
                            .map(function(cb) {{ return cb.value; }});
            var opacity = document.getElementById('opacity-slider').value / 100;

            // remove existing surface first
            if (currentSurface !== null) {{
                viewer.removeSurface(currentSurface);
                currentSurface = null;
            }}

            // nothing checked — leave surface off
            var selectedRadio = document.querySelector('input[name="group4-radio"]:checked');
            if (checked.length === 0 && !selectedRadio) {{
                viewer.render();
                return;
            }}

            var isConvolved = document.getElementById('surface-toggle').checked;
            var hasClinvar  = checked.includes('a1');
            var hasHomo     = checked.includes('b1');
            var hasHet      = checked.includes('c1');
            var hasHemi     = {sex_fill[3]}

            var chargeScale = {{
                'ASP': -1.0, 'GLU': -1.0,
                'LYS':  1.0, 'ARG':  1.0,
                'HIS':  0.5,
                'CYS': -0.2, 'TYR': -0.2,
                'SER': -0.1, 'THR': -0.1
            }};
            var atoms = viewer.getModel().selectedAtoms({{}});
            atoms.forEach(function(atom) {{
                atom.charge = chargeScale[atom.resn] !== undefined ? chargeScale[atom.resn] : 0;
            }});

            var hydScale = {{
                'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5,
                'MET': 1.9, 'ALA': 1.8, 'GLY': -0.4, 'THR': -0.7, 'SER': -0.8,
                'TRP': -0.9, 'TYR': -1.3, 'PRO': -1.6, 'HIS': -3.2, 'GLU': -3.5,
                'GLN': -3.5, 'ASP': -3.5, 'ASN': -3.5, 'LYS': -3.9, 'ARG': -4.5
            }};
            var atoms = viewer.getModel().selectedAtoms({{}});
            atoms.forEach(function(atom) {{
                atom.hydrophobicity = hydScale[atom.resn] !== undefined ? hydScale[atom.resn] : 0;
            }});

            var style = {{opacity : opacity}}
            if (isConvolved){{
                if (hasClinvar && !hasHomo && !hasHet && !hasHemi){{  //Clinvar
                    style.colorscheme = {{prop:'resi', map: clinvar_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && hasHomo && !hasHet && !hasHemi){{   //ClinVar Homo
                    style.colorscheme = {{prop:'resi', map: clinvar_homo_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && hasHomo && hasHet && !hasHemi){{     // Clinvar Homo Het
                    style.colorscheme = {{prop:'resi', map: clinvar_homo_het_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && hasHomo && hasHet && hasHemi){{      // Clinvar Homo Het Hemi
                    style.colorscheme = {{prop:'resi', map: clinvar_homo_het_hemi_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && hasHomo && !hasHet && hasHemi){{      // Clinvar Homo Hemi
                    style.colorscheme = {{prop:'resi', map: clinvar_homo_hemi_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && !hasHomo && hasHet && !hasHemi){{     // Clinvar Het
                    style.colorscheme = {{prop:'resi', map: clinvar_het_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && !hasHomo && hasHet && hasHemi){{    // Clinvar Het Hemi
                    style.colorscheme = {{prop:'resi', map: clinvar_het_hemi_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && !hasHomo && !hasHet && hasHemi){{   //  Clinvar Hemi
                    style.colorscheme = {{prop:'resi', map: clinvar_hemi_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && hasHomo && !hasHet && !hasHemi){{   //  Homo
                    style.colorscheme = {{prop:'resi', map: homo_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && hasHomo && hasHet && !hasHemi){{   //  Homo Het
                    style.colorscheme = {{prop:'resi', map: homo_het_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && hasHomo && hasHet && hasHemi){{   //  Homo Het Hemi
                    style.colorscheme = {{prop:'resi', map: homo_het_hemi_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && hasHomo && !hasHet && hasHemi){{   //  Homo Hemi
                    style.colorscheme = {{prop:'resi', map: homo_hemi_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && !hasHomo && hasHet && !hasHemi){{   //  Het
                    style.colorscheme = {{prop:'resi', map: het_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && !hasHomo && hasHet && hasHemi){{   //  Het Hemi
                    style.colorscheme = {{prop:'resi', map: het_hemi_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && !hasHomo && !hasHet && hasHemi){{   //  Hemi
                    style.colorscheme = {{prop:'resi', map: hemi_ConvolvedColors }}
                    style.opacity = opacity
                }}
                else {{ style.color = 'white'
                        style.opacity = opacity}}
            }}
            else{{
                if (hasClinvar && !hasHomo && !hasHet && !hasHemi){{  //Clinvar
                    style.colorscheme = {{prop:'resi', map: clinvar_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && hasHomo && !hasHet && !hasHemi){{   //ClinVar Homo
                    style.colorscheme = {{prop:'resi', map: clinvar_homo_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && hasHomo && hasHet && !hasHemi){{     // Clinvar Homo Het
                    style.colorscheme = {{prop:'resi', map: clinvar_homo_het_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && hasHomo && hasHet && hasHemi){{      // Clinvar Homo Het Hemi
                    style.colorscheme = {{prop:'resi', map: clinvar_homo_het_hemi_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && hasHomo && !hasHet && hasHemi){{      // Clinvar Homo Hemi
                    style.colorscheme = {{prop:'resi', map: clinvar_homo_hemi_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && !hasHomo && hasHet && !hasHemi){{     // Clinvar Het
                    style.colorscheme = {{prop:'resi', map: clinvar_het_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && !hasHomo && hasHet && hasHemi){{    // Clinvar Het Hemi
                    style.colorscheme = {{prop:'resi', map: clinvar_het_hemi_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (hasClinvar && !hasHomo && !hasHet && hasHemi){{   //  Clinvar Hemi
                    style.colorscheme = {{prop:'resi', map: clinvar_hemi_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && hasHomo && !hasHet && !hasHemi){{   //  Homo
                    style.colorscheme = {{prop:'resi', map: homo_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && hasHomo && hasHet && !hasHemi){{   //  Homo Het
                    style.colorscheme = {{prop:'resi', map: homo_het_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && hasHomo && hasHet && hasHemi){{   //  Homo Het Hemi
                    style.colorscheme = {{prop:'resi', map: homo_het_hemi_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && hasHomo && !hasHet && hasHemi){{   //  Homo Hemi
                    style.colorscheme = {{prop:'resi', map: homo_hemi_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && !hasHomo && hasHet && !hasHemi){{   //  Het
                    style.colorscheme = {{prop:'resi', map: het_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && !hasHomo && hasHet && hasHemi){{   //  Het Hemi
                    style.colorscheme = {{prop:'resi', map: het_hemi_InstanceColors }}
                    style.opacity = opacity
                }}
                else if (!hasClinvar && !hasHomo && !hasHet && hasHemi){{   //  Hemi
                    style.colorscheme = {{prop:'resi', map: hemi_InstanceColors }}
                    style.opacity = opacity
                }}
                else {{ style.color = 'white'
                style.opacity = opacity}}
            }}

            var selected = document.querySelector('input[name="group4-radio"]:checked');
            if (selected) {{
                if (selected.value === 'white') {{style.color = '#FFFDD0'}}
                if (selected.value === 'binding') {{ 
                    style.colorscheme = {{prop: 'resi', map: bindingColors}} 
                    style.opacity = opacity }}
                if (selected.value === 'disorder') {{ 
                    style.colorscheme = {{prop: 'resi', map: disorderColors}} 
                    style.opacity = opacity }}
                if (selected.value === 'residue') {{ style.colorscheme = 'amino' }}
                if (selected.value === 'chain') {{ style.colorscheme = 'chainHetatm' }}
                if (selected.value === 'ss') {{ style.colorscheme = 'ssJmol' }}
                if (selected.value == 'hydropho') {{
                    style.colorscheme = {{prop: 'hydrophobicity', gradient: 'rwb', min: -4.5, max: 4.5}}}}
                if (selected.value == 'charges'){{
                    style.colorscheme = {{prop: 'charge', gradient: 'rwb', min: -1, max: 1}}}}
            }}

            viewer.addSurface($3Dmol.SurfaceType.VDW, style).then(function(id) {{
                    currentSurface = id;
                    viewer.render();
            }});
        }}

    function toggleCategory(id, header) {{
        var el = document.getElementById(id);
        var isOpen = el.style.display === 'block';
        el.style.display = isOpen ? 'none' : 'block';
        header.querySelector('.arrow').textContent = isOpen ? '▶' : '▼';}}

    function onToggle() {{
                onGroup1Check();
        }}

    function onMissenseCheck() {{
        const anyMissense = ['a','b','c','d'].some(v => {{
            const el = document.querySelector(`input[value="${{v}}"]`);
            return el && el.checked;
        }});
        if (anyMissense) {{
            ['e','f','g','h'].forEach(v => {{
                const el = document.querySelector(`input[value="${{v}}"]`);
                if (el) el.checked = false;
            }});
            document.getElementById('cartoon-toggle').disabled = false;
        }}
        onGroup1Check();
    }}

    function onTruncationCheck() {{
        const anyTrunc = ['e','f','g','h'].some(v => {{
            const el = document.querySelector(`input[value="${{v}}"]`);
            return el && el.checked;
        }});
        if (anyTrunc) {{
            ['a','b','c','d'].forEach(v => {{
                const el = document.querySelector(`input[value="${{v}}"]`);
                if (el) el.checked = false;
            }});
            document.getElementById('cartoon-toggle').checked = false;
            document.getElementById('cartoon-toggle').disabled = true;
        }}
        onGroup1Check();
    }}


    function onGroup1Check() {{
            document.querySelectorAll('input[name="group2-radio"]')
                    .forEach(function(r) {{ r.checked = false; }});
            var checked = Array.from(document.querySelectorAll('#group1 input:checked'))
                            .map(function(cb) {{ return cb.value; }});
            if (!viewer) return;

            var isConvolved = document.getElementById('cartoon-toggle').checked;
            var hasClinvar  = checked.includes('a');
            var hasHomo     = checked.includes('b');
            var hasHet      = checked.includes('c');
            var hasHemi     = {sex_fill[4]}



            var color = 'spectrum';
            if (isConvolved) {{
                if (hasClinvar && !hasHomo && !hasHet && !hasHemi){{  //Clinvar
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && hasHomo && !hasHet && !hasHemi){{   //ClinVar Homo
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_homo_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && hasHomo && hasHet && !hasHemi){{     // Clinvar Homo Het
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_homo_het_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && hasHomo && hasHet && hasHemi){{      // Clinvar Homo Het Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_homo_het_hemi_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && hasHomo && !hasHet && hasHemi){{      // Clinvar Homo Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_homo_hemi_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && !hasHomo && hasHet && !hasHemi){{     // Clinvar Het
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_het_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && !hasHomo && hasHet && hasHemi){{    // Clinvar Het Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_het_hemi_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && !hasHomo && !hasHet && hasHemi){{   //  Clinvar Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_hemi_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && hasHomo && !hasHet && !hasHemi){{   //  Homo
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: homo_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && hasHomo && hasHet && !hasHemi){{   //  Homo Het
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: homo_het_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && hasHomo && hasHet && hasHemi){{   //  Homo Het Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: homo_het_hemi_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && hasHomo && !hasHet && hasHemi){{   //  Homo Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: homo_hemi_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && !hasHomo && hasHet && !hasHemi){{   //  Het
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: het_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && !hasHomo && hasHet && hasHemi){{   //  Het Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: het_hemi_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && !hasHomo && !hasHet && hasHemi){{   //  Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: hemi_ConvolvedColors}}}}}});
                    viewer.render(); return;}}
            }}
            else {{
                if (hasClinvar && !hasHomo && !hasHet && !hasHemi){{  //Clinvar
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && hasHomo && !hasHet && !hasHemi){{   //ClinVar Homo
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_homo_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && hasHomo && hasHet && !hasHemi){{     // Clinvar Homo Het
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_homo_het_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && hasHomo && hasHet && hasHemi){{      // Clinvar Homo Het Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_homo_het_hemi_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && hasHomo && !hasHet && hasHemi){{      // Clinvar Homo Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_homo_hemi_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && !hasHomo && hasHet && !hasHemi){{     // Clinvar Het
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_het_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && !hasHomo && hasHet && hasHemi){{    // Clinvar Het Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_het_hemi_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (hasClinvar && !hasHomo && !hasHet && hasHemi){{   //  Clinvar Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_hemi_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && hasHomo && !hasHet && !hasHemi){{   //  Homo
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: homo_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && hasHomo && hasHet && !hasHemi){{   //  Homo Het
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: homo_het_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && hasHomo && hasHet && hasHemi){{   //  Homo Het Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: homo_het_hemi_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && hasHomo && !hasHet && hasHemi){{   //  Homo Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: homo_hemi_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && !hasHomo && hasHet && !hasHemi){{   //  Het
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: het_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && !hasHomo && hasHet && hasHemi){{   //  Het Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: het_hemi_InstanceColors}}}}}});
                    viewer.render(); return;}}
                if (!hasClinvar && !hasHomo && !hasHet && hasHemi){{   //  Hemi
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: hemi_InstanceColors}}}}}});
                    viewer.render(); return;}}
            }}

            var hasClinvarTrunc = checked.includes('e');
            var hasHomoTrunc = checked.includes('f');
            var hasHetTrunc = checked.includes('g');
            var hasHemiTrunc = {sex_fill[5]}

            if (hasClinvarTrunc && !hasHomoTrunc && !hasHetTrunc && !hasHemiTrunc){{  //Clinvar
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_TruncColors}}}}}});
                viewer.render(); return;}}
            if (hasClinvarTrunc && hasHomoTrunc && !hasHetTrunc && !hasHemiTrunc){{   //ClinVar Homo
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_homo_TruncColors}}}}}});
                viewer.render(); return;}}
            if (hasClinvarTrunc && hasHomoTrunc && hasHetTrunc && !hasHemiTrunc){{     // Clinvar Homo Het
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_homo_het_TruncColors}}}}}});
                viewer.render(); return;}}
            if (hasClinvarTrunc && hasHomoTrunc && hasHetTrunc && hasHemiTrunc){{      // Clinvar Homo Het Hemi
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_homo_het_hemi_TruncColors}}}}}});
                viewer.render(); return;}}
            if (hasClinvarTrunc && hasHomoTrunc && !hasHetTrunc && hasHemiTrunc){{      // Clinvar Homo Hemi
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_homo_hemi_TruncColors}}}}}});
                viewer.render(); return;}}
            if (hasClinvarTrunc && !hasHomoTrunc && hasHetTrunc && !hasHemiTrunc){{     // Clinvar Het
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_het_TruncColors}}}}}});
                viewer.render(); return;}}
            if (hasClinvarTrunc && !hasHomoTrunc && hasHetTrunc && hasHemiTrunc){{    // Clinvar Het Hemi
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_het_hemi_TruncColors}}}}}});
                viewer.render(); return;}}
            if (hasClinvarTrunc && !hasHomoTrunc && !hasHetTrunc && hasHemiTrunc){{   //  Clinvar Hemi
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: clinvar_hemi_TruncColors}}}}}});
                viewer.render(); return;}}
            if (!hasClinvarTrunc && hasHomoTrunc && !hasHetTrunc && !hasHemiTrunc){{   //  Homo
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: homo_TruncColors}}}}}});
                viewer.render(); return;}}
            if (!hasClinvarTrunc && hasHomoTrunc && hasHetTrunc && !hasHemiTrunc){{   //  Homo Het
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: homo_het_TruncColors}}}}}});
                viewer.render(); return;}}
            if (!hasClinvarTrunc && hasHomoTrunc && hasHetTrunc && hasHemiTrunc){{   //  Homo Het Hemi
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: homo_het_hemi_TruncColors}}}}}});
                viewer.render(); return;}}
            if (!hasClinvarTrunc && hasHomoTrunc && !hasHetTrunc && hasHemiTrunc){{   //  Homo Hemi
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: homo_hemi_TruncColors}}}}}});
                viewer.render(); return;}}
            if (!hasClinvarTrunc && !hasHomoTrunc && hasHetTrunc && !hasHemiTrunc){{   //  Het
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: het_TruncColors}}}}}});
                viewer.render(); return;}}
            if (!hasClinvarTrunc && !hasHomoTrunc && hasHetTrunc && hasHemiTrunc){{   //  Het Hemi
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: het_hemi_TruncColors}}}}}});
                viewer.render(); return;}}
            if (!hasClinvarTrunc && !hasHomoTrunc && !hasHetTrunc && hasHemiTrunc){{   //  Hemi
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{prop:'resi', map: hemi_TruncColors}}}}}});
                viewer.render(); return;}}


            console.log('hasClinvar:', hasClinvar, 'hasHomo:', hasHomo, 'hasHet:', hasHet, 'hasHemi:', hasHemi);
            console.log('clinvar_InstanceColors keys:', Object.keys(clinvar_InstanceColors).length);

            viewer.setStyle({{}}, {{cartoon: {{color: color}}}});
            viewer.render();
            console.log('isConvolved:', isConvolved);
            console.log('hasClinvar:', hasClinvar);
            console.log('clinvar_InstanceColors keys:', Object.keys(clinvar_InstanceColors).length);
        }}


        function onGroup2Radio() {{
            document.querySelectorAll('#group1 input[type="checkbox"]').forEach(function(cb) {{ cb.checked = false; }}); 
                if (!viewer) return;

            var selected = document.querySelector('input[name="group2-radio"]:checked');
                if (!selected) return;

            var chargeScale = {{
                'ASP': -1.0, 'GLU': -1.0,
                'LYS':  1.0, 'ARG':  1.0,
                'HIS':  0.5,
                'CYS': -0.2, 'TYR': -0.2,
                'SER': -0.1, 'THR': -0.1
            }};
            var atoms = viewer.getModel().selectedAtoms({{}});
            atoms.forEach(function(atom) {{
                atom.charge = chargeScale[atom.resn] !== undefined ? chargeScale[atom.resn] : 0;
            }});

            var hydScale = {{
                'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5,
                'MET': 1.9, 'ALA': 1.8, 'GLY': -0.4, 'THR': -0.7, 'SER': -0.8,
                'TRP': -0.9, 'TYR': -1.3, 'PRO': -1.6, 'HIS': -3.2, 'GLU': -3.5,
                'GLN': -3.5, 'ASP': -3.5, 'ASN': -3.5, 'LYS': -3.9, 'ARG': -4.5
            }};
            var atoms = viewer.getModel().selectedAtoms({{}});
            atoms.forEach(function(atom) {{
                atom.hydrophobicity = hydScale[atom.resn] !== undefined ? hydScale[atom.resn] : 0;
            }});
                if (selected.value === 'spectrum') {{
                    viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
                    viewer.render();
                }}
                if (selected.value === 'binding') {{
                    for (var resi in bindingColors) {{
                    viewer.setStyle({{resi: parseInt(resi)}}, 
                    {{cartoon: {{color: bindingColors[resi]}}}});
                    }}viewer.render(); return;}}

                if (selected.value === 'disorder') {{
                    for (var resi in disorderColors) {{
                    viewer.setStyle({{resi: parseInt(resi)}}, 
                    {{cartoon: {{color: disorderColors[resi]}}}});
                    }}viewer.render(); return;}}

                if (selected.value === 'chain') {{
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: 'chainHetatm'}}}});
                    viewer.render();
                    return;
                }}
                if (selected.value === 'residue') {{
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: 'amino'}}}});
                    viewer.render();
                    return;
                }}
                if (selected.value === 'ss') {{
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: 'ssJmol'}}}});
                    viewer.render();
                    return;
                }}
                if (selected.value === 'hydropho') {{
                    for (var resi in window.hydrophobicityCartoonColors) {{
                        viewer.setStyle({{resi: parseInt(resi)}}, {{cartoon: {{color: window.hydrophobicityCartoonColors[resi]}}}});
                    }}
                    viewer.render(); return;
                }}
                if (selected.value === 'charges') {{
                    for (var resi in window.chargeCartoonColors) {{
                        viewer.setStyle({{resi: parseInt(resi)}}, {{cartoon: {{color: window.chargeCartoonColors[resi]}}}});
                    }}
                    viewer.render(); return;
                }}

            }}
        function onGroup3Check() {{ 
            document.querySelectorAll('input[name="group4-radio"]')
                .forEach(function(r) {{ r.checked = false; }}); 
                updateSurfaceMesh();
            }}
        function onGroup4Radio() {{ 
            document.querySelectorAll('#group3 input[type="checkbox"]').forEach(function(cb) {{ cb.checked = false; }}); 
            updateSurfaceMesh(); 
            }}

        document.head.appendChild(script);

        function loadPDB() {{
            if (!viewer) {{ alert('Viewer still loading'); return; }}
            var pdbId = document.getElementById('pdb-input').value.trim().toUpperCase();
            var status = document.getElementById('pdb-status');
            status.textContent = 'Loading...';
            status.style.color = 'grey';
            if (!pdbId || pdbId == "AF") {{
                url = 'https://alphafold.ebi.ac.uk/files/AF-' + {get_uniprot_id(transcript_id)} + '-F1-model_v4.cif';
            }} else {{
                url = 'https://files.rcsb.org/download/' + pdbId + '.cif';
            }}
            fetch(url)
                .then(function(r) {{
                    if (!r.ok) throw new Error('Not found');
                    return r.text();
                }})
                .then(function(data) {{
                    viewer.removeAllModels();
                    viewer.removeAllSurfaces();
                    viewer.addModel(data, 'cif');
                    var hydScale = {{
                        'ILE':4.5,'VAL':4.2,'LEU':3.8,'PHE':2.8,'CYS':2.5,'MET':1.9,'ALA':1.8,
                        'GLY':-0.4,'THR':-0.7,'SER':-0.8,'TRP':-0.9,'TYR':-1.3,'PRO':-1.6,
                        'HIS':-3.2,'GLU':-3.5,'GLN':-3.5,'ASP':-3.5,'ASN':-3.5,'LYS':-3.9,'ARG':-4.5
                    }};
                    var chargeScale = {{
                        'ASP':-1.0,'GLU':-1.0,'LYS':1.0,'ARG':1.0,'HIS':0.5
                    }};

                    function scaleToColor(val, min, max) {{
                        var t = Math.max(0, Math.min(1, (val - min) / (max - min)));
                        var r = Math.round(255 * Math.min(1, 2*t));
                        var b = Math.round(255 * Math.min(1, 2*(1-t)));
                        var g = Math.round(255 * (1 - Math.abs(2*t - 1)));
                        return '#' + [r,g,b].map(function(x){{return x.toString(16).padStart(2,'0')}}).join('');
                    }}

                    window.hydrophobicityCartoonColors = {{}};
                    window.chargeCartoonColors = {{}};
                    viewer.getModel().selectedAtoms({{}}).forEach(function(atom) {{
                        window.hydrophobicityCartoonColors[atom.resi] = scaleToColor(hydScale[atom.resn] || 0, -4.5, 4.5);
                        window.chargeCartoonColors[atom.resi] = scaleToColor(chargeScale[atom.resn] || 0, -1, 1);
                    }});

                    viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
                    viewer.zoomTo();
                    viewer.render();
                    status.textContent = pdbId + ' loaded';
                    status.style.color = 'green';
                    viewer.setHoverable({{}}, true,
                        function(atom, viewer) {{
                            if (!atom.label) {{
                                var val = toleranceValues[atom.resi] !== undefined
                                    ? toleranceValues[atom.resi].toFixed(3) : "N/A";
                                atom.label = viewer.addLabel(
                                    atom.resn + atom.resi + " | tolerance: " + val,
                                    {{position: atom, backgroundColor: "black", fontColor: "white", fontSize: 12}}
                                );
                            }}
                        }},
                        function(atom, viewer) {{
                            if (atom.label) {{
                                viewer.removeLabel(atom.label);
                                delete atom.label;
                            }}
                        }}
                    );
                    viewer.render();

                }})
        }}

        document.getElementById('pdb-input').addEventListener('keydown', function(e) {{
            if (e.key === 'Enter') loadPDB();
        }});

    var idleTimer = null;
    function resetIdleTimer() {{
        clearTimeout(idleTimer);
        viewer.spin(false);  // stop spinning on interaction
        idleTimer = setTimeout(function() {{
            viewer.spin('y', 0.3);  // start spinning after 6.5 seconds idle
        }}, 6500);
    }}

    // Stop spin and reset timer on any user interaction
    document.getElementById('{viewer_id}').addEventListener('mousedown', resetIdleTimer);
    document.getElementById('{viewer_id}').addEventListener('wheel',     resetIdleTimer);
    document.getElementById('{viewer_id}').addEventListener('touchstart', resetIdleTimer);

    resetIdleTimer();  // start the timer on load

    </script>
    </body>
    </html>"""


    return html_content