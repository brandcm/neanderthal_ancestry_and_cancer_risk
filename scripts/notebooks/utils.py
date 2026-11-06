"""Shared helpers for analysis.ipynb."""

import pandas as pd
import pybedtools as pbt

# Columns of the ENCODE SCREEN cCRE BED files, in file order.
SCREEN_CCRE_COLUMNS = [
    'cre_chr', 'cre_start', 'cre_end',
    'dhs_accession', 'ccre_accession', 'ccre_class',
]


def prefer_current(path):
    """Return the `_current` version of a pipeline output if it exists.

    The fine-mapping pipeline writes `<name>_current.<ext>` beside the stable
    file after a partial rerun, and the R figure helpers prefer it, so the
    notebooks do the same.
    """
    current = path.with_name(f'{path.stem}_current{path.suffix}')
    return current if current.exists() else path


def clean_cancer_names(cancers, ancestry_groups):
    """Strip the ancestry, author and accession from pipeline cancer keys.

    'basal_cell_carcinoma_EUR_Seviiri_GCST90137411' -> 'basal_cell_carcinoma'
    """
    return [
        cancer.split(f'_{ancestry}_')[0]
        for cancer, ancestry in zip(cancers, ancestry_groups)
    ]


def get_credible_set(variants, cancer, locus_id):
    """Rows of `variants` for one cancer and fine-mapped locus."""
    return variants[
        (variants['cancer'] == cancer) & (variants['locus_id'] == locus_id)
    ]


def intersect_with_cres(
    variants, cre_path, cre_columns=SCREEN_CCRE_COLUMNS, chrom_prefix='chr'
):
    """Return variants that overlap a CRE, with the overlapping CRE's columns.

    Expects 'chr' and 'pos' (1-based) columns in `variants`. A variant
    overlapping several CREs appears once per CRE. CRE start positions are
    0-based, as in the BED file.
    """
    variants = variants.reset_index(drop=True)

    variant_bed = pbt.BedTool.from_dataframe(pd.DataFrame({
        'chrom': chrom_prefix + variants['chr'].astype(str),
        'start': variants['pos'] - 1,  # BED is 0-based, half-open
        'end': variants['pos'],
        'row': variants.index,
    }))

    hits = variant_bed.intersect(pbt.BedTool(str(cre_path)), wa=True, wb=True)
    if hits.count() == 0:
        return variants.iloc[0:0]

    hits = hits.to_dataframe(header=None, disable_auto_names=True)
    cres = hits.iloc[:, 3:]
    cres.columns = ['row'] + list(cre_columns)

    return variants.merge(cres, left_index=True, right_on='row').drop(
        columns='row'
    )


def summarize_phenotypes(open_targets):
    """One row per Open Targets trait, with the lead-variant beta direction.

    Betas are for each credible set's lead variant, so the direction reflects
    the Neanderthal allele only when the lead's alt allele lies on the
    Neanderthal haplotype.
    """
    summary = (
        open_targets
        .drop_duplicates('Open Targets Study Locus ID')
        .groupby(['Mapped Trait', 'Reported Trait'])
        .agg(
            n_studies=('Study ID', 'nunique'),
            beta=('Lead Variant Beta', 'median'),
            min_p=('Lead Variant P-Value', 'min'),
        )
        .reset_index()
    )
    summary['direction'] = summary['beta'].gt(0).map(
        {True: 'increase', False: 'decrease'}
    )
    return summary.sort_values('min_p')
