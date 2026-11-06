#!/bin/bash
# Download the annotation files the fine-mapping pipeline and the notebooks read.
#
#   GENCODE v50 basic, GRCh37 mapping   gene tracks in the locuszoom figures,
#                                       and nearest and overlapping genes
#                                       (fine_mapping/annotate_credible_sets.R)
#   ENCODE V3 cCREs, GRCh38             candidate cis-regulatory elements, lifted
#                                       to hg19 by lift_regulatory_annotations.sh
#   Rinker et al. 2020 supplement       reintroduced ancestral allele calls
#
# The Iasi et al. 2025 adaptively introgressed regions are prepared by hand; see
# the README. Run on a dev node, since compute nodes have no internet access.
# Files already present are not downloaded again.
#
# Usage:
#   ./retrieve_annotations.sh

set -euo pipefail

PROJECT_DIR="/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk"
ANNOTATION_DIR="${PROJECT_DIR}/data/annotations"
RAA_DIR="${PROJECT_DIR}/data/reintroduced_ancestral_alleles"
GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"

mkdir -p "$ANNOTATION_DIR" "$RAA_DIR"

# download <url> <directory>
download() {
    local url="$1" directory="$2"
    echo "$(basename "$url")"
    wget --quiet --no-clobber --directory-prefix "$directory" "$url"
}

download "${GENCODE_URL}/release_50/GRCh37_mapping/gencode.v50lift37.basic.annotation.gff3.gz" "$ANNOTATION_DIR"
download "https://downloads.wenglab.org/V3/GRCh38-cCREs.bed" "$ANNOTATION_DIR"
download "https://github.com/DaRinker/Neanderthal.reintroduction/raw/master/Rinker_Suppl.File1.tar.gz" "$RAA_DIR"

echo "Done. Next: ./lift_regulatory_annotations.sh to lift the cCREs to hg19."
