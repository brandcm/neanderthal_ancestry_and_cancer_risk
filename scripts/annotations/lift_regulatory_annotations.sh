#!/bin/bash
# Lift the ENCODE cCREs from hg38 to hg19 for the chromosomes of the locus
# examples, writing GRCh37-cCREs_<chromosome>.bed.
#
# liftOver only carries the first four BED columns through reliably, so any
# further columns are joined into the name field with "|" beforehand and split
# back out afterwards. Splitting on "|" restores them without the script needing
# to know how many there were.
#
# Usage:
#   ./lift_regulatory_annotations.sh [annotation directory]
#
# Defaults to the project's data/annotations. To cover another locus, add its
# chromosome to CHROMOSOMES.

set -euo pipefail

LIFTOVER_DIR="/wynton/group/capra/bin/liftOver"
CHAIN="${LIFTOVER_DIR}/hg38ToHg19.over.chain.gz"
PROJECT_DIR="/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk"
ANNOTATION_DIR="${1:-${PROJECT_DIR}/data/annotations}"
INPUT="GRCh38-cCREs.bed"

# chr1 PABPC4, chr9 BNC2, chr12 PTHLH and the KRT cluster, chr14 SLC24A4
CHROMOSOMES=(chr1 chr9 chr12 chr14)

for required in "${LIFTOVER_DIR}/liftOver" "${CHAIN}"; do
    [[ -e "$required" ]] || { echo "Missing: $required" >&2; exit 1; }
done

# lift_chromosome <input hg38 bed> <chromosome> <output hg19 bed>
lift_chromosome() {
    local input="$1" chromosome="$2" output="$3"
    local stem="${output%.bed}"

    awk -v OFS='\t' -v chromosome="$chromosome" '$1 == chromosome {
        name = $4
        for (i = 5; i <= NF; i++) name = name "|" $i
        print $1, $2, $3, name
    }' "$input" > "${stem}_hg38.bed"

    local n_input
    n_input=$(wc -l < "${stem}_hg38.bed")
    if [[ "$n_input" -eq 0 ]]; then
        echo "No intervals on ${chromosome} in $(basename "$input") -- check the chromosome naming" >&2
        return 1
    fi

    "${LIFTOVER_DIR}/liftOver" \
        "${stem}_hg38.bed" \
        "${CHAIN}" \
        "${stem}.tmp" \
        "${stem}.unlifted"

    tr '|' '\t' < "${stem}.tmp" > "$output"
    rm -f "${stem}.tmp"

    local n_output n_unlifted
    n_output=$(wc -l < "$output")
    n_unlifted=$(grep -vc '^#' "${stem}.unlifted" || true)
    echo "$(basename "$output"): ${n_output} of ${n_input} intervals lifted (${n_unlifted} unlifted)"
}

[[ -d "$ANNOTATION_DIR" ]] || { echo "No such directory: $ANNOTATION_DIR" >&2; exit 1; }
cd "$ANNOTATION_DIR"
[[ -e "$INPUT" ]] || { echo "Missing input: ${ANNOTATION_DIR}/${INPUT}; run retrieve_annotations.sh first" >&2; exit 1; }
echo "Working in ${ANNOTATION_DIR}"

for chromosome in "${CHROMOSOMES[@]}"; do
    lift_chromosome "$INPUT" "$chromosome" "GRCh37-cCREs_${chromosome}.bed"
done
