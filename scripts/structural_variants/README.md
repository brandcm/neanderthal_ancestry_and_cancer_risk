# `Structural Variants`

- `analyze_PTHLH_insertion.py` — Characterizes a polymorphic LINE-1 insertion near *PTHLH* from long-read genome assemblies (578 haplotypes; personal communication). Classifies haplotypes by size, locates breakpoints with k-mer anchoring, extracts the inserted sequence, identifies target site duplication and orientation. Input data (personal communication) are not redistributed. Usage:

  ```
  python3 scripts/structural_variants/analyze_PTHLH_insertion.py \
      --fasta-dir data/structural_variants/PTHLH_long_reads/fasta \
      --reference data/structural_variants/PTHLH_long_reads/reference/PTHLH.t2t.fa \
      --outdir data/structural_variants \
      --exemplar HG00512_h1
  ```