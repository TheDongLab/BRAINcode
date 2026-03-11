#!/bin/bash

set -euo pipefail
module load Kent_tools

GTF="${1:-gencode.v49.annotation.gtf}"
OUT="${2:-refFlat.txt}"
OUT_CIRS7="${OUT%.txt}_cIRS7.txt"
PREFIX="${GTF%.gtf}"

# 1. Build transcript_id -> label lookup from GTF
# Label format: ENSG___transcript_type___gene_name
awk '$3=="transcript"{
    match($0,/transcript_id "([^"]+)"/,tid)
    match($0,/gene_id "([^"]+)"/,gid)
    match($0,/transcript_type "([^"]+)"/,tt)
    match($0,/gene_name "([^"]+)"/,gn)
    print tid[1], gid[1]"___"tt[1]"___"gn[1]
}' "$GTF" > tmp.label_lookup.tsv

# 2. GTF -> GenePred
gtfToGenePred -genePredExt -geneNameAsName2 "$GTF" tmp.gp

# 3. GenePred -> BED12 (original, unmodified name col)
genePredToBed tmp.gp "${PREFIX}.bed12"

# 4. BED12 with rich name col (ENST___ENSG___transcript_type___gene_name)
awk 'BEGIN{OFS="\t"}
NR==FNR{ lut[$1]=$2; next }
{ $4 = ($4 in lut) ? $4"___"lut[$4] : $4; print }
' tmp.label_lookup.tsv "${PREFIX}.bed12" > "${PREFIX}.labeled.bed12"

# 5. GenePred -> refFlat column order (standard, no label modification)
awk 'BEGIN{OFS="\t"} {print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' tmp.gp > refFlat.tmp

# 6. Enforce trailing commas on exon fields (CIRCexplorer2 requirement)
awk -F'\t' 'BEGIN{OFS="\t"} {
    gsub(/,+$/, "", $10); $10=$10",";
    gsub(/,+$/, "", $11); $11=$11",";
    print
}' refFlat.tmp > "$OUT"

# 7. Create ciRS-7 patched version
cp "$OUT" "$OUT_CIRS7"
printf "CDR1as\tciRS-7\tchrX\t-\t140783174\t140784660\t140783174\t140784660\t1\t140783174,\t140784660,\n" >> "$OUT_CIRS7"

rm tmp.gp tmp.label_lookup.tsv refFlat.tmp

# VERIFICATION
echo "--- Verification ---"
grep "ENST00000674533" "$OUT"      || echo "WARNING: Linear CDR1 not found in $OUT"
grep "CDR1as"         "$OUT_CIRS7" || echo "WARNING: Circular CDR1as not found in $OUT_CIRS7"
