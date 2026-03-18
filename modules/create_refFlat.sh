#!/bin/bash

set -euo pipefail
module load Kent_tools

GTF="gencode.v49.annotation.gtf"
PREFIX="${GTF%.gtf}"
OUT="refFlat.txt"
OUT_CIRS7="${OUT%.txt}_cIRS7.txt"

# 1. Build transcript_id -> label lookup from GTF
awk '$3=="transcript"{
    match($0,/transcript_id "([^"]+)"/,tid)
    match($0,/gene_id "([^"]+)"/,gid)
    match($0,/transcript_type "([^"]+)"/,tt)
    match($0,/gene_name "([^"]+)"/,gn)
    print tid[1], gid[1]"___"tt[1]"___"gn[1]
}' "$GTF" > tmp.label_lookup.tsv

# 2. GTF -> GenePred
gtfToGenePred -genePredExt -geneNameAsName2 "$GTF" tmp.gp

# 3. GenePred -> BED12 (v49, original unmodified name col)
genePredToBed tmp.gp "${PREFIX}.transcript.bed12"

# 4. BED12 with rich name col (ENST___ENSG___transcript_type___gene_name)
awk 'BEGIN{OFS="\t"}
NR==FNR{ lut[$1]=$2; next }
{ $4 = ($4 in lut) ? $4"___"lut[$4] : $4; print }
' tmp.label_lookup.tsv "${PREFIX}.transcript.bed12" > "${PREFIX}.labeled.transcript.bed12"

# 5. GenePred -> refFlat column order
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

# 8. Generate labeled BED12 for v47
echo ""
echo "--- Generating v47 labeled BED12 ---"
awk '$3=="transcript"{
    match($0,/transcript_id "([^"]+)"/,tid)
    match($0,/gene_id "([^"]+)"/,gid)
    match($0,/transcript_type "([^"]+)"/,tt)
    match($0,/gene_name "([^"]+)"/,gn)
    print tid[1], gid[1]"___"tt[1]"___"gn[1]
}' gencode.v47.annotation.gtf > tmp.v47.label_lookup.tsv

gtfToGenePred -genePredExt -geneNameAsName2 gencode.v47.annotation.gtf tmp.v47.gp
genePredToBed tmp.v47.gp gencode.v47.annotation.transcript.bed12

awk 'BEGIN{OFS="\t"}
NR==FNR{ lut[$1]=$2; next }
{ $4 = ($4 in lut) ? $4"___"lut[$4] : $4; print }
' tmp.v47.label_lookup.tsv gencode.v47.annotation.transcript.bed12 > gencode.v47.annotation.labeled.transcript.bed12

rm tmp.v47.gp tmp.v47.label_lookup.tsv gencode.v47.annotation.transcript.bed12

echo "--- Done ---"
