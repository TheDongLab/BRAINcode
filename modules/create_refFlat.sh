#!/bin/bash
module load Kent_tools

GTF="gencode.v49.annotation.gtf"
OUT="refFlat_cIRS7.txt"

# 1. Convert GTF to GenePred with GeneName as Name2
gtfToGenePred -genePredExt -geneNameAsName2 "$GTF" tmp.gp

# 2. Reorder to refFlat format:
# GeneName (12), TranscriptID (1), Chrom (2), Strand (3), TxStart (4), TxEnd (5), 
# CdsStart (6), CdsEnd (7), ExonCount (8), ExonStarts (9), ExonEnds (10)
awk 'BEGIN{OFS="\t"} {print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' tmp.gp > refFlat.tmp

# 3. Cleanup: Ensure every coordinate list ends with a comma (strict CIRCexplorer2 requirement)
awk -F'\t' 'BEGIN{OFS="\t"} {
    gsub(/,+$/, "", $10); $10=$10","; 
    gsub(/,+$/, "", $11); $11=$11","; 
    print
}' refFlat.tmp > "$OUT"

# 4. PATCH: Manually append the ciRS-7 (CDR1as) entry with literal tabs
printf "CDR1as\tciRS-7\tchrX\t-\t140783174\t140784660\t140783174\t140784660\t1\t140783174,\t140784660,\n" >> "$OUT"

rm tmp.gp refFlat.tmp

# VERIFICATION
echo "--- Verification of CDR1 locus (Linear vs Circular) ---"
# Check for the standard linear transcript (CDR1)
grep "ENST00000674533" "$OUT" || echo "WARNING: Linear CDR1 not found in $OUT"

# Check for the newly added circular entry (CDR1as)
grep "CDR1as" "$OUT" || echo "WARNING: Circular CDR1as not found in $OUT"
