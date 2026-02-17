#!/bin/bash
module load Kent_tools

GTF="gencode.v49.annotation.gtf"
OUT="refFlat.txt"

# 1. Use -geneNameAsName2 to put "CDR1" in the GeneName column
# This handles the coordinate math internally (no more 0-length exons)
gtfToGenePred -genePredExt -geneNameAsName2 "$GTF" tmp.gp

# 2. Reorder to refFlat format:
# Col 12: GeneName (CDR1)
# Col 1:  TranscriptID (ENST...)
# Col 2-10: The rest of the coordinates
awk 'BEGIN{OFS="\t"} {print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' tmp.gp > refFlat.tmp

# 3. Final cleanup: Ensure every coordinate in the list ends with a comma
# CIRCexplorer2 is extremely strict about the trailing comma.
awk -F'\t' 'BEGIN{OFS="\t"} {
    gsub(/,+$/, "", $10); $10=$10","; 
    gsub(/,+$/, "", $11); $11=$11","; 
    print
}' refFlat.tmp > "$OUT"

rm tmp.gp refFlat.tmp

# VERIFICATION
echo "--- Verification of CDR1 locus ---"
grep "ENST00000674533" "$OUT"
