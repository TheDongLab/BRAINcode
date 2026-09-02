#!/bin/bash
#SBATCH --job-name=GAPDH_strand_test
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/GAPDH_strand_test.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/GAPDH_strand_test.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

module load BEDTools
module load SAMtools

python - <<'PY'
import re,subprocess,tempfile,os
from pathlib import Path

# ============================================================
# PARAMETERS
# ============================================================

SD=Path("/home/zw529/donglab/data/target_ALS/Cerebellum/RNAseq/Processed/CGND_HRA_00070")
GTF=Path("/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf")
BAM=SD/"STAR.Aligned.sortedByCoord.out.bam"
NORM=SD/"normalization.tab"

GENE_NAME="GAPDH"

# ============================================================
# 1. Resolve GAPDH from GENCODE
# ============================================================

gene=None
chrom=start=end=strand=None

with open(GTF) as f:
    for l in f:
        if l.startswith("#"):
            continue
        x=l.rstrip().split("\t")
        if len(x)<9 or x[2]!="gene":
            continue

        name=re.search(r'gene_name "([^"]+)"',x[8])
        gid=re.search(r'gene_id "([^"]+)"',x[8])

        if name and name.group(1)==GENE_NAME:
            gene=gid.group(1)
            chrom=x[0]
            start=int(x[3])-1
            end=int(x[4])
            strand=x[6]
            break

if gene is None:
    raise RuntimeError("GAPDH not found in GENCODE GTF")

print(f"GENCODE GAPDH: {gene} | {chrom}:{start+1}-{end} | strand={strand}")

# ============================================================
# 2. Get GAPDH HTSeq raw count from normalization.tab
# ============================================================

raw_count=None

with open(NORM) as f:
    header=f.readline().rstrip().split("\t")
    gi=header.index("gene_id")
    ci=header.index("raw_count")

    for l in f:
        x=l.rstrip().split("\t")
        if len(x)<=max(gi,ci):
            continue
        if x[gi]==gene:
            raw_count=int(float(x[ci]))
            break

if raw_count is None:
    raise RuntimeError(f"{gene} absent from normalization.tab")

# ============================================================
# 3. Build merged GAPDH exon BED
# ============================================================

exons=[]

with open(GTF) as f:
    for l in f:
        if l.startswith("#"):
            continue
        x=l.rstrip().split("\t")
        if len(x)<9 or x[2]!="exon":
            continue

        if re.search(r'gene_id "'+re.escape(gene)+r'"',x[8]):
            exons.append((x[0],int(x[3])-1,int(x[4])))

if not exons:
    raise RuntimeError("No GAPDH exons found")

exons.sort()
merged=[]

for c,s,e in exons:
    if merged and c==merged[-1][0] and s<=merged[-1][2]:
        merged[-1][2]=max(merged[-1][2],e)
    else:
        merged.append([c,s,e])

with tempfile.NamedTemporaryFile("w",delete=False,suffix=".bed") as tmp:
    bed=tmp.name
    for c,s,e in merged:
        tmp.write(f"{c}\t{s}\t{e}\n")

# ============================================================
# 4. Count four paired-end orientation classes
#
# Remove:
#   0x100 = secondary
#   0x800 = supplementary
# => 0x900
# ============================================================

def count(required,excluded=0x900):
    cmd=[
        "samtools","view","-c",
        "-f",hex(required),
        "-F",hex(excluded),
        "-L",bed,
        str(BAM)
    ]
    return int(subprocess.check_output(cmd,text=True).strip())

# R1 forward = read1 flag, NOT reverse
r1_fwd=count(0x40,0x910)

# R1 reverse = read1 + reverse
r1_rev=count(0x50,0x900)

# R2 forward = read2 flag, NOT reverse
r2_fwd=count(0x80,0x910)

# R2 reverse = read2 + reverse
r2_rev=count(0x90,0x900)

os.unlink(bed)

# ============================================================
# 5. Evaluate strand orientation and bam2bigwig settings
# ============================================================

# For fr-firststrand / reverse-stranded paired-end RNA-seq:
#
#   PLUS transcript strand:
#     R1 reverse + R2 forward
#
#   MINUS transcript strand:
#     R1 forward + R2 reverse
#
# This corresponds to:
#
#   PLUS  = R2 forward + R1 reverse
#   MINUS = R1 forward + R2 reverse

plus_orientation  = r1_rev + r2_fwd
minus_orientation = r1_fwd + r2_rev

print()
print("="*72)
print("TARGETALS STRAND-ORIENTATION TEST — GAPDH")
print("="*72)
print(f"Sample        : {SD.name}")
print(f"Gene          : {GENE_NAME} ({gene})")
print(f"Coordinates   : {chrom}:{start+1}-{end}")
print(f"GENCODE strand: {strand}")
print(f"HTSeq count   : {raw_count}")
print()

print("Observed alignments overlapping GAPDH exons:")
print(f"  R1 forward  : {r1_fwd}")
print(f"  R1 reverse  : {r1_rev}")
print(f"  R2 forward  : {r2_fwd}")
print(f"  R2 reverse  : {r2_rev}")
print()

print("fr-firststrand / reverse-stranded interpretation:")
print(f"  PLUS  transcript orientation = R1 reverse + R2 forward = {plus_orientation}")
print(f"  MINUS transcript orientation = R1 forward + R2 reverse = {minus_orientation}")
print()

if strand == "+":
    compatible = plus_orientation
    opposite   = minus_orientation
    expected   = "R1 reverse + R2 forward"
elif strand == "-":
    compatible = minus_orientation
    opposite   = plus_orientation
    expected   = "R1 forward + R2 reverse"
else:
    compatible = opposite = 0
    expected = "UNKNOWN"

print(f"For this GENCODE '{strand}' transcript, expected orientation is:")
print(f"  {expected}")
print()
print(f"Transcript-compatible reads = {compatible}")
print(f"Opposite-orientation reads  = {opposite}")
print(f"Orientation ratio           = {compatible/max(opposite,1):.2f}x")
print()

print("Interpretation:")
if compatible > opposite:
    print(f"  PASS: {GENE_NAME} shows the expected {strand}-strand orientation.")
    print("  This is consistent with fr-firststrand / reverse-stranded RNA-seq.")
    print("  bam2bigwig.sh should use:")
    print("    PLUS  = R2 forward + R1 reverse")
    print("      R2 forward : samtools view -bh -f 0x80 -F 0x10")
    print("      R1 reverse : samtools view -bh -f 0x40 -f 0x10")
    print("    MINUS = R1 forward + R2 reverse")
    print("      R1 forward : samtools view -bh -f 0x40 -F 0x10")
    print("      R2 reverse : samtools view -bh -f 0x80 -f 0x10")
else:
    print(f"  WARNING: {GENE_NAME} does not show the expected {strand}-strand orientation.")

print("="*72)
PY
