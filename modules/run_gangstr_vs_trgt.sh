#!/bin/bash
#SBATCH --job-name=GangSTR_vs_TRGT
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --partition=week
#SBATCH --output=/home/zw529/donglab/data/target_ALS/WGS_LR/repeat_comparison_gangSTR_vs_TRGT/logs/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/WGS_LR/repeat_comparison_gangSTR_vs_TRGT/logs/%x_%j.err

set -euo pipefail

###############################################################################
# PARAMETERS
###############################################################################

BASE="$HOME/donglab/data/target_ALS/WGS_LR"
OUT="$BASE/repeat_comparison_gangSTR_vs_TRGT"

CATALOG="$BASE/str_screen_20260911/inputs/hg38_ver13.bed.gz"

REF="$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"

SHORT_BAM_ORIG="$BASE/NEUAD700YFB.SD-029-24-CBLL.2.bam"
SHORT_BAI_ORIG="$BASE/wgs_bam_NEUAD700YFB_NEUAD700YFB.SD-029-24-CBLL.2.bam.bai"

LONG_BAM_ORIG="$BASE/NEUAD700YFB.SD-029-24-CBLL.3.mapped.bam"
LONG_BAI_ORIG="$BASE/lr-wgs_bam_NEUAD700YFB_NEUAD700YFB.SD-029-24-CBLL.3.mapped.bam.bai"

GANGSTR_ENV="$HOME/donglab/pipelines/modules/gangstr-env"
TRGT="$HOME/donglab/pipelines/modules/trgt-5.1.0/trgt"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
FLANK=250
SEED=20260911

INPUTS="$OUT/inputs"
GANGSTR_OUT="$OUT/gangstr"
TRGT_OUT="$OUT/trgt"
SUMMARY="$OUT/summary"

MASTER="$INPUTS/master.tsv"
GANGSTR_BED="$INPUTS/master.gangstr.bed"
TRGT_BED="$INPUTS/master.trgt.bed"

mkdir -p "$INPUTS" "$GANGSTR_OUT" "$TRGT_OUT" "$SUMMARY" "$OUT/logs"

###############################################################################
# SOFTWARE
###############################################################################

module purge
module load SAMtools

PYTHON="$HOME/donglab/pipelines/modules/miniconda3/bin/python"
GANGSTR="$HOME/donglab/pipelines/modules/gangstr-env/bin/GangSTR"
TRGT="$HOME/donglab/pipelines/modules/trgt-5.1.0/trgt"
SAMTOOLS="$(command -v samtools)"

echo "============================================================"
echo "SOFTWARE"
echo "============================================================"

echo "GangSTR: $GANGSTR"
"$GANGSTR" --version

echo
echo "TRGT: $TRGT"
"$TRGT" --version

echo
echo "samtools: $SAMTOOLS"
"$SAMTOOLS" --version | head -2

echo
echo "Python: $PYTHON"
"$PYTHON" -c 'import sys,pysam; print(sys.executable); print("pysam",pysam.__version__)'

###############################################################################
# INPUT CHECKS
###############################################################################

echo "============================================================"
echo "INPUT CHECKS"
echo "============================================================"

for f in \
    "$CATALOG" \
    "$REF" \
    "$REF.fai" \
    "$SHORT_BAM_ORIG" \
    "$SHORT_BAI_ORIG" \
    "$LONG_BAM_ORIG" \
    "$LONG_BAI_ORIG"
do
    [[ -s "$f" ]] || {
        echo "ERROR: Missing/empty input: $f" >&2
        exit 1
    }
    echo "OK: $f"
done

[[ -x "$TRGT" ]] || {
    echo "ERROR: TRGT is not executable: $TRGT" >&2
    exit 1
}

###############################################################################
# CREATE LOCAL SYMLINKS WITH STANDARD BAM/BAI NAMES
###############################################################################

ln -sfn "$SHORT_BAM_ORIG" "$INPUTS/short.bam"
ln -sfn "$SHORT_BAI_ORIG" "$INPUTS/short.bam.bai"

ln -sfn "$LONG_BAM_ORIG" "$INPUTS/long.bam"
ln -sfn "$LONG_BAI_ORIG" "$INPUTS/long.bam.bai"

ln -sfn "$REF" "$INPUTS/hg38.fa"
ln -sfn "$REF.fai" "$INPUTS/hg38.fa.fai"

SHORT_BAM="$INPUTS/short.bam"
LONG_BAM="$INPUTS/long.bam"
RUN_REF="$INPUTS/hg38.fa"

###############################################################################
# VERIFY PRIMARY AUTOSOME REFERENCES
###############################################################################

echo
echo "============================================================"
echo "VERIFY chr1-chr22 BAM/FASTA CONTIGS"
echo "============================================================"

"$SAMTOOLS" view -H "$SHORT_BAM" |
awk -F'\t' '
$1=="@SQ" {
    sn=""; ln="";
    for(i=2;i<=NF;i++) {
        if($i~/^SN:/) {sn=$i; sub(/^SN:/,"",sn)}
        if($i~/^LN:/) {ln=$i; sub(/^LN:/,"",ln)}
    }
    if(sn ~ /^chr([1-9]|1[0-9]|2[0-2])$/)
        print sn"\t"ln
}' |
sort -V > "$SUMMARY/short_bam.autosomes.tsv"

"$SAMTOOLS" view -H "$LONG_BAM" |
awk -F'\t' '
$1=="@SQ" {
    sn=""; ln="";
    for(i=2;i<=NF;i++) {
        if($i~/^SN:/) {sn=$i; sub(/^SN:/,"",sn)}
        if($i~/^LN:/) {ln=$i; sub(/^LN:/,"",ln)}
    }
    if(sn ~ /^chr([1-9]|1[0-9]|2[0-2])$/)
        print sn"\t"ln
}' |
sort -V > "$SUMMARY/long_bam.autosomes.tsv"

awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/{print $1"\t"$2}' "$REF.fai" |
sort -V > "$SUMMARY/reference.autosomes.tsv"

if ! diff -u \
    "$SUMMARY/reference.autosomes.tsv" \
    "$SUMMARY/short_bam.autosomes.tsv" \
    > "$SUMMARY/reference_vs_short.diff"
then
    echo "ERROR: short-read BAM chr1-chr22 do not match reference." >&2
    cat "$SUMMARY/reference_vs_short.diff" >&2
    exit 1
fi

if ! diff -u \
    "$SUMMARY/reference.autosomes.tsv" \
    "$SUMMARY/long_bam.autosomes.tsv" \
    > "$SUMMARY/reference_vs_long.diff"
then
    echo "ERROR: long-read BAM chr1-chr22 do not match reference." >&2
    cat "$SUMMARY/reference_vs_long.diff" >&2
    exit 1
fi

echo "PASS: short BAM, long BAM, and reference chr1-chr22 names/lengths match."

###############################################################################
# BUILD MASTER CATALOG
#
# IMPORTANT:
# GangSTR hg38 v13 source coordinates are treated exactly as in the previous
# validated script:
#
#     GangSTR = 1-based inclusive
#     TRGT    = 0-based half-open
#
# Rule 4 only:
#   1. chr1-chr22
#   2. reference interval must be an exact uninterrupted repeat
#      allowing cyclic motif rotations and reverse-complement rotations
#   3. no N in repeat + 250 bp left + 250 bp right
#
# NO motif-size filter.
# NO repeat-length filter.
# NO sampling/binning.
# NO MAPQ/coverage/mappability filter.
###############################################################################

echo
echo "============================================================"
echo "BUILD SHARED MASTER CATALOG"
echo "============================================================"

export CATALOG REF MASTER GANGSTR_BED TRGT_BED FLANK SUMMARY

"$PYTHON" <<'PY'
import gzip
import json
import os
from collections import Counter
from pathlib import Path

import pysam

catalog = Path(os.environ["CATALOG"])
ref_path = Path(os.environ["REF"])
master_path = Path(os.environ["MASTER"])
gangstr_path = Path(os.environ["GANGSTR_BED"])
trgt_path = Path(os.environ["TRGT_BED"])
summary_dir = Path(os.environ["SUMMARY"])
flank = int(os.environ["FLANK"])

autosomes = {f"chr{i}" for i in range(1, 23)}
comp = str.maketrans("ACGTacgt", "TGCAtgca")

fasta = pysam.FastaFile(str(ref_path))
contig_lengths = dict(zip(fasta.references, fasta.lengths))

counts = Counter()

with gzip.open(catalog, "rt") as src, \
     master_path.open("w") as master, \
     gangstr_path.open("w") as gang, \
     trgt_path.open("w") as trgt:

    master.write(
        "locus_id\tchrom\tgangstr_start_1based\tgangstr_end_1based\t"
        "trgt_start_0based\ttrgt_end_0based\tmotif_length\tmotif\t"
        "reference_length_bp\n"
    )

    for line_no, line in enumerate(src, 1):
        if not line.strip() or line.startswith("#"):
            continue

        counts["source_records"] += 1

        fields = line.rstrip().split()
        if len(fields) != 5:
            counts["bad_column_count"] += 1
            continue

        chrom, start_s, end_s, k_s, motif = fields

        if chrom not in autosomes:
            counts["non_autosomal"] += 1
            continue

        counts["autosomal_records"] += 1

        try:
            start = int(start_s)
            end = int(end_s)
            k = int(k_s)
        except ValueError:
            counts["invalid_numeric_field"] += 1
            continue

        motif = motif.upper()

        if chrom not in contig_lengths:
            counts["missing_reference_contig"] += 1
            continue

        chrom_len = contig_lengths[chrom]

        # GangSTR v13 coordinates: 1-based inclusive.
        if start < 1 or end < start or end > chrom_len:
            counts["invalid_coordinates"] += 1
            continue

        reference_length = end - start + 1

        # Exact repeat purity inherently requires whole motif copies.
        if k <= 0 or len(motif) != k:
            counts["motif_length_mismatch"] += 1
            continue

        if reference_length % k != 0:
            counts["not_whole_motif_copies"] += 1
            continue

        # Require full 250 bp of sequence on both sides.
        if start - flank < 1 or end + flank > chrom_len:
            counts["insufficient_250bp_flank"] += 1
            continue

        # pysam.fetch = 0-based half-open.
        seq = fasta.fetch(chrom, start - 1, end).upper()

        rev = motif.translate(comp)[::-1]
        rotations = {
            m[i:] + m[:i]
            for m in (motif, rev)
            for i in range(k)
        }

        if not any(seq == m * (len(seq) // k) for m in rotations):
            counts["imperfect_reference"] += 1
            continue

        # Exactly 250 bp left + repeat + exactly 250 bp right.
        context = fasta.fetch(
            chrom,
            start - 1 - flank,
            end + flank
        ).upper()

        expected_context_len = flank + reference_length + flank
        if len(context) != expected_context_len:
            counts["unexpected_context_length"] += 1
            continue

        if "N" in context:
            counts["N_flank"] += 1
            continue

        trgt_start = start - 1
        trgt_end = end

        locus_id = f"{chrom}_{start}_{end}_{motif}"

        master.write(
            f"{locus_id}\t{chrom}\t{start}\t{end}\t"
            f"{trgt_start}\t{trgt_end}\t{k}\t{motif}\t"
            f"{reference_length}\n"
        )

        # Preserve original GangSTR coordinate convention exactly.
        gang.write(
            f"{chrom}\t{start}\t{end}\t{k}\t{motif}\n"
        )

        # Exact same biological interval converted to proper TRGT BED coords.
        trgt.write(
            f"{chrom}\t{trgt_start}\t{trgt_end}\t"
            f"ID={locus_id};MOTIFS={motif};STRUC=<TR>\n"
        )

        counts["master_pass"] += 1

fasta.close()

summary = {
    "source_catalog": str(catalog),
    "reference": str(ref_path),
    "rule": (
        "autosomes chr1-chr22; exact uninterrupted reference repeat "
        "allowing cyclic rotations and reverse-complement rotations; "
        f"no N in repeat plus {flank} bp on each side"
    ),
    "coordinate_conventions": {
        "GangSTR": "1-based inclusive, preserved from hg38_ver13 catalog",
        "TRGT": "0-based half-open; start = GangSTR start - 1, end unchanged",
    },
    "counts": dict(counts),
}

(summary_dir / "catalog_filter_summary.json").write_text(
    json.dumps(summary, indent=2) + "\n"
)

print(json.dumps(summary, indent=2))
PY

###############################################################################
# VERIFY DERIVED CATALOGS ARE EXACTLY 1:1
###############################################################################

echo
echo "============================================================"
echo "VERIFY SHARED LOCUS SET"
echo "============================================================"

MASTER_N=$(( $(wc -l < "$MASTER") - 1 ))
GANGSTR_N=$(wc -l < "$GANGSTR_BED")
TRGT_N=$(wc -l < "$TRGT_BED")

echo "Master loci : $MASTER_N"
echo "GangSTR loci: $GANGSTR_N"
echo "TRGT loci   : $TRGT_N"

if [[ "$MASTER_N" -ne "$GANGSTR_N" || "$MASTER_N" -ne "$TRGT_N" ]]; then
    echo "ERROR: derived catalog counts differ." >&2
    exit 1
fi

# Verify every GangSTR row corresponds exactly to its master row.
awk 'BEGIN{FS=OFS="\t"} NR>1{print $2,$3,$4,$7,$8}' "$MASTER" \
    > "$SUMMARY/master_as_gangstr.bed"

if ! cmp -s "$SUMMARY/master_as_gangstr.bed" "$GANGSTR_BED"; then
    echo "ERROR: GangSTR BED differs from master TSV." >&2
    exit 1
fi

# Verify every TRGT interval/ID is generated from same master row.
"$PYTHON" <<'PY'
import os

master = os.environ["MASTER"]
trgt = os.environ["TRGT_BED"]

with open(master) as m, open(trgt) as t:
    header = next(m)

    for i, (mr, tr) in enumerate(zip(m, t), 1):
        x = mr.rstrip().split("\t")
        y = tr.rstrip().split("\t")

        locus_id, chrom, gs, ge, ts, te, k, motif, reflen = x

        expected = [
            chrom,
            ts,
            te,
            f"ID={locus_id};MOTIFS={motif};STRUC=<TR>"
        ]

        if y != expected:
            raise SystemExit(
                f"TRGT/master mismatch at record {i}\n"
                f"MASTER={x}\nTRGT={y}\nEXPECTED={expected}"
            )

print("PASS: every TRGT record is an exact deterministic conversion of master.tsv")
PY

echo "PASS: GangSTR and TRGT will interrogate the same $MASTER_N loci."

###############################################################################
# SAVE MANIFEST BEFORE CALLING
###############################################################################

echo
echo "============================================================"
echo "MANIFEST"
echo "============================================================"

{
    echo -e "item\tpath\tsha256"
    for f in "$CATALOG" "$REF" "$MASTER" "$GANGSTR_BED" "$TRGT_BED"; do
        printf "%s\t%s\t" "$(basename "$f")" "$f"
        sha256sum "$f" | awk '{print $1}'
    done
} > "$SUMMARY/input_manifest.tsv"

cat "$SUMMARY/input_manifest.tsv"

###############################################################################
# RUN GANGSTR — SHORT READS
###############################################################################

echo
echo "============================================================"
echo "RUN GANGSTR"
echo "Start: $(date)"
echo "============================================================"

"$GANGSTR" \
    --bam "$SHORT_BAM" \
    --ref "$RUN_REF" \
    --regions "$GANGSTR_BED" \
    --out "$GANGSTR_OUT/NEUAD700YFB.gangstr" \
    --seed "$SEED"

echo "GangSTR completed: $(date)"

###############################################################################
# RUN TRGT — LONG READS
###############################################################################

echo
echo "============================================================"
echo "RUN TRGT"
echo "Start: $(date)"
echo "============================================================"

"$TRGT" genotype \
    --genome "$RUN_REF" \
    --reads "$LONG_BAM" \
    --repeats "$TRGT_BED" \
    --output-prefix "$TRGT_OUT/NEUAD700YFB.trgt" \
    --preset wgs \
    --threads "$THREADS"

echo "TRGT completed: $(date)"

###############################################################################
# BASIC OUTPUT AUDIT
###############################################################################

echo
echo "============================================================"
echo "OUTPUT AUDIT"
echo "============================================================"

find "$GANGSTR_OUT" "$TRGT_OUT" \
    -maxdepth 1 \
    -type f \
    -printf '%p\t%s bytes\n' \
    | sort \
    > "$SUMMARY/output_files.tsv"

cat "$SUMMARY/output_files.tsv"

# Count non-header GangSTR VCF records if generated.
GVCF="$GANGSTR_OUT/NEUAD700YFB.gangstr.vcf"

if [[ -s "$GVCF" ]]; then
    grep -vc '^#' "$GVCF" \
        > "$SUMMARY/gangstr_vcf_record_count.txt"
    echo "GangSTR VCF records: $(cat "$SUMMARY/gangstr_vcf_record_count.txt")"
else
    echo "WARNING: expected GangSTR VCF not found at $GVCF" >&2
fi

# Count TRGT records.
TVCF="$TRGT_OUT/NEUAD700YFB.trgt.vcf.gz"

if [[ -s "$TVCF" ]]; then
    zgrep -vc '^#' "$TVCF" \
        > "$SUMMARY/trgt_vcf_record_count.txt"
    echo "TRGT VCF records: $(cat "$SUMMARY/trgt_vcf_record_count.txt")"
else
    echo "WARNING: expected TRGT VCF not found at $TVCF" >&2
fi

###############################################################################
# FINAL SUMMARY
###############################################################################

echo
echo "============================================================"
echo "DONE"
echo "============================================================"
echo "Shared input loci : $MASTER_N"
echo "Master TSV        : $MASTER"
echo "GangSTR BED       : $GANGSTR_BED"
echo "TRGT BED          : $TRGT_BED"
echo "GangSTR outputs   : $GANGSTR_OUT"
echo "TRGT outputs      : $TRGT_OUT"
echo "Summary           : $SUMMARY"
echo "Slurm logs        : $OUT/logs"
echo "Finished          : $(date)"
