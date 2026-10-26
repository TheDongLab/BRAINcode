#!/bin/bash
#############################################
# Strand-specific BigWig generation.
# Default reference supports existing TargetALS calls.
# An optional third argument supplies AnswerALS chromosome sizes and restricts analysis to the contigs listed in that file.
#############################################

set -euo pipefail

DEFAULT_CHROMINFO="$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/STAR/genome.sizes"

###########################################
# Inputs
###########################################
inputfile="${1:-}"
split="${2:--nosplit}"
CHROMINFO="${3:-$DEFAULT_CHROMINFO}"

# Supplying a third argument enables strict contig-whitelist mode.
CANONICAL_ONLY=false
if [[ $# -ge 3 ]]; then
    CANONICAL_ONLY=true
fi

if [[ -z "$inputfile" || ! -e "$inputfile" ]]; then
    echo "Usage: bam2bigwig.sh <in.bam|cram> <-split|-nosplit> [genome.sizes]" >&2
    exit 1
fi

if [[ ! -s "$CHROMINFO" ]]; then
    echo "[ERROR] Chromosome sizes file is missing or empty: $CHROMINFO" >&2
    exit 1
fi

ext="${inputfile##*.}"
bname="${inputfile%.*}"

case "$ext" in
    bam|cram)
        echo "Input is ${ext^^}. Generating normalized BigWig..."
        ;;
    *)
        echo "[ERROR] Unsupported input format: $ext" >&2
        exit 1
        ;;
esac

case "$split" in
    -split|-nosplit)
        ;;
    *)
        echo "[ERROR] Second argument must be -split or -nosplit." >&2
        exit 1
        ;;
esac

###########################################
# Temporary directory
###########################################
WORK_TMPDIR=$(mktemp -d "/tmp/bam2bigwig_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-$$}_XXXXXX")
trap 'rm -rf "$WORK_TMPDIR"' EXIT

###########################################
# Read the BAM/CRAM reference dictionary
###########################################
BAM_REFS="$WORK_TMPDIR/bam.references.tsv"

samtools view -H "$inputfile" |
    awk -F'\t' '$1 == "@SQ" {
        sn=""; ln="";
        for (i=2; i<=NF; i++) {
            if ($i ~ /^SN:/) sn=substr($i,4);
            if ($i ~ /^LN:/) ln=substr($i,4);
        }
        if (sn != "" && ln != "") print sn "\t" ln;
    }' > "$BAM_REFS"

if [[ ! -s "$BAM_REFS" ]]; then
    echo "[ERROR] Could not read any reference contigs from $inputfile" >&2
    exit 1
fi

###########################################
# Validate references and build whitelist BED
###########################################
FILTER_ARGS=()
CANONICAL_BED=""

if "$CANONICAL_ONLY"; then
    # Region-based selection requires an index.
    if [[ "$ext" == "bam" ]]; then
        if [[ ! -s "${inputfile}.bai" && ! -s "${bname}.bai" ]]; then
            echo "[ERROR] Canonical-contig filtering requires a BAM index for $inputfile" >&2
            exit 1
        fi
    else
        if [[ ! -s "${inputfile}.crai" && ! -s "${bname}.crai" ]]; then
            echo "[ERROR] Canonical-contig filtering requires a CRAM index for $inputfile" >&2
            exit 1
        fi
    fi

    # Every contig in the supplied whitelist must exist in the BAM/CRAM
    # with the same length. Extra BAM/CRAM contigs are allowed and excluded.
    while IFS=$'\t ' read -r chrom length rest; do
        [[ -n "$chrom" && -n "$length" ]] || continue

        if [[ ! "$length" =~ ^[0-9]+$ ]]; then
            echo "[ERROR] Invalid chromosome length for '$chrom' in $CHROMINFO" >&2
            exit 1
        fi

        bam_length=$(
            awk -F'\t' -v chrom="$chrom" '$1 == chrom { print $2; exit }' "$BAM_REFS"
        )

        if [[ -z "$bam_length" ]]; then
            echo "[ERROR] Whitelisted contig '$chrom' is absent from $inputfile" >&2
            exit 1
        fi

        if [[ "$bam_length" != "$length" ]]; then
            echo "[ERROR] Length mismatch for '$chrom': BAM/CRAM=$bam_length, sizes=$length" >&2
            exit 1
        fi
    done < "$CHROMINFO"

    CANONICAL_BED="$WORK_TMPDIR/canonical_contigs.bed"
    awk 'BEGIN { FS="[ \t]+"; OFS="\t" } NF >= 2 { print $1, 0, $2 }' \
        "$CHROMINFO" > "$CANONICAL_BED"

    if [[ ! -s "$CANONICAL_BED" ]]; then
        echo "[ERROR] Could not build a contig whitelist from $CHROMINFO" >&2
        exit 1
    fi

    FILTER_ARGS=(-L "$CANONICAL_BED")
    echo "Contig mode         : whitelist ($CHROMINFO)"
else
    # Preserve the existing TargetALS check: every BAM/CRAM contig must
    # appear in the default chromosome sizes file.
    while IFS=$'\t' read -r chrom length; do
        [[ -n "$chrom" && -n "$length" ]] || continue
        if ! awk -F'\t' -v chrom="$chrom" '$1 == chrom { found=1; exit } END { exit !found }' "$CHROMINFO"; then
            echo "[ERROR] BAM/CRAM contig '$chrom' is absent from $CHROMINFO" >&2
            exit 1
        fi
    done < "$BAM_REFS"

    echo "Contig mode         : existing default behavior"
fi

###########################################
# RPM normalization: primary alignments only
###########################################
PRIMARY_READS=$(samtools view "${FILTER_ARGS[@]}" -F 0x100 -c "$inputfile")

if [[ "$PRIMARY_READS" -eq 0 ]]; then
    echo "[ERROR] No primary alignments found in the selected contigs of $inputfile" >&2
    exit 1
fi

RPMscale=$(bc <<< "scale=6;1000000/$PRIMARY_READS")

echo "Primary alignments : $PRIMARY_READS"
echo "RPM scale          : $RPMscale"
echo "Chromosome sizes   : $CHROMINFO"

###########################################
# Strand-specific BigWigs
#
# For reverse-stranded / fr-firststrand libraries:
#   PLUS  = R2 forward + R1 reverse
#   MINUS = R1 forward + R2 reverse
###########################################
if [[ "$split" == "-split" ]]; then
    echo "Generating strand-specific BigWigs (reverse-stranded / fr-firststrand)."

    echo "[PLUS] Extracting R2 forward reads..."
    samtools view -bh "${FILTER_ARGS[@]}" -f 0x80 -F 0x10 "$inputfile" \
        > "$WORK_TMPDIR/plus.part1.bam"

    echo "[PLUS] Extracting R1 reverse reads..."
    samtools view -bh "${FILTER_ARGS[@]}" -f 0x40 -f 0x10 "$inputfile" \
        > "$WORK_TMPDIR/plus.part2.bam"

    echo "[PLUS] Merging read orientations..."
    samtools merge -f "$WORK_TMPDIR/plus.bam" \
        "$WORK_TMPDIR/plus.part1.bam" \
        "$WORK_TMPDIR/plus.part2.bam"

    echo "[PLUS] Creating bedGraph..."
    bedtools genomecov \
        -ibam "$WORK_TMPDIR/plus.bam" \
        -bg \
        -scale "$RPMscale" \
        -split |
        LC_COLLATE=C sort -k1,1 -k2,2n \
        > "$bname.plus.normalized.bedGraph"

    echo "[PLUS] Creating BigWig..."
    bedGraphToBigWig \
        "$bname.plus.normalized.bedGraph" \
        "$CHROMINFO" \
        "$bname.plus.normalized.bw"

    echo "[MINUS] Extracting R1 forward reads..."
    samtools view -bh "${FILTER_ARGS[@]}" -f 0x40 -F 0x10 "$inputfile" \
        > "$WORK_TMPDIR/minus.part1.bam"

    echo "[MINUS] Extracting R2 reverse reads..."
    samtools view -bh "${FILTER_ARGS[@]}" -f 0x80 -f 0x10 "$inputfile" \
        > "$WORK_TMPDIR/minus.part2.bam"

    echo "[MINUS] Merging read orientations..."
    samtools merge -f "$WORK_TMPDIR/minus.bam" \
        "$WORK_TMPDIR/minus.part1.bam" \
        "$WORK_TMPDIR/minus.part2.bam"

    echo "[MINUS] Creating bedGraph..."
    bedtools genomecov \
        -ibam "$WORK_TMPDIR/minus.bam" \
        -bg \
        -scale "$RPMscale" \
        -split |
        LC_COLLATE=C sort -k1,1 -k2,2n \
        > "$bname.minus.normalized.bedGraph"

    echo "[MINUS] Creating BigWig..."
    bedGraphToBigWig \
        "$bname.minus.normalized.bedGraph" \
        "$CHROMINFO" \
        "$bname.minus.normalized.bw"

    echo "Strand-specific BigWigs completed:"
    echo "  $bname.plus.normalized.bw"
    echo "  $bname.minus.normalized.bw"

###########################################
# Unsplit BigWig
###########################################
else
    echo "Generating unsplit BigWig."

    if [[ "$ext" == "cram" || "$CANONICAL_ONLY" == true ]]; then
        samtools view -bh "${FILTER_ARGS[@]}" "$inputfile" |
            bedtools genomecov \
                -ibam stdin \
                -bg \
                -scale "$RPMscale" \
                -split |
            LC_COLLATE=C sort -k1,1 -k2,2n \
            > "$bname.normalized.bedGraph"
    else
        bedtools genomecov \
            -ibam "$inputfile" \
            -bg \
            -scale "$RPMscale" \
            -split |
            LC_COLLATE=C sort -k1,1 -k2,2n \
            > "$bname.normalized.bedGraph"
    fi

    bedGraphToBigWig \
        "$bname.normalized.bedGraph" \
        "$CHROMINFO" \
        "$bname.normalized.bw"

    echo "BigWig completed:"
    echo "  $bname.normalized.bw"
fi
