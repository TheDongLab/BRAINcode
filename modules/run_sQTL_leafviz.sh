#!/bin/bash
#SBATCH --job-name=sQTL_leafviz
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/leafviz/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/leafviz/%x_%j.err

set -euo pipefail

# ============================================================
# ENVIRONMENT
# ============================================================

# Initialize conda for a non-interactive SLURM shell.
eval "$(conda shell.bash hook)"
conda activate RNAseq

# Suppress normal Lmod loading/reloading messages.
module load R >/dev/null 2>&1
module load BCFtools >/dev/null 2>&1

usage() {
  cat <<'EOF'
Usage:
  sbatch run_sQTL_leafviz.sh TISSUE ANCHOR_JUNCTION SNP GENE

Arguments:
  TISSUE           Canonical tissue name:
                   Motor_Cortex
                   Cervical_Spinal_Cord
                   Lumbar_Spinal_Cord
                   Thoracic_Spinal_Cord
                   Frontal_Cortex
                   Cerebellum

  ANCHOR_JUNCTION  Stable LeafCutter PSI junction ID:
                   chr19:-:17641556-17642845

  SNP              One of:
                   rs12345
                   chr19:123456
                   19:123456
                   chr19:123456:C:G
                   19:123456:C:G

                   rsID input:
                     local dbSNP rsID
                       -> GRCh38 coordinate
                       -> Target ALS VCF variant
                       -> genotypes

                   coordinate input:
                     Target ALS VCF variant
                       -> local dbSNP
                       -> rsID

  GENE             Gene symbol, e.g. UNC13A

Examples:

  sbatch run_sQTL_leafviz.sh \
    Frontal_Cortex \
    chr19:-:17636157-17639083 \
    rs8106014 \
    UNC13A

  sbatch run_sQTL_leafviz.sh \
    Frontal_Cortex \
    chr19:-:17636157-17639083 \
    chr19:17015573:C:G \
    UNC13A
EOF
}

if [[ $# -ne 4 ]]; then
  usage
  exit 1
fi

TISSUE="$1"
ANCHOR="$2"
SNP="$3"
GENE="$4"

# ============================================================
# PATHS
# ============================================================

OUTDIR="$HOME/donglab/data/target_ALS/QTL/leafviz"
mkdir -p "$OUTDIR"

VCF="$HOME/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
META="$HOME/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
DATA="$HOME/donglab/data/target_ALS"

ANNOT="$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode"
GTF="$ANNOT/gencode.v49.annotation.gtf"

# Local dbSNP build 157, GRCh38.p14
DBSNP="$ANNOT/GCF_000001405.40.gz"
CHR_MAP="$ANNOT/ncbi_to_ucsc.txt"

# Under sbatch, BASH_SOURCE[0] points to Slurm's temporary spool copy.
# SLURM_SUBMIT_DIR points to the directory where sbatch was called.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
R_SCRIPT="$SCRIPT_DIR/plot_sQTL_leafviz.R"

for f in \
  "$VCF" \
  "$VCF.tbi" \
  "$META" \
  "$GTF" \
  "$DBSNP" \
  "$DBSNP.tbi" \
  "$CHR_MAP" \
  "$R_SCRIPT"
do
  if [[ ! -e "$f" ]]; then
    echo "ERROR: required file not found: $f" >&2
    exit 1
  fi
done

for cmd in bcftools Rscript awk sed; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: $cmd is not available in PATH." >&2
    exit 1
  fi
done

echo "Environment:"
echo "  Rscript:  $(command -v Rscript)"
echo "  bcftools: $(command -v bcftools)"
echo

# ============================================================
# TISSUE MAPPING
# ============================================================

case "$TISSUE" in
  Motor_Cortex)
    TISSUE_REGEX='Motor Cortex Lateral|Motor Cortex Medial|Lateral Motor Cortex|Medial Motor Cortex|Primary Motor Cortex L|Primary Motor Cortex M|Cortex_Motor_Unspecified|Cortex_Motor_BA4|BA4 Motor Cortex|Lateral_motor_cortex|Motor Cortex|BA4'
    ;;
  Cervical_Spinal_Cord)
    TISSUE_REGEX='Spinal_Cord_Cervical|Cervical Spinal Cord|Cervical_spinal_cord|Spinal_cord_Cervical|Cervical'
    ;;
  Lumbar_Spinal_Cord)
    TISSUE_REGEX='Lumbar Spinal Cord|Spinal_Cord_Lumbosacral|Lumbosacral_Spinal_Cord|Lumbar_spinal_cord|Lumbar|Lumbosacral'
    ;;
  Thoracic_Spinal_Cord)
    TISSUE_REGEX='Thoracic Spinal Cord|Thoracic'
    ;;
  Frontal_Cortex)
    TISSUE_REGEX='Frontal Cortex|Frontal'
    ;;
  Cerebellum)
    TISSUE_REGEX='Cerebellum'
    ;;
  *)
    echo "ERROR: unknown tissue: $TISSUE" >&2
    usage
    exit 1
    ;;
esac

# ============================================================
# HELPERS
# ============================================================

sanitize() {
  printf '%s' "$1" |
    sed 's#[/: ]#_#g; s#[^A-Za-z0-9_.+-]#_#g'
}

normalize_ucsc_chr() {
  local chr="$1"

  if [[ "$chr" == chr* ]]; then
    printf '%s\n' "$chr"
  else
    printf 'chr%s\n' "$chr"
  fi
}

# UCSC chr19 -> NC_000019.10
ucsc_to_ncbi_chr() {
  local ucsc_chr="$1"

  awk -v c="$ucsc_chr" '
    $2 == c {
      print $1
      exit
    }
  ' "$CHR_MAP"
}

# NC_000019.10 -> chr19
ncbi_to_ucsc_chr() {
  local ncbi_chr="$1"

  awk -v c="$ncbi_chr" '
    $1 == c {
      print $2
      exit
    }
  ' "$CHR_MAP"
}

# ============================================================
# rsID -> LOCAL dbSNP RECORD
#
# Returns:
# NCBI_CHR  POS  RSID  REF  ALT
#
# Example:
# NC_000019.10  17015573  rs8106014  C  G,T
# ============================================================

lookup_rsid_dbsnp() {
  local rsid="$1"

  bcftools query \
    -i "ID=\"$rsid\"" \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
    "$DBSNP" |
    head -1
}

# ============================================================
# Coordinate + allele -> rsID
#
# Input uses Target ALS / UCSC coordinates:
# chr19 17015573 C G
#
# dbSNP may have:
# NC_000019.10 17015573 rs8106014 C G,T
#
# Therefore target ALT only has to be one member of dbSNP ALT.
# ============================================================

lookup_coordinate_dbsnp() {
  local ucsc_chr="$1"
  local pos="$2"
  local ref="$3"
  local alt="$4"

  local ncbi_chr
  ncbi_chr="$(ucsc_to_ncbi_chr "$ucsc_chr")"

  if [[ -z "$ncbi_chr" ]]; then
    echo "ERROR: could not convert $ucsc_chr using $CHR_MAP" >&2
    return 1
  fi

  local region="${ncbi_chr}:${pos}-${pos}"

  bcftools query \
    -r "$region" \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
    "$DBSNP" |
    awk -F'\t' \
      -v ref="$ref" \
      -v alt="$alt" '
        $4 == ref && $3 ~ /^rs/ {
          n = split($5, a, ",")
          for (i = 1; i <= n; i++) {
            if (a[i] == alt) {
              print
              exit
            }
          }
        }
      '
}

# ============================================================
# FIND TARGET ALS VCF VARIANT FOR AN rsID
#
# dbSNP can be multiallelic:
#   rs8106014 C -> G,T
#
# Target ALS may contain only:
#   C -> G
#
# We therefore:
#   1. map rsID -> coordinate + dbSNP alleles
#   2. query all Target ALS records at that coordinate
#   3. require same REF
#   4. accept a Target ALS ALT if it occurs in the dbSNP ALT set
#
# If >1 Target ALS allele matches, stop rather than guessing.
# ============================================================

resolve_rsid_in_target_vcf() {
  local rsid="$1"

  local db_line
  db_line="$(lookup_rsid_dbsnp "$rsid")"

  if [[ -z "$db_line" ]]; then
    echo "ERROR: rsID not found in local dbSNP: $rsid" >&2
    return 1
  fi

  local db_chr db_pos db_id db_ref db_alt
  IFS=$'\t' read -r db_chr db_pos db_id db_ref db_alt <<< "$db_line"

  local ucsc_chr
  ucsc_chr="$(ncbi_to_ucsc_chr "$db_chr")"

  if [[ -z "$ucsc_chr" ]]; then
    echo "ERROR: could not map dbSNP chromosome $db_chr to UCSC." >&2
    return 1
  fi

  local region="${ucsc_chr}:${db_pos}-${db_pos}"

  echo "dbSNP mapping:"
  echo "  rsID:       $db_id"
  echo "  Coordinate: ${ucsc_chr}:${db_pos}"
  echo "  dbSNP REF:  $db_ref"
  echo "  dbSNP ALT:  $db_alt"
  echo

  local matches
  matches="$(
    bcftools query \
      -r "$region" \
      -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
      "$VCF" |
    awk -F'\t' \
      -v ref="$db_ref" \
      -v dbalt="$db_alt" '
        BEGIN {
          n = split(dbalt, allowed, ",")
          for (i = 1; i <= n; i++) {
            ok[allowed[i]] = 1
          }
        }

        $4 == ref {
          m = split($5, target_alt, ",")

          matched = 0

          for (j = 1; j <= m; j++) {
            if (target_alt[j] in ok) {
              matched = 1
            }
          }

          if (matched) {
            print
          }
        }
      '
  )"

  local n_matches
  n_matches="$(printf '%s\n' "$matches" | sed '/^$/d' | wc -l)"

  if [[ "$n_matches" -eq 0 ]]; then
    echo "ERROR: rsID mapped to ${ucsc_chr}:${db_pos}," >&2
    echo "but no REF/ALT-compatible variant was found in Target ALS VCF." >&2
    echo >&2
    echo "dbSNP:" >&2
    echo "  REF=$db_ref ALT=$db_alt" >&2
    echo >&2
    echo "Target ALS records at this position:" >&2

    bcftools query \
      -r "$region" \
      -f '  %CHROM:%POS:%REF:%ALT\n' \
      "$VCF" >&2 || true

    return 1
  fi

  if [[ "$n_matches" -gt 1 ]]; then
    echo "ERROR: multiple Target ALS variants match $rsid." >&2
    echo "Refusing to guess:" >&2
    echo >&2

    printf '%s\n' "$matches" |
      awk -F'\t' '{
        print "  " $1 ":" $2 ":" $4 ":" $5
      }' >&2

    return 1
  fi

  printf '%s\n' "$matches"
}

# ============================================================
# HEADER
# ============================================================

echo "============================================================"
echo "sQTL LeafViz"
echo "============================================================"
echo "Tissue:            $TISSUE"
echo "Anchor junction:   $ANCHOR"
echo "Requested SNP:     $SNP"
echo "Gene:              $GENE"
echo "Output directory:  $OUTDIR"
echo

# ============================================================
# LOAD TARGET ALS VCF SAMPLE ORDER
# ============================================================

mapfile -t VCF_SAMPLES < <(
  bcftools query -l "$VCF"
)

if [[ ${#VCF_SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: no samples found in Target ALS VCF." >&2
  exit 1
fi

# ============================================================
# RESOLVE SNP INPUT
# ============================================================

VARIANT_LINE=""
REQUESTED_RSID=""
RESOLVED_RSID=""

SNP_CHR=""
SNP_POS=""
SNP_REF=""
SNP_ALT=""

# ------------------------------------------------------------
# CASE 1: rsID
# ------------------------------------------------------------

if [[ "$SNP" =~ ^rs[0-9]+$ ]]; then

  REQUESTED_RSID="$SNP"
  RESOLVED_RSID="$SNP"

  echo "Input type: rsID"
  echo "Resolving through local dbSNP..."
  echo

  VARIANT_LINE="$(
    resolve_rsid_in_target_vcf "$SNP"
  )"

# ------------------------------------------------------------
# CASE 2: coordinate
# ------------------------------------------------------------

else

  echo "Input type: genomic coordinate"
  echo

  IFS=':' read -r -a SNP_PARTS <<< "$SNP"

  if [[ ${#SNP_PARTS[@]} -lt 2 ]]; then
    echo "ERROR: SNP must be one of:" >&2
    echo "  rs12345" >&2
    echo "  chr19:123456" >&2
    echo "  chr19:123456:C:G" >&2
    exit 1
  fi

  SNP_CHR="$(normalize_ucsc_chr "${SNP_PARTS[0]}")"
  SNP_POS="${SNP_PARTS[1]}"

  if ! [[ "$SNP_POS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: invalid SNP position: $SNP_POS" >&2
    exit 1
  fi

  REGION="${SNP_CHR}:${SNP_POS}-${SNP_POS}"

  # ----------------------------------------------------------
  # Coordinate + REF + ALT
  # ----------------------------------------------------------

  if [[ ${#SNP_PARTS[@]} -ge 4 ]]; then

    SNP_REF="${SNP_PARTS[2]}"
    SNP_ALT="${SNP_PARTS[3]}"

    VARIANT_LINE="$(
      bcftools query \
        -r "$REGION" \
        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
        "$VCF" |
      awk -F'\t' \
        -v ref="$SNP_REF" \
        -v alt="$SNP_ALT" '
          $4 == ref && $5 == alt {
            print
            exit
          }
        '
    )"

  # ----------------------------------------------------------
  # Coordinate only
  # ----------------------------------------------------------

  else

    mapfile -t COORD_VARIANTS < <(
      bcftools query \
        -r "$REGION" \
        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
        "$VCF"
    )

    if [[ ${#COORD_VARIANTS[@]} -eq 0 ]]; then
      echo "ERROR: no Target ALS variant found at $REGION" >&2
      exit 1
    fi

    if [[ ${#COORD_VARIANTS[@]} -gt 1 ]]; then
      echo "ERROR: multiple variants exist at ${SNP_CHR}:${SNP_POS}." >&2
      echo "Specify REF and ALT explicitly:" >&2
      echo >&2

      for line in "${COORD_VARIANTS[@]}"; do
        IFS=$'\t' read -r -a tmp <<< "$line"

        echo "  ${tmp[0]}:${tmp[1]}:${tmp[3]}:${tmp[4]}" >&2
      done

      exit 1
    fi

    VARIANT_LINE="${COORD_VARIANTS[0]}"

  fi

fi

if [[ -z "$VARIANT_LINE" ]]; then
  echo "ERROR: requested variant could not be resolved in Target ALS VCF." >&2
  echo "Input: $SNP" >&2
  exit 1
fi

# ============================================================
# PARSE TARGET ALS VARIANT
# ============================================================

IFS=$'\t' read -r -a VAR_FIELDS <<< "$VARIANT_LINE"

if [[ ${#VAR_FIELDS[@]} -lt 6 ]]; then
  echo "ERROR: malformed Target ALS VCF query result." >&2
  exit 1
fi

VAR_CHR="${VAR_FIELDS[0]}"
VAR_POS="${VAR_FIELDS[1]}"
VAR_ID="${VAR_FIELDS[2]}"
VAR_REF="${VAR_FIELDS[3]}"
VAR_ALT="${VAR_FIELDS[4]}"

VAR_CHR="$(normalize_ucsc_chr "$VAR_CHR")"

EXPECTED_FIELDS=$((5 + ${#VCF_SAMPLES[@]}))

if [[ ${#VAR_FIELDS[@]} -ne "$EXPECTED_FIELDS" ]]; then
  echo "ERROR: genotype field count does not match VCF sample count." >&2
  echo "Expected: $EXPECTED_FIELDS" >&2
  echo "Observed: ${#VAR_FIELDS[@]}" >&2
  exit 1
fi

# ============================================================
# COORDINATE -> rsID
#
# Even if input was an rsID, validate/report it using the actual
# Target ALS REF/ALT that was ultimately selected.
# ============================================================

DB_MATCH="$(
  lookup_coordinate_dbsnp \
    "$VAR_CHR" \
    "$VAR_POS" \
    "$VAR_REF" \
    "$VAR_ALT" || true
)"

if [[ -n "$DB_MATCH" ]]; then

  IFS=$'\t' read -r \
    DB_CHR \
    DB_POS \
    DB_RSID \
    DB_REF \
    DB_ALT <<< "$DB_MATCH"

  RESOLVED_RSID="$DB_RSID"

elif [[ -n "$REQUESTED_RSID" ]]; then

  # We successfully mapped the supplied rsID to this variant,
  # but reverse lookup unexpectedly failed.
  RESOLVED_RSID="$REQUESTED_RSID"

  echo "WARNING: reverse dbSNP lookup did not reproduce the rsID." >&2

else

  RESOLVED_RSID="."

  echo "WARNING: no rsID found in local dbSNP for:" >&2
  echo "  ${VAR_CHR}:${VAR_POS}:${VAR_REF}:${VAR_ALT}" >&2

fi

# ============================================================
# REPORT FINAL RESOLUTION
# ============================================================

echo "============================================================"
echo "RESOLVED VARIANT"
echo "============================================================"
echo "Input:       $SNP"
echo "rsID:        $RESOLVED_RSID"
echo "Coordinate:  ${VAR_CHR}:${VAR_POS}"
echo "REF:         $VAR_REF"
echo "ALT:         $VAR_ALT"

if [[ -n "${DB_ALT:-}" ]]; then
  echo "dbSNP ALT:   $DB_ALT"
fi

echo

# ============================================================
# OUTPUT PREFIX
# ============================================================

if [[ "$RESOLVED_RSID" =~ ^rs[0-9]+$ ]]; then
  OUTPUT_SNP="$RESOLVED_RSID"
else
  OUTPUT_SNP="${VAR_CHR}:${VAR_POS}:${VAR_REF}:${VAR_ALT}"
fi

SAFE_SNP="$(sanitize "$OUTPUT_SNP")"
SAFE_ANCHOR="$(sanitize "$ANCHOR")"

PREFIX="${TISSUE}_${GENE}_${SAFE_SNP}_${SAFE_ANCHOR}"

VCF_GT_TSV="$OUTDIR/${PREFIX}.vcf_genotypes.tsv"

# ============================================================
# WRITE SUBJECT GENOTYPES
# ============================================================

{
  printf 'externalsubjectid\tGT\n'

  for ((i=0; i<${#VCF_SAMPLES[@]}; i++)); do

    printf '%s\t%s\n' \
      "${VCF_SAMPLES[$i]}" \
      "${VAR_FIELDS[$((5+i))]}"

  done

} > "$VCF_GT_TSV"

echo "Genotype table:"
echo "  $VCF_GT_TSV"
echo

echo "Genotype counts in full VCF:"
awk '
  NR > 1 {
    gt=$2

    gsub(/\|/, "/", gt)

    if (gt == "0/0") ref++
    else if (gt == "0/1" || gt == "1/0") het++
    else if (gt == "1/1") alt++
    else missing++
  }

  END {
    print "  Ref/Ref: " ref+0
    print "  Het:     " het+0
    print "  Hom Alt: " alt+0
    print "  Missing: " missing+0
  }
' "$VCF_GT_TSV"

echo

# ============================================================
# RUN STATIC LEAFVIZ-STYLE PLOT
# ============================================================

Rscript --vanilla "$R_SCRIPT" \
  --tissue "$TISSUE" \
  --tissue-regex "$TISSUE_REGEX" \
  --anchor "$ANCHOR" \
  --gene "$GENE" \
  --variant-chr "$VAR_CHR" \
  --variant-pos "$VAR_POS" \
  --variant-id "$RESOLVED_RSID" \
  --variant-ref "$VAR_REF" \
  --variant-alt "$VAR_ALT" \
  --genotypes "$VCF_GT_TSV" \
  --metadata "$META" \
  --data-root "$DATA" \
  --gtf "$GTF" \
  --outdir "$OUTDIR" \
  --prefix "$PREFIX"

echo
echo "============================================================"
echo "DONE"
echo "============================================================"
echo "Resolved SNP:"
echo "  rsID:       $RESOLVED_RSID"
echo "  Coordinate: ${VAR_CHR}:${VAR_POS}:${VAR_REF}:${VAR_ALT}"
echo
echo "Results:"
ls -lh "$OUTDIR"/"$PREFIX"* 2>/dev/null || true
