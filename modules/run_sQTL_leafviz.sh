#!/bin/bash
#SBATCH --job-name=sQTL_leafviz
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/leafviz/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/leafviz/%x_%j.err

set -euo pipefail

module load R
module load BCFtools

conda activate RNAseq

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
                   chr19:123456:C:G
                   19:123456:C:G

  GENE             Gene symbol, e.g. UNC13A

Example:
  sbatch run_sQTL_leafviz.sh \
    Frontal_Cortex \
    chr19:-:17641556-17642845 \
    19:123456:C:G \
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

OUTDIR="$HOME/donglab/data/target_ALS/QTL/leafviz"
mkdir -p "$OUTDIR"

VCF="$HOME/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
META="$HOME/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
DATA="$HOME/donglab/data/target_ALS"
GTF="$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="$SCRIPT_DIR/plot_sQTL_leafviz.R"

for f in "$VCF" "$META" "$GTF" "$R_SCRIPT"; do
  if [[ ! -e "$f" ]]; then
    echo "ERROR: required file not found: $f" >&2
    exit 1
  fi
done

# Use the same environment as the RNA-seq/LeafCutter workflow.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate RNAseq

command -v bcftools >/dev/null 2>&1 || {
  echo "ERROR: bcftools is not available in PATH." >&2
  exit 1
}

# Tissue regex is passed to R and applied to metadata$tissue.
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

sanitize() {
  printf '%s' "$1" | sed 's#[/: ]#_#g; s#[^A-Za-z0-9_.+-]#_#g'
}

SAFE_SNP="$(sanitize "$SNP")"
SAFE_ANCHOR="$(sanitize "$ANCHOR")"
PREFIX="${TISSUE}_${GENE}_${SAFE_SNP}_${SAFE_ANCHOR}"

VCF_GT_TSV="$OUTDIR/${PREFIX}.vcf_genotypes.tsv"

echo "============================================================"
echo "sQTL LeafViz"
echo "============================================================"
echo "Tissue:  $TISSUE"
echo "Anchor:  $ANCHOR"
echo "SNP:     $SNP"
echo "Gene:    $GENE"
echo "Output:  $OUTDIR"
echo

# Resolve the requested SNP and collect GT for every VCF subject.
mapfile -t VCF_SAMPLES < <(bcftools query -l "$VCF")

if [[ ${#VCF_SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: no samples found in VCF." >&2
  exit 1
fi

VARIANT_LINE=""

if [[ "$SNP" == rs* ]]; then
  VARIANT_LINE="$(bcftools query \
    -i "ID=\"$SNP\"" \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
    "$VCF" | head -1 || true)"
else
  IFS=':' read -r -a SNP_PARTS <<< "$SNP"

  if [[ ${#SNP_PARTS[@]} -lt 2 ]]; then
    echo "ERROR: SNP must be rsID, chr:pos, or chr:pos:REF:ALT." >&2
    exit 1
  fi

  SNP_CHR="${SNP_PARTS[0]}"
  SNP_POS="${SNP_PARTS[1]}"

  [[ "$SNP_CHR" == chr* ]] || SNP_CHR="chr${SNP_CHR}"
  REGION="${SNP_CHR}:${SNP_POS}-${SNP_POS}"

  if [[ ${#SNP_PARTS[@]} -ge 4 ]]; then
    SNP_REF="${SNP_PARTS[2]}"
    SNP_ALT="${SNP_PARTS[3]}"
    VARIANT_LINE="$(bcftools query \
      -r "$REGION" \
      -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
      "$VCF" |
      awk -F'\t' -v ref="$SNP_REF" -v alt="$SNP_ALT" \
        '$4==ref && $5==alt {print; exit}' || true)"
  else
    VARIANT_LINE="$(bcftools query \
      -r "$REGION" \
      -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
      "$VCF" | head -1 || true)"
  fi
fi

if [[ -z "$VARIANT_LINE" ]]; then
  echo "ERROR: SNP not found in VCF: $SNP" >&2
  exit 1
fi

IFS=$'\t' read -r -a VAR_FIELDS <<< "$VARIANT_LINE"

if [[ ${#VAR_FIELDS[@]} -lt 6 ]]; then
  echo "ERROR: malformed bcftools query result for SNP: $SNP" >&2
  exit 1
fi

VAR_CHR="${VAR_FIELDS[0]}"
VAR_POS="${VAR_FIELDS[1]}"
VAR_ID="${VAR_FIELDS[2]}"
VAR_REF="${VAR_FIELDS[3]}"
VAR_ALT="${VAR_FIELDS[4]}"

EXPECTED_FIELDS=$((5 + ${#VCF_SAMPLES[@]}))
if [[ ${#VAR_FIELDS[@]} -ne $EXPECTED_FIELDS ]]; then
  echo "ERROR: genotype field count does not match VCF sample count." >&2
  echo "Expected $EXPECTED_FIELDS fields, got ${#VAR_FIELDS[@]}." >&2
  exit 1
fi

{
  printf 'externalsubjectid\tGT\n'
  for ((i=0; i<${#VCF_SAMPLES[@]}; i++)); do
    printf '%s\t%s\n' "${VCF_SAMPLES[$i]}" "${VAR_FIELDS[$((5+i))]}"
  done
} > "$VCF_GT_TSV"

echo "Resolved variant:"
echo "  ${VAR_CHR}:${VAR_POS}:${VAR_REF}:${VAR_ALT}"
echo "  ID: ${VAR_ID}"
echo "Genotype table:"
echo "  $VCF_GT_TSV"
echo

# The plotting script itself does not require the leafcutter R package.
Rscript "$R_SCRIPT" \
  --tissue "$TISSUE" \
  --tissue-regex "$TISSUE_REGEX" \
  --anchor "$ANCHOR" \
  --gene "$GENE" \
  --variant-chr "$VAR_CHR" \
  --variant-pos "$VAR_POS" \
  --variant-id "$VAR_ID" \
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
echo "Results:"
ls -lh "$OUTDIR"/"$PREFIX"* 2>/dev/null || true
