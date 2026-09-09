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
eval "$(conda shell.bash hook)"
conda activate RNAseq
module load R >/dev/null 2>&1
module load BCFtools >/dev/null 2>&1

usage() {
  cat <<'EOF'
Usage:
  sbatch run_sQTL_leafviz.sh TISSUE ANCHOR_JUNCTION SNP GENE

Arguments:
  TISSUE
    Motor_Cortex
    Cervical_Spinal_Cord
    Lumbar_Spinal_Cord
    Thoracic_Spinal_Cord
    Frontal_Cortex
    Cerebellum

  ANCHOR_JUNCTION
    Stable LeafCutter junction, e.g.
    chr19:-:17636157-17639083

  SNP
    One of:
      rs8106014
      chr19:17015573
      19:17015573
      chr19:17015573:C:G
      19:17015573:C:G

  GENE
    Gene symbol, e.g. UNC13A

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

[[ $# -eq 4 ]] || { usage; exit 1; }

TISSUE="$1"
ANCHOR="$2"
SNP="$3"
GENE="$4"

# ============================================================
# PATHS
# ============================================================
OUTDIR="$HOME/donglab/data/target_ALS/QTL/leafviz"
DATA="$HOME/donglab/data/target_ALS"
META="$DATA/targetALS_rnaseq_metadata.csv"
RAW="$DATA/QTL/plink/joint_all_chrs_matrixEQTL.raw"
SNP_MAT="$DATA/$TISSUE/sQTL/snp_${TISSUE}.txt"
SNP_LOC="$DATA/$TISSUE/sQTL/snp_location.txt"

ANNOT="$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode"
GTF="$ANNOT/gencode.v49.annotation.gtf"
DBSNP="$ANNOT/GCF_000001405.40.gz"
CHR_MAP="$ANNOT/ncbi_to_ucsc.txt"

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
R_SCRIPT="$SCRIPT_DIR/plot_sQTL_leafviz.R"

mkdir -p "$OUTDIR"

for f in "$META" "$RAW" "$SNP_MAT" "$SNP_LOC" "$GTF" "$DBSNP" "$DBSNP.tbi" "$CHR_MAP" "$R_SCRIPT"; do
  [[ -e "$f" ]] || { echo "ERROR: required file not found: $f" >&2; exit 1; }
done

for cmd in bcftools Rscript python awk sed grep head; do
  command -v "$cmd" >/dev/null 2>&1 || {
    echo "ERROR: $cmd is not available in PATH." >&2
    exit 1
  }
done

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
  printf '%s' "$1" | sed 's#[/: ]#_#g; s#[^A-Za-z0-9_.+-]#_#g'
}

normalize_ucsc_chr() {
  local chr="$1"
  [[ "$chr" == chr* ]] && printf '%s\n' "$chr" || printf 'chr%s\n' "$chr"
}

strip_chr() {
  printf '%s\n' "${1#chr}"
}

ucsc_to_ncbi_chr() {
  awk -v c="$1" '$2==c {print $1; exit}' "$CHR_MAP"
}

ncbi_to_ucsc_chr() {
  awk -v c="$1" '$1==c {print $2; exit}' "$CHR_MAP"
}

lookup_rsid_dbsnp() {
  local rsid="$1"
  bcftools query \
    -i "ID=\"$rsid\"" \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
    "$DBSNP" | head -1
}

lookup_coordinate_dbsnp() {
  local ucsc_chr="$1"
  local pos="$2"
  local ref="$3"
  local alt="$4"
  local ncbi_chr
  ncbi_chr="$(ucsc_to_ncbi_chr "$ucsc_chr")"

  [[ -n "$ncbi_chr" ]] || {
    echo "ERROR: could not convert $ucsc_chr using $CHR_MAP" >&2
    return 1
  }

  bcftools query \
    -r "${ncbi_chr}:${pos}-${pos}" \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
    "$DBSNP" |
  awk -F'\t' -v ref="$ref" -v alt="$alt" '
    $4==ref && $3 ~ /^rs/ {
      n=split($5,a,",")
      for(i=1;i<=n;i++) {
        if(a[i]==alt) {
          print
          exit
        }
      }
    }
  '
}

raw_columns_at_position() {
  local chr_no_chr="$1"
  local pos="$2"
  head -1 "$RAW" |
    tr '\t' '\n' |
    awk -F'[:_]' -v chr="$chr_no_chr" -v pos="$pos" '
      NF==5 {
        c=$1
        sub(/^chr/,"",c)
        if(c==chr && $2==pos) print
      }
    '
}

# ============================================================
# HEADER
# ============================================================
echo "============================================================"
echo "sQTL LeafViz"
echo "============================================================"
echo "Tissue:           $TISSUE"
echo "Anchor junction:  $ANCHOR"
echo "Requested SNP:    $SNP"
echo "Gene:             $GENE"
echo "Output directory: $OUTDIR"
echo

# ============================================================
# RESOLVE INPUT TO EXACT PLINK RAW VARIANT
#
# The RAW header is authoritative for:
#   chr / pos / REF / ALT / PLINK-counted allele
#
# Example:
#   19:17015573:C:G_C
#
# means:
#   REF=C
#   ALT=G
#   PLINK dosage counts C
# ============================================================
REQUESTED_RSID=""
RESOLVED_RSID=""
RAW_COL=""
VAR_CHR=""
VAR_CHR_NOCHR=""
VAR_POS=""
VAR_REF=""
VAR_ALT=""
COUNTED=""

if [[ "$SNP" =~ ^rs[0-9]+$ ]]; then
  REQUESTED_RSID="$SNP"
  echo "Input type: rsID"
  DB_LINE="$(lookup_rsid_dbsnp "$SNP")"

  [[ -n "$DB_LINE" ]] || {
    echo "ERROR: rsID not found in local dbSNP: $SNP" >&2
    exit 1
  }

  IFS=$'\t' read -r DB_CHR DB_POS DB_RSID DB_REF DB_ALT <<< "$DB_LINE"
  VAR_CHR="$(ncbi_to_ucsc_chr "$DB_CHR")"

  [[ -n "$VAR_CHR" ]] || {
    echo "ERROR: could not map dbSNP chromosome $DB_CHR to UCSC." >&2
    exit 1
  }

  VAR_CHR_NOCHR="$(strip_chr "$VAR_CHR")"
  VAR_POS="$DB_POS"

  echo "dbSNP:"
  echo "  rsID:       $DB_RSID"
  echo "  Coordinate: ${VAR_CHR}:${VAR_POS}"
  echo "  REF:        $DB_REF"
  echo "  ALT:        $DB_ALT"

  RAW_MATCHES="$(
    raw_columns_at_position "$VAR_CHR_NOCHR" "$VAR_POS" |
    awk -F'[:_]' -v ref="$DB_REF" -v dbalt="$DB_ALT" '
      BEGIN {
        n=split(dbalt,a,",")
        for(i=1;i<=n;i++) allowed[a[i]]=1
      }
      $3==ref && ($4 in allowed) {print}
    '
  )"
else
  echo "Input type: genomic coordinate"
  IFS=':' read -r -a SNP_PARTS <<< "$SNP"

  [[ ${#SNP_PARTS[@]} -ge 2 ]] || {
    echo "ERROR: coordinate must be chr:pos or chr:pos:REF:ALT." >&2
    exit 1
  }

  VAR_CHR="$(normalize_ucsc_chr "${SNP_PARTS[0]}")"
  VAR_CHR_NOCHR="$(strip_chr "$VAR_CHR")"
  VAR_POS="${SNP_PARTS[1]}"

  [[ "$VAR_POS" =~ ^[0-9]+$ ]] || {
    echo "ERROR: invalid position: $VAR_POS" >&2
    exit 1
  }

  RAW_MATCHES="$(raw_columns_at_position "$VAR_CHR_NOCHR" "$VAR_POS")"

  if [[ ${#SNP_PARTS[@]} -ge 4 ]]; then
    INPUT_REF="${SNP_PARTS[2]}"
    INPUT_ALT="${SNP_PARTS[3]}"
    RAW_MATCHES="$(
      printf '%s\n' "$RAW_MATCHES" |
      awk -F'[:_]' -v ref="$INPUT_REF" -v alt="$INPUT_ALT" '
        $3==ref && $4==alt {print}
      '
    )"
  fi
fi

N_RAW_MATCHES="$(printf '%s\n' "$RAW_MATCHES" | sed '/^$/d' | wc -l)"

if [[ "$N_RAW_MATCHES" -eq 0 ]]; then
  echo "ERROR: no matching variant found in PLINK RAW header." >&2
  echo "Input: $SNP" >&2
  echo "Resolved position: ${VAR_CHR}:${VAR_POS}" >&2
  exit 1
fi

if [[ "$N_RAW_MATCHES" -gt 1 ]]; then
  echo "ERROR: multiple matching variants found in PLINK RAW header:" >&2
  printf '%s\n' "$RAW_MATCHES" | sed 's/^/  /' >&2
  echo "Specify chr:pos:REF:ALT explicitly." >&2
  exit 1
fi

RAW_COL="$(printf '%s\n' "$RAW_MATCHES" | head -1)"
IFS=':_ ' read -r RAW_CHR VAR_POS VAR_REF VAR_ALT COUNTED <<< "$RAW_COL"
VAR_CHR="$(normalize_ucsc_chr "$RAW_CHR")"
VAR_CHR_NOCHR="$(strip_chr "$VAR_CHR")"

[[ "$COUNTED" == "$VAR_REF" || "$COUNTED" == "$VAR_ALT" ]] || {
  echo "ERROR: PLINK counted allele $COUNTED is neither REF=$VAR_REF nor ALT=$VAR_ALT." >&2
  exit 1
}

# ============================================================
# RESOLVE EXACT VARIANT BACK TO rsID
# ============================================================
DB_MATCH="$(lookup_coordinate_dbsnp "$VAR_CHR" "$VAR_POS" "$VAR_REF" "$VAR_ALT" || true)"

if [[ -n "$DB_MATCH" ]]; then
  IFS=$'\t' read -r DB_CHR DB_POS DB_RSID DB_REF DB_ALT <<< "$DB_MATCH"
  RESOLVED_RSID="$DB_RSID"
elif [[ -n "$REQUESTED_RSID" ]]; then
  RESOLVED_RSID="$REQUESTED_RSID"
  echo "WARNING: exact reverse dbSNP lookup failed; retaining requested rsID." >&2
else
  RESOLVED_RSID="."
  echo "WARNING: no rsID found for ${VAR_CHR}:${VAR_POS}:${VAR_REF}:${VAR_ALT}" >&2
fi

echo
echo "============================================================"
echo "RESOLVED VARIANT"
echo "============================================================"
echo "Input:          $SNP"
echo "rsID:           $RESOLVED_RSID"
echo "Coordinate:     ${VAR_CHR}:${VAR_POS}"
echo "REF:            $VAR_REF"
echo "ALT:            $VAR_ALT"
echo "PLINK counted:  $COUNTED"
echo "RAW column:     $RAW_COL"
[[ -n "${DB_ALT:-}" ]] && echo "dbSNP ALT set:  $DB_ALT"
echo

# ============================================================
# FIND TISSUE-SPECIFIC MATRIXEQTL SNP ID
# ============================================================
mapfile -t SNPID_MATCHES < <(
  awk -v chr="$VAR_CHR_NOCHR" -v pos="$VAR_POS" '
    BEGIN {FS="[ \t]+"}
    NR==1 {
      for(i=1;i<=NF;i++) {
        if($i=="snpid") sid=i
        if($i=="chr") ci=i
        if($i=="pos") pi=i
      }
      if(!sid || !ci || !pi) exit 2
      next
    }
    {
      c=$ci
      sub(/^chr/,"",c)
      if(c==chr && $pi==pos) print $sid
    }
  ' "$SNP_LOC"
)

if [[ ${#SNPID_MATCHES[@]} -eq 0 ]]; then
  echo "ERROR: no SNP found at ${VAR_CHR}:${VAR_POS} in $SNP_LOC" >&2
  exit 1
fi

EXACT_VARIANT="${VAR_CHR_NOCHR}:${VAR_POS}:${VAR_REF}:${VAR_ALT}"
SNPID=""

for candidate in "${SNPID_MATCHES[@]}"; do
  candidate_nochr="${candidate#chr}"
  if [[ "$candidate_nochr" == "$EXACT_VARIANT" ]]; then
    SNPID="$candidate"
    break
  fi
done

if [[ -z "$SNPID" && ${#SNPID_MATCHES[@]} -eq 1 ]]; then
  SNPID="${SNPID_MATCHES[0]}"
fi

if [[ -z "$SNPID" ]]; then
  echo "ERROR: multiple SNP IDs found at ${VAR_CHR}:${VAR_POS}, and exact allele match was ambiguous:" >&2
  printf '  %s\n' "${SNPID_MATCHES[@]}" >&2
  exit 1
fi

echo "MatrixEQTL SNP ID: $SNPID"
echo "Genotype matrix:   $SNP_MAT"
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
GENOTYPE_TSV="$OUTDIR/${PREFIX}.genotypes.tsv"

# ============================================================
# EXTRACT EXACT TISSUE-SPECIFIC GENOTYPES
#
# This reproduces _sQTL_boxplot.R:
#
#   snp_vals = tissue-specific MatrixEQTL genotype matrix
#
#   if PLINK counted allele == REF:
#       ALT dosage = 2 - PLINK dosage
#   else:
#       ALT dosage = PLINK dosage
#
# Then:
#   0 -> 0/0 -> Ref/Ref
#   1 -> 0/1 -> Het
#   2 -> 1/1 -> Hom Alt
#
# Matrix columns are mapped through targetALS metadata so the
# output always uses externalsubjectid for downstream joining.
# ============================================================
python - "$SNP_MAT" "$SNPID" "$COUNTED" "$VAR_REF" "$VAR_ALT" "$META" "$GENOTYPE_TSV" <<'PY'
import csv
import math
import sys

snp_file, snpid, counted, ref, alt, meta_file, out_file = sys.argv[1:]

subject_ids = set()
sample_to_subject = {}

with open(meta_file, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample = (row.get("externalsampleid") or "").strip()
        subject = (row.get("externalsubjectid") or "").strip()
        if not subject:
            continue
        subject_ids.add(subject)
        if sample:
            sample_to_subject[sample] = subject
            sample_to_subject[sample.replace("-", "_")] = subject

header = None
values = None

with open(snp_file) as f:
    for line in f:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split()
        if header is None:
            header = fields
            continue
        if fields[0] == snpid:
            values = fields
            break

if header is None:
    raise SystemExit(f"ERROR: empty SNP matrix: {snp_file}")
if values is None:
    raise SystemExit(f"ERROR: SNP {snpid} not found in genotype matrix: {snp_file}")
if len(values) != len(header):
    raise SystemExit(
        f"ERROR: SNP row/header length mismatch: "
        f"{len(values)} values vs {len(header)} columns"
    )

if counted not in (ref, alt):
    raise SystemExit(
        f"ERROR: counted allele {counted} is not REF={ref} or ALT={alt}"
    )

records = {}
unmapped = []

for matrix_id, raw_value in zip(header[1:], values[1:]):
    if raw_value.upper() in ("NA", "NAN", ".", ""):
        continue
    try:
        dosage = float(raw_value)
    except ValueError:
        continue

    alt_dosage = 2.0 - dosage if counted == ref else dosage

    nearest = round(alt_dosage)
    if not math.isfinite(alt_dosage) or abs(alt_dosage - nearest) > 1e-6:
        raise SystemExit(
            f"ERROR: non-integer genotype dosage for {matrix_id}: "
            f"raw={dosage}, ALT dosage={alt_dosage}"
        )

    alt_dosage = int(nearest)
    if alt_dosage == 0:
        gt = "0/0"
    elif alt_dosage == 1:
        gt = "0/1"
    elif alt_dosage == 2:
        gt = "1/1"
    else:
        raise SystemExit(
            f"ERROR: invalid ALT dosage {alt_dosage} for {matrix_id}"
        )

    if matrix_id in subject_ids:
        subject = matrix_id
    elif matrix_id in sample_to_subject:
        subject = sample_to_subject[matrix_id]
    elif matrix_id.replace("_", "-") in sample_to_subject:
        subject = sample_to_subject[matrix_id.replace("_", "-")]
    else:
        unmapped.append(matrix_id)
        continue

    if subject in records and records[subject] != gt:
        raise SystemExit(
            f"ERROR: conflicting genotypes for subject {subject}: "
            f"{records[subject]} vs {gt}"
        )

    records[subject] = gt

with open(out_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["externalsubjectid", "GT"])
    for subject in sorted(records):
        writer.writerow([subject, records[subject]])

counts = {"0/0": 0, "0/1": 0, "1/1": 0}
for gt in records.values():
    counts[gt] += 1

print(f"Mapped genotype subjects: {len(records)}")
print(f"Ref/Ref: {counts['0/0']}")
print(f"Het:     {counts['0/1']}")
print(f"Hom Alt: {counts['1/1']}")
print(f"Unmapped matrix IDs: {len(unmapped)}")

if unmapped:
    print("First unmapped IDs:", ", ".join(unmapped[:10]))
PY

echo
echo "Genotype table:"
echo "  $GENOTYPE_TSV"
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
  --genotypes "$GENOTYPE_TSV" \
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
echo "  rsID:          $RESOLVED_RSID"
echo "  Coordinate:    ${VAR_CHR}:${VAR_POS}:${VAR_REF}:${VAR_ALT}"
echo "  MatrixEQTL ID: $SNPID"
echo "  Counted allele:$COUNTED"
echo
echo "Results:"
ls -lh "$OUTDIR"/"$PREFIX"* 2>/dev/null || true
