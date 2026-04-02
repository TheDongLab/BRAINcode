#!/bin/bash
#SBATCH --job-name=run_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/run_eQTL_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/run_eQTL_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

#########################################
# run_eQTL.sh
# Purpose: Run Matrix eQTL (cis + trans) for a given tissue type.
#           Expects prep_eQTL.sh to have been run first.
#
# Covariate encoding applied before R:
#   - sex              -> binary (Male=1, Female=0)
#   - ethnicity        -> one-hot dummy columns (drop first level)
#   - subject_group    -> one-hot dummy columns (drop first level)
#   - subject_group_subcategory -> one-hot dummy columns (drop first level)
#   - tissue           -> DROPPED (constant within tissue subset)
#   - age_at_death, age_at_symptom_onset, rin -> numeric as-is
#   - externalsampleid, externalsubjectid -> DROPPED (ID columns)
#   - all other non-numeric free-text columns -> DROPPED with warning
#
# Usage:
#   sbatch run_eQTL.sh "Cerebellum"
#   sbatch run_eQTL.sh "Frontal Cortex"
#   bash   run_eQTL.sh "Cerebellum"
#########################################

set -euo pipefail

# ---- Tissue argument ------------------------------------------------
if [ $# -lt 1 ]; then
    echo "ERROR: No tissue specified."
    echo "Usage: sbatch run_eQTL.sh \"Cerebellum\""
    exit 1
fi

TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

echo "============================================"
echo "  Matrix eQTL run for tissue : $TISSUE"
echo "  SLURM job ID               : ${SLURM_JOB_ID:-local}"
echo "  $(date)"
echo "============================================"

# ---- Paths ----------------------------------------------------------
BASE=/home/zw529/donglab/data/target_ALS
PIPELINE=/home/zw529/donglab/pipelines/scripts/eQTL

INDIR=$BASE/$TISSUE_DIR/eQTL
OUTDIR=$INDIR/results
mkdir -p $OUTDIR

SNP_FILE=$INDIR/snp_${TISSUE_DIR}.txt
EXPR_FILE=$INDIR/expression_${TISSUE_DIR}.txt
COV_FILE=$INDIR/covariates_${TISSUE_DIR}.txt
COV_ENCODED=$INDIR/covariates_${TISSUE_DIR}_encoded.txt   # numeric-only, written here
GENE_LOC=$INDIR/gene_location.txt
SNP_LOC=$INDIR/snp_location.txt
OUTPUT_PREFIX=$OUTDIR/${TISSUE_DIR}_eQTL

# ---- Pre-flight: check input files exist ----------------------------
echo ""
echo "[0] Checking input files exist..."

MISSING=0
for f in "$SNP_FILE" "$EXPR_FILE" "$COV_FILE" "$GENE_LOC" "$SNP_LOC"; do
    if [ ! -f "$f" ]; then
        echo "  MISSING: $f"
        MISSING=1
    else
        SIZE=$(du -sh "$f" | cut -f1)
        echo "  OK ($SIZE): $f"
    fi
done

if [ "$MISSING" -eq 1 ]; then
    echo ""
    echo "ERROR: One or more input files are missing."
    echo "Run: sbatch prep_eQTL.sh \"$TISSUE\" first."
    exit 1
fi

# ---- Pre-flight: verify sample counts match -------------------------
echo ""
echo "[0b] Verifying sample count alignment..."

python3 << EOF
import sys

def col_count(fpath):
    with open(fpath) as f:
        return len(f.readline().rstrip('\n').split('\t')) - 1  # minus row-label col

n_expr = col_count("$EXPR_FILE")
n_snp  = col_count("$SNP_FILE")
n_cov  = col_count("$COV_FILE")

print(f"  Expression samples : {n_expr}", flush=True)
print(f"  SNP samples        : {n_snp}", flush=True)
print(f"  Covariate samples  : {n_cov}", flush=True)

ok = True
if n_expr != n_snp:
    print(f"  ERROR: Expression ({n_expr}) and SNP ({n_snp}) sample counts differ.", flush=True)
    ok = False
if n_cov != n_expr:
    print(f"  ERROR: Covariate ({n_cov}) and expression ({n_expr}) sample counts differ.", flush=True)
    ok = False
if ok:
    print(f"  All sample counts match: {n_expr}", flush=True)
else:
    sys.exit(1)
EOF

# ---- Encode covariates to numeric -----------------------------------
echo ""
echo "[1] Encoding covariates to numeric for Matrix eQTL..."

python3 << EOF
import csv

cov_file = "$COV_FILE"
out_file = "$COV_ENCODED"

# Columns to drop entirely
DROP_COLS = {'externalsampleid', 'externalsubjectid', 'tissue',
             'subject_group_subcategory'}   # too many levels, sparse dummies

# Binary encoding map: col -> {string_value -> numeric}
BINARY_COLS = {
    'sex': {'male': 1, 'female': 0, 'm': 1, 'f': 0,
            'Male': 1, 'Female': 0, 'M': 1, 'F': 0}
}

# Columns to treat as direct numeric (pass through, coerce to float)
NUMERIC_COLS = {'age_at_death', 'age_at_symptom_onset', 'rin'}

# Columns to one-hot encode (drop first level to avoid multicollinearity)
ONEHOT_COLS = {'ethnicity', 'subject_group'}

# Read wide-format covariate file: rows = covariates, cols = samples
# Header row: covariate  sample1  sample2 ...
with open(cov_file) as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)   # ['covariate', sample1, sample2, ...]
    sample_ids = header[1:]

    rows = {}   # covariate_name -> [values per sample]
    for row in reader:
        rows[row[0]] = row[1:]

out_rows = []   # list of (row_name, [numeric values])
dropped  = []
warnings = []

for cov_name, values in rows.items():
    cov_lower = cov_name.lower()

    # Drop columns
    if cov_name in DROP_COLS or cov_lower in DROP_COLS:
        dropped.append(cov_name)
        continue

    # Binary encoding
    if cov_name in BINARY_COLS or cov_lower in {k.lower() for k in BINARY_COLS}:
        mapping = BINARY_COLS.get(cov_name) or BINARY_COLS.get(cov_lower, {})
        encoded = []
        for v in values:
            v_stripped = v.strip()
            mapped = mapping.get(v_stripped, mapping.get(v_stripped.lower(), None))
            if mapped is not None:
                encoded.append(str(mapped))
            elif v_stripped in ('NA', '', 'nan', 'NaN'):
                encoded.append('NA')
            else:
                encoded.append('NA')
                warnings.append(f"  WARNING: Unknown value '{v_stripped}' in binary col '{cov_name}' -> NA")
        out_rows.append((cov_name, encoded))
        continue

    # Direct numeric
    if cov_name in NUMERIC_COLS or cov_lower in NUMERIC_COLS:
        encoded = []
        for v in values:
            v_stripped = v.strip()
            try:
                encoded.append(str(float(v_stripped)))
            except ValueError:
                encoded.append('NA')
        out_rows.append((cov_name, encoded))
        continue

    # One-hot encoding
    if cov_name in ONEHOT_COLS or cov_lower in ONEHOT_COLS:
        unique_levels = sorted(set(
            v.strip() for v in values
            if v.strip() not in ('NA', '', 'nan', 'NaN')
        ))
        reference = unique_levels[0]   # drop first level (reference)
        for level in unique_levels[1:]:
            dummy_name = f"{cov_name}_{level.replace(' ', '_').replace('/', '_')}"
            encoded = []
            for v in values:
                v_stripped = v.strip()
                if v_stripped in ('NA', '', 'nan', 'NaN'):
                    encoded.append('NA')
                else:
                    encoded.append('1' if v_stripped == level else '0')
            out_rows.append((dummy_name, encoded))
        print(f"  One-hot: '{cov_name}' -> {len(unique_levels)-1} dummy cols (ref='{reference}')", flush=True)
        continue

    # Unknown column: try numeric coercion, drop if all NA or non-numeric
    numeric_vals = []
    n_valid = 0
    for v in values:
        v_stripped = v.strip()
        try:
            numeric_vals.append(str(float(v_stripped)))
            n_valid += 1
        except ValueError:
            numeric_vals.append('NA')
    if n_valid > 0:
        out_rows.append((cov_name, numeric_vals))
        print(f"  Numeric coercion: '{cov_name}' ({n_valid}/{len(values)} valid)", flush=True)
    else:
        dropped.append(cov_name)
        warnings.append(f"  DROPPED non-numeric col: '{cov_name}'")

# Write encoded file: rows = covariates, cols = samples
with open(out_file, 'w') as f:
    f.write('covariate\t' + '\t'.join(sample_ids) + '\n')
    for row_name, vals in out_rows:
        f.write(row_name + '\t' + '\t'.join(vals) + '\n')

print(f"", flush=True)
print(f"  Dropped columns    : {dropped}", flush=True)
for w in warnings:
    print(w, flush=True)
print(f"  Output covariates  : {len(out_rows)} rows x {len(sample_ids)} samples", flush=True)
print(f"  Written to         : {out_file}", flush=True)
EOF

echo "  Covariate encoding done."

# ---- Run Matrix eQTL ------------------------------------------------
echo ""
echo "[2] Running Matrix eQTL (cis + trans)..."
echo "  SNP file    : $SNP_FILE"
echo "  Expr file   : $EXPR_FILE"
echo "  Cov file    : $COV_ENCODED"
echo "  Gene loc    : $GENE_LOC"
echo "  SNP loc     : $SNP_LOC"
echo "  Output      : $OUTPUT_PREFIX"
echo ""

Rscript $PIPELINE/_eQTL.R \
    "$SNP_FILE" \
    "$EXPR_FILE" \
    "$COV_ENCODED" \
    "$OUTPUT_PREFIX" \
    "$GENE_LOC" \
    "$SNP_LOC"

# ---- Post-run summary -----------------------------------------------
echo ""
echo "[3] Results summary for $TISSUE..."

CIS_FILE="${OUTPUT_PREFIX}.cis.txt"
TRANS_FILE="${OUTPUT_PREFIX}.trans.txt"
PDF_FILE="${OUTPUT_PREFIX}.pdf"

if [ -f "$CIS_FILE" ]; then
    N_CIS=$(tail -n +2 "$CIS_FILE" | wc -l)
    echo "  Cis eQTLs detected  : $N_CIS"
    echo "  Top 5 cis eQTLs:"
    head -6 "$CIS_FILE" | column -t
else
    echo "  WARNING: Cis output file not found: $CIS_FILE"
fi

if [ -f "$TRANS_FILE" ]; then
    N_TRANS=$(tail -n +2 "$TRANS_FILE" | wc -l)
    echo "  Trans eQTLs detected: $N_TRANS"
    echo "  Top 5 trans eQTLs:"
    head -6 "$TRANS_FILE" | column -t
else
    echo "  WARNING: Trans output file not found: $TRANS_FILE"
fi

[ -f "$PDF_FILE" ] && echo "  QQ plot saved to    : $PDF_FILE"

echo ""
echo "============================================"
echo "  Matrix eQTL complete for : $TISSUE"
echo "  Results in               : $OUTDIR"
echo "  $(date)"
echo "============================================"
