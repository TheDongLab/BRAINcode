#!/bin/bash
#SBATCH --job-name=prep_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=499G

###########################################
# prep_eQTL.sh
# Purpose: Prepare all Matrix eQTL inputs for a given tissue type.
#
# Steps:
#   1.  Build sample mapping (HRA <-> subjectid <-> HDA)
#   2.  Subset expression matrix to tissue samples
#   3.  Subset + transpose SNP matrix to tissue samples -> SNP matrix is the AUTHORITATIVE sample set after this step
#   4.  Subset covariates to match SNP-confirmed samples
#   5.  Align expression + covariates to exact SNP sample order/count
#   6.  Encode covariates to numeric (Matrix eQTL requirement)
#   7.  Generate SNP location file -> SNP names taken from .raw header (authority) matched positionally to chr/pos from .bim
#   8.  Generate gene location file from gencode v49 BED6 -> Strips ___type___name decorator, keeps versioned ENSG ID
#
# Logic: Keeps subjects with BOTH:
#   - RNAseq for the specified tissue in expression_sample_metadata.csv
#   - WGS data (any tissue) in wgs_samples_for_vcf_merge.csv
#   Subjects in the WGS manifest but absent from the .raw genotype
#   file are automatically excluded during alignment (Step 5).
#
# Covariate encoding (Step 6):
#   - sex                       -> binary (Male=1, Female=0)
#   - ethnicity                 -> one-hot dummies (drop first level)
#                                  normalised to lowercase before encoding
#                                  to avoid duplicate levels from typos
#   - subject_group             -> one-hot dummies (drop first level)
#                                  normalised: lowercase, hyphens->underscores
#   - subject_group_subcategory -> DROPPED (too many sparse levels)
#   - tissue                    -> DROPPED (constant within tissue)
#   - age_at_death, age_at_symptom_onset, rin -> numeric as-is
#   - externalsampleid, externalsubjectid     -> DROPPED (ID cols)
#   - other non-numeric columns -> DROPPED with warning
#
# Usage:
#   sbatch prep_eQTL.sh "Cerebellum"
#   sbatch prep_eQTL.sh "Frontal Cortex"
#   bash   prep_eQTL.sh "Cerebellum"
#
# Output: $BASE/[TISSUE_DIR]/eQTL/
###########################################

set -euo pipefail

# ── Tissue argument ───────────────────────────────────────────────────
if [ $# -lt 1 ]; then
    echo "ERROR: No tissue specified."
    echo "Usage: sbatch prep_eQTL.sh \"Cerebellum\""
    echo "       sbatch prep_eQTL.sh \"Frontal Cortex\""
    exit 1
fi

TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

echo "============================================"
echo "  eQTL Prep for tissue : $TISSUE"
echo "  Directory label      : $TISSUE_DIR"
echo "  $(date)"
echo "============================================"

# ── Paths ─────────────────────────────────────────────────────────────
BASE=/home/zw529/donglab/data/target_ALS
PLINK=$BASE/eQTL/plink
REFS=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode

EXPR=$BASE/eQTL/expression_matrix.txt
RAW=$PLINK/joint_autosomes_matrixEQTL.raw
META=$BASE/eQTL/expression_sample_metadata.csv
WGS_MAP=$BASE/eQTL/wgs_samples_for_vcf_merge.csv
COV=$BASE/eQTL/covariates.tsv
BIM=$PLINK/joint_autosomes_filtered_bed.bim
GTF_BED6=$REFS/gencode.v49.annotation.gene.bed6

OUTDIR=$BASE/$TISSUE_DIR/eQTL
mkdir -p $OUTDIR

# Temp files in OUTDIR so they persist across compute nodes
TMP_META=$OUTDIR/tmp_meta.txt
TMP_WGS=$OUTDIR/tmp_wgs_map.txt
TMP_MATCHED=$OUTDIR/tmp_matched.txt
TMP_HRA=$OUTDIR/tmp_HRA_ids.txt
TMP_HRA_UNDER=$OUTDIR/tmp_HRA_ids_underscore.txt
TMP_HDA=$OUTDIR/tmp_HDA_ids.txt

##############################################
# STEP 1: Build sample mapping table
##############################################
echo ""
echo "[1] Building sample mapping (HRA <-> subjectid <-> HDA) for: $TISSUE"

tail -n +2 $META \
    | awk -F',' -v tissue="$TISSUE" '$3 == tissue {print $2","$1}' \
    | sort -t',' -k1,1 \
    > $TMP_META

N_META=$(wc -l < $TMP_META)
echo "  RNA samples in metadata for $TISSUE : $N_META"

if [ "$N_META" -eq 0 ]; then
    echo "ERROR: No samples found for tissue '$TISSUE'."
    echo "Available tissue labels:"
    tail -n +2 $META | awk -F',' '{print $3}' | sort | uniq -c | sort -rn
    exit 1
fi

tail -n +2 $WGS_MAP \
    | awk -F',' '{print $1","$2}' \
    | sort -t',' -k1,1 \
    > $TMP_WGS

N_WGS=$(wc -l < $TMP_WGS)
echo "  Subjects with WGS data             : $N_WGS"

join -t',' -1 1 -2 1 $TMP_META $TMP_WGS > $TMP_MATCHED

N_MATCHED=$(wc -l < $TMP_MATCHED)
echo "  Subjects with $TISSUE RNA + WGS    : $N_MATCHED"

if [ "$N_MATCHED" -eq 0 ]; then
    echo "ERROR: No matched samples found. Check subjectid formatting."
    exit 1
fi

awk -F',' '{print $2}' $TMP_MATCHED > $TMP_HRA
awk -F',' '{print $3}' $TMP_MATCHED > $TMP_HDA
sed 's/-/_/g' $TMP_HRA > $TMP_HRA_UNDER

echo "  Sample HRA IDs (expression): $(head -3 $TMP_HRA | tr '\n' ' ')..."
echo "  Sample HDA IDs (SNP):        $(head -3 $TMP_HDA | tr '\n' ' ')..."

##############################################
# STEP 2: Subset expression matrix
##############################################
echo ""
echo "[2] Subsetting expression matrix..."

python3 << EOF
ids_file  = "$TMP_HRA_UNDER"
expr_file = "$EXPR"
out_file  = "$OUTDIR/expression_${TISSUE_DIR}.txt"

with open(ids_file) as f:
    keep_ids = set(line.strip() for line in f if line.strip())

with open(expr_file) as f, open(out_file, 'w') as out:
    header = f.readline().rstrip('\n').split('\t')
    keep_idx = [0] + [i for i, h in enumerate(header) if h in keep_ids]
    found = [header[i] for i in keep_idx[1:]]
    missing = keep_ids - set(found)

    print(f"  Requested : {len(keep_ids)}", flush=True)
    print(f"  Found     : {len(found)}", flush=True)
    if missing:
        print(f"  WARNING: {len(missing)} IDs not in expression matrix: {list(missing)[:5]}", flush=True)

    out.write('\t'.join(header[i] for i in keep_idx) + '\n')
    for line in f:
        fields = line.rstrip('\n').split('\t')
        out.write('\t'.join(fields[i] for i in keep_idx) + '\n')

print(f"  Written to: {out_file}", flush=True)
EOF

echo "  Expression matrix done."

##############################################
# STEP 3: Subset + transpose SNP matrix
# After this step the SNP file is the AUTHORITATIVE sample set.
# SNP row names are taken directly from the .raw header — these are
# in PLINK's SNPID_ALLELE format (e.g. ._AC) and must be used
# consistently in the SNP location file (Step 7).
##############################################
echo ""
echo "[3] Subsetting and transposing SNP matrix (chunked)..."

python3 << EOF
import re

ids_file = "$TMP_HDA"
raw_file = "$RAW"
out_file = "$OUTDIR/snp_${TISSUE_DIR}.txt"

with open(ids_file) as f:
    keep_ids = set(line.strip() for line in f if line.strip())

print(f"  HDA IDs requested: {len(keep_ids)}", flush=True)

def normalize_hda(iid):
    return re.sub(r'-b\d+$', '', iid)

rows = []
with open(raw_file) as f:
    header = f.readline().split()
    # SNP names as PLINK wrote them (SNPID_ALLELE format) — used as-is
    snp_names = header[6:]
    for line in f:
        fields = line.split()
        iid_raw  = fields[1]
        iid_norm = normalize_hda(iid_raw)
        matched  = iid_norm if iid_norm in keep_ids else (iid_raw if iid_raw in keep_ids else None)
        if matched:
            rows.append((matched, fields[6:]))

print(f"  HDA IDs found in .raw : {len(rows)}", flush=True)
missing = keep_ids - set(r[0] for r in rows)
if missing:
    print(f"  NOTE: {len(missing)} HDA IDs absent from .raw (excluded from genotype data): {list(missing)[:5]}", flush=True)

sample_ids = [r[0] for r in rows]
n_samples  = len(rows)
n_snps     = len(snp_names)
CHUNK      = 50000

with open(out_file, 'w') as out:
    out.write('snpid\t' + '\t'.join(sample_ids) + '\n')
    for chunk_start in range(0, n_snps, CHUNK):
        chunk_end = min(chunk_start + CHUNK, n_snps)
        for j in range(chunk_start, chunk_end):
            out.write(snp_names[j] + '\t' + '\t'.join(rows[i][1][j] for i in range(n_samples)) + '\n')
        print(f"  Written SNPs {chunk_start}-{chunk_end} of {n_snps}", flush=True)

print(f"  Done: {n_snps} SNPs x {n_samples} samples -> {out_file}", flush=True)
EOF

echo "  SNP matrix done."

##############################################
# STEP 4: Subset covariates (long -> wide)
##############################################
echo ""
echo "[4] Subsetting and transposing covariates..."

python3 << EOF
ids_file = "$TMP_HRA"
cov_file = "$COV"
out_file = "$OUTDIR/covariates_${TISSUE_DIR}.txt"

with open(ids_file) as f:
    keep_ids = set(line.strip() for line in f if line.strip())

header = None
kept_rows = {}

with open(cov_file) as f:
    header = f.readline().rstrip('\n').split('\t')
    cov_cols = header[1:]
    for line in f:
        fields = line.rstrip('\n').split('\t')
        if fields[0] in keep_ids:
            kept_rows[fields[0]] = fields[1:]

found   = len(kept_rows)
missing = keep_ids - set(kept_rows.keys())
print(f"  Requested : {len(keep_ids)}", flush=True)
print(f"  Found     : {found}", flush=True)
if missing:
    print(f"  WARNING: {len(missing)} IDs not in covariates: {list(missing)[:5]}", flush=True)

sample_order = sorted(kept_rows.keys())

with open(out_file, 'w') as out:
    out.write('covariate\t' + '\t'.join(sample_order) + '\n')
    for i, cov_name in enumerate(cov_cols):
        vals = [kept_rows[s][i] if s in kept_rows else 'NA' for s in sample_order]
        out.write(cov_name + '\t' + '\t'.join(vals) + '\n')

print(f"  Transposed {len(cov_cols)} covariates x {len(sample_order)} samples", flush=True)
print(f"  Written to: {out_file}", flush=True)
EOF

echo "  Covariates done."

##############################################
# STEP 5: Align expression + covariates to SNP sample set
##############################################
echo ""
echo "[5] Aligning expression and covariates to SNP sample set..."

python3 << EOF
snp_file  = "$OUTDIR/snp_${TISSUE_DIR}.txt"
expr_file = "$OUTDIR/expression_${TISSUE_DIR}.txt"
cov_file  = "$OUTDIR/covariates_${TISSUE_DIR}.txt"
meta_file = "$META"
wgs_file  = "$WGS_MAP"
tissue    = "$TISSUE"

with open(snp_file) as f:
    hda_ordered = f.readline().rstrip('\n').split('\t')[1:]
hda_set = set(hda_ordered)
print(f"  SNP matrix samples (authoritative): {len(hda_ordered)}", flush=True)

subj_to_hda = {}
with open(wgs_file) as f:
    f.readline()
    for line in f:
        parts = line.strip().split(',')
        if len(parts) >= 2:
            subj_to_hda[parts[0]] = parts[1]

subj_to_hra = {}
with open(meta_file) as f:
    f.readline()
    for line in f:
        parts = line.strip().split(',')
        if len(parts) >= 3 and parts[2] == tissue:
            subj_to_hra[parts[1]] = parts[0]

hda_to_hra_under  = {}
hda_to_hra_hyphen = {}
for subj, hda in subj_to_hda.items():
    if hda in hda_set and subj in subj_to_hra:
        hra_hyphen = subj_to_hra[subj]
        hda_to_hra_under[hda]  = hra_hyphen.replace('-', '_')
        hda_to_hra_hyphen[hda] = hra_hyphen

unmapped = [h for h in hda_ordered if h not in hda_to_hra_under]
if unmapped:
    print(f"  WARNING: {len(unmapped)} HDA IDs have no HRA mapping: {unmapped[:5]}", flush=True)

hra_under_ordered  = [hda_to_hra_under.get(h)  for h in hda_ordered]
hra_hyphen_ordered = [hda_to_hra_hyphen.get(h) for h in hda_ordered]

print(f"  Re-ordering expression matrix...", flush=True)
with open(expr_file) as f:
    header = f.readline().rstrip('\n').split('\t')
    idx_map  = {h: i for i, h in enumerate(header)}
    keep_idx = [0] + [idx_map[h] for h in hra_under_ordered if h and h in idx_map]
    rows = ['\t'.join(header[i] for i in keep_idx)]
    for line in f:
        fields = line.rstrip('\n').split('\t')
        rows.append('\t'.join(fields[i] for i in keep_idx))
with open(expr_file, 'w') as f:
    f.write('\n'.join(rows) + '\n')
print(f"  Expression aligned: {len(keep_idx)-1} samples", flush=True)

print(f"  Re-ordering covariate matrix...", flush=True)
with open(cov_file) as f:
    header = f.readline().rstrip('\n').split('\t')
    idx_map  = {h: i for i, h in enumerate(header)}
    keep_idx = [0] + [idx_map[h] for h in hra_hyphen_ordered if h and h in idx_map]
    rows = ['\t'.join(header[i] for i in keep_idx)]
    for line in f:
        fields = line.rstrip('\n').split('\t')
        rows.append('\t'.join(fields[i] for i in keep_idx))
with open(cov_file, 'w') as f:
    f.write('\n'.join(rows) + '\n')
print(f"  Covariates aligned: {len(keep_idx)-1} samples", flush=True)

def col_count(fpath):
    with open(fpath) as f:
        return len(f.readline().rstrip('\n').split('\t')) - 1

n_snp  = col_count(snp_file)
n_expr = col_count(expr_file)
n_cov  = col_count(cov_file)
print(f"  Final sample counts -> SNP: {n_snp}  Expr: {n_expr}  Cov: {n_cov}", flush=True)

if n_snp == n_expr == n_cov:
    print(f"  All counts match. Alignment complete.", flush=True)
else:
    print(f"  ERROR: Counts still differ after alignment.", flush=True)
    exit(1)
EOF

echo "  Alignment done."

##############################################
# STEP 6: Encode covariates to numeric
# Categorical variables are normalised before encoding to prevent
# duplicate dummy columns from capitalisation/punctuation inconsistencies:
#   - ethnicity: lowercased before level detection
#   - subject_group: lowercased + hyphens replaced with underscores
##############################################
echo ""
echo "[6] Encoding covariates to numeric..."

python3 << EOF
import csv

cov_file = "$OUTDIR/covariates_${TISSUE_DIR}.txt"
out_file = "$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt"

DROP_COLS    = {'externalsampleid', 'externalsubjectid', 'tissue',
                'subject_group_subcategory'}
BINARY_COLS  = {'sex': {'male': 1, 'female': 0, 'm': 1, 'f': 0,
                        'Male': 1, 'Female': 0, 'M': 1, 'F': 0}}
NUMERIC_COLS = {'age_at_death', 'age_at_symptom_onset', 'rin'}
ONEHOT_COLS  = {'ethnicity', 'subject_group'}

def normalise_level(col, val):
    """Normalise categorical values to prevent duplicate dummy columns."""
    v = val.strip()
    if col == 'ethnicity':
        return v.lower()
    if col == 'subject_group':
        return v.lower().replace('-', '_').replace(',', '').replace('/', '_')
    return v

with open(cov_file) as f:
    reader     = csv.reader(f, delimiter='\t')
    header     = next(reader)
    sample_ids = header[1:]
    rows       = {row[0]: row[1:] for row in reader}

out_rows = []
dropped  = []
warnings = []

for cov_name, values in rows.items():
    cov_lower = cov_name.lower()

    if cov_name in DROP_COLS or cov_lower in DROP_COLS:
        dropped.append(cov_name)
        continue

    if cov_name in BINARY_COLS or cov_lower in BINARY_COLS:
        mapping = BINARY_COLS.get(cov_name) or BINARY_COLS.get(cov_lower, {})
        encoded = []
        for v in values:
            v = v.strip()
            mapped = mapping.get(v, mapping.get(v.lower()))
            encoded.append(str(mapped) if mapped is not None else 'NA')
            if mapped is None and v not in ('NA', '', 'nan', 'NaN'):
                warnings.append(f"  WARNING: Unknown value '{v}' in '{cov_name}' -> NA")
        out_rows.append((cov_name, encoded))
        continue

    if cov_name in NUMERIC_COLS or cov_lower in NUMERIC_COLS:
        encoded = []
        for v in values:
            try:
                encoded.append(str(float(v.strip())))
            except ValueError:
                encoded.append('NA')
        out_rows.append((cov_name, encoded))
        continue

    if cov_name in ONEHOT_COLS or cov_lower in ONEHOT_COLS:
        # Normalise values first to collapse duplicates from typos/capitalisation
        norm_values = [normalise_level(cov_lower, v) for v in values]
        levels = sorted(set(
            v for v in norm_values
            if v not in ('na', '', 'nan')
        ))
        ref = levels[0]
        for level in levels[1:]:
            dummy = f"{cov_name}_{level.replace(' ','_').replace('/','_').replace(',','')}"
            encoded = ['1' if v == level else ('NA' if v in ('na','','nan') else '0')
                       for v in norm_values]
            out_rows.append((dummy, encoded))
        print(f"  One-hot: '{cov_name}' -> {len(levels)-1} dummies (ref='{ref}')", flush=True)
        continue

    # Unknown column: attempt numeric coercion
    numeric_vals = []
    n_valid = 0
    for v in values:
        try:
            numeric_vals.append(str(float(v.strip())))
            n_valid += 1
        except ValueError:
            numeric_vals.append('NA')
    if n_valid > 0:
        out_rows.append((cov_name, numeric_vals))
        print(f"  Numeric coercion: '{cov_name}' ({n_valid}/{len(values)} valid)", flush=True)
    else:
        dropped.append(cov_name)
        warnings.append(f"  DROPPED non-numeric col: '{cov_name}'")

with open(out_file, 'w') as f:
    f.write('covariate\t' + '\t'.join(sample_ids) + '\n')
    for row_name, vals in out_rows:
        f.write(row_name + '\t' + '\t'.join(vals) + '\n')

print(f"  Dropped            : {dropped}", flush=True)
for w in warnings:
    print(w, flush=True)
print(f"  Output covariates  : {len(out_rows)} rows x {len(sample_ids)} samples", flush=True)
print(f"  Written to         : {out_file}", flush=True)
EOF

echo "  Covariate encoding done."

##############################################
# STEP 7: SNP location file
# SNP names MUST match the row labels in the SNP matrix exactly.
# The .raw file names SNPs in PLINK's SNPID_ALLELE format (e.g. ._AC).
# We read the SNP names from the .raw header and pair them positionally
# with chr/pos from the .bim file (same row order guaranteed by PLINK).
##############################################
echo ""
echo "[7] Generating SNP location file..."

SNP_LOC=$OUTDIR/snp_location.txt

python3 << EOF
raw_file = "$RAW"
bim_file = "$BIM"
out_file = "$SNP_LOC"

# Read SNP names from .raw header (columns 7 onward, space-delimited)
with open(raw_file) as f:
    raw_header = f.readline().split()
snp_names = raw_header[6:]   # skip FID IID PAT MAT SEX PHENOTYPE
print(f"  SNP names from .raw header: {len(snp_names)}", flush=True)

# Read chr/pos from .bim (col 1=chr, col 4=pos), same row order as .raw SNPs
bim_positions = []
with open(bim_file) as f:
    for line in f:
        fields = line.split()
        bim_positions.append(('chr' + fields[0], fields[3]))
print(f"  Positions from .bim       : {len(bim_positions)}", flush=True)

if len(snp_names) != len(bim_positions):
    print(f"  ERROR: SNP count mismatch between .raw ({len(snp_names)}) and .bim ({len(bim_positions)})", flush=True)
    exit(1)

with open(out_file, 'w') as out:
    out.write('snpid\tchr\tpos\n')
    for name, (chrom, pos) in zip(snp_names, bim_positions):
        out.write(f'{name}\t{chrom}\t{pos}\n')

print(f"  Written {len(snp_names)} SNPs to: {out_file}", flush=True)
EOF

echo "  SNP location file done."

##############################################
# STEP 8: Gene location file from gencode v49 BED6
# BED6 field 4 format: ENSG00000123.7___gene_type___gene_name
# Strip ___* suffix to recover versioned ENSG ID.
# Versioned IDs match expression matrix directly (same gencode release).
##############################################
echo ""
echo "[8] Generating gene location file..."

GENE_LOC=$OUTDIR/gene_location.txt
echo -e "geneid\tchr\tleft\tright" > $GENE_LOC
awk 'BEGIN{OFS="\t"} {
    gene=$4; sub(/___.*$/, "", gene)
    print gene, $1, $2, $3
}' $GTF_BED6 >> $GENE_LOC

echo "  Genes in location file : $(tail -n +2 $GENE_LOC | wc -l)"
echo "  Written to: $GENE_LOC"

python3 << EOF
expr_file = "$OUTDIR/expression_${TISSUE_DIR}.txt"
loc_file  = "$OUTDIR/gene_location.txt"

with open(loc_file) as f:
    f.readline()
    loc_genes = set(line.split('\t')[0] for line in f)
with open(expr_file) as f:
    f.readline()
    expr_genes = [line.split('\t')[0] for line in f]

found   = sum(1 for g in expr_genes if g in loc_genes)
missing = [g for g in expr_genes if g not in loc_genes]
print(f"  Expression genes            : {len(expr_genes)}", flush=True)
print(f"  Found in gene location file : {found}", flush=True)
print(f"  Missing                     : {len(missing)}", flush=True)
if missing:
    print(f"  First 5 missing             : {missing[:5]}", flush=True)
EOF

##############################################
# STEP 9: Cleanup temp files
##############################################
echo ""
echo "[9] Cleaning up temp files..."
rm -f $TMP_META $TMP_WGS $TMP_MATCHED $TMP_HRA $TMP_HRA_UNDER $TMP_HDA
echo "  Done."

##############################################
# STEP 10: Final summary
##############################################
echo ""
echo "============================================"
echo "  Prep complete for  : $TISSUE"
echo "  Output directory   : $OUTDIR"
echo ""
echo "  snp_${TISSUE_DIR}.txt                  -> args[1] SNP_file_name"
echo "  expression_${TISSUE_DIR}.txt           -> args[2] expression_file_name"
echo "  covariates_${TISSUE_DIR}_encoded.txt   -> args[3] covariates_file_name"
echo "  gene_location.txt                      -> args[5] gene_location_file_name"
echo "  snp_location.txt                       -> args[6] snp_location_file_name"
echo ""
echo "  Next step: sbatch run_eQTL.sh \"$TISSUE\""
echo "  $(date)"
echo "============================================"
