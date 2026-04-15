#!/bin/bash
#SBATCH --job-name=prep_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=499G

set -euo pipefail

### Arguments ###
if [ $# -lt 1 ]; then
    echo "ERROR: Usage: sbatch prep_eQTL.sh \"Tissue Name\""
    exit 1
fi

TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

# Writing directly to the tissue's eQTL directory
OUTDIR=/home/zw529/donglab/data/target_ALS/$TISSUE_DIR/eQTL
mkdir -p $OUTDIR

### Verified Paths from ls -lA ###
BASE=/home/zw529/donglab/data/target_ALS
QTL_DIR=$BASE/QTL
PLINK=$QTL_DIR/plink
REFS=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode

# Input Files
EXPR=$QTL_DIR/expression_matrix.txt
META=$QTL_DIR/expression_sample_metadata.csv
SAMPLE_MAP=$QTL_DIR/sample_map.txt        # Updated from your ls
COV=$QTL_DIR/covariates.tsv               # Updated from your ls
RAW=$PLINK/joint_autosomes_matrixEQTL.raw
BIM=$PLINK/joint_autosomes_filtered_bed.bim
PCA=$PLINK/joint_pca.eigenvec
GTF_BED6=$REFS/gencode.v49.annotation.gene.bed6

# Temp files for alignment
TMP_MATCHED=$OUTDIR/tmp_matched.txt
TMP_HDA=$OUTDIR/tmp_HDA_ids.txt
CLEAN_COV=$OUTDIR/covariates_filtered.tsv

echo "============================================"
echo "  eQTL Prep for $TISSUE"
echo "  $(date)"
echo "============================================"

##############################################
# STEP 1: Build sample mapping
##############################################
echo "[1] Filtering metadata and building ID maps..."

python3 << EOF
import pandas as pd
df = pd.read_csv("$COV", sep='\t')

# Quality Filter: Drop samples with low RIN and high RNA skew
bad_samples = df[df['rin'].isna() & (df['rna_skew'] > 1.0)]['externalsampleid'].tolist()
df = df[~df['externalsampleid'].isin(bad_samples)]

# Case/Control Encoding
def collapse(g):
    if 'ALS' in str(g): return 1
    if 'Non-Neurological' in str(g): return 0
    return None

df['is_als'] = df['subject_group'].apply(collapse)
df = df[df['is_als'].notna()]
df.to_csv("$CLEAN_COV", sep='\t', index=False)

with open("$OUTDIR/keep_samples.txt", 'w') as f:
    for s in df['externalsampleid'].tolist(): f.write(s + "\n")
EOF

# Mapping HRA (RNA) to HDA (DNA) using sample_map.txt
# sample_map.txt format: SubjectID,HRA_ID,HDA_ID
tail -n +2 $META | awk -F',' -v tissue="$TISSUE" '$3 == tissue {print $2","$1}' | grep -Ff $OUTDIR/keep_samples.txt | sort -t',' -k1,1 > $OUTDIR/tmp_meta.txt
tail -n +2 $SAMPLE_MAP | awk -F',' '{print $1","$2","$3}' | sort -t',' -k1,1 > $OUTDIR/tmp_map_sorted.txt

# Join on SubjectID (Column 1)
join -t',' -1 1 -2 1 $OUTDIR/tmp_meta.txt $OUTDIR/tmp_map_sorted.txt > $TMP_MATCHED
# Extract DNA IDs (HDA)
awk -F',' '{print $4}' $TMP_MATCHED > $TMP_HDA
# Extract RNA IDs with underscore fix (HRA_001)
awk -F',' '{print $2}' $TMP_MATCHED | sed 's/-/_/g' > $OUTDIR/tmp_hra_under.txt

##############################################
# STEP 2: Subset Expression
##############################################
echo "[2] Subsetting expression matrix..."
python3 << EOF
ids_file, expr_file, out_file = "$OUTDIR/tmp_hra_under.txt", "$EXPR", "$OUTDIR/expression_${TISSUE_DIR}.txt"
with open(ids_file) as f: keep_ids = set(line.strip() for line in f if line.strip())
with open(expr_file) as f, open(out_file, 'w') as out:
    header = f.readline().rstrip('\n').split('\t')
    keep_idx = [0] + [i for i, h in enumerate(header) if h in keep_ids]
    out.write('\t'.join(header[i] for i in keep_idx) + '\n')
    for line in f:
        fields = line.rstrip('\n').split('\t')
        out.write('\t'.join(fields[i] for i in keep_idx) + '\n')
EOF

##############################################
# STEP 3: Subset + Transpose SNPs
##############################################
echo "[3] Processing SNPs (MAF > 0.05)..."
python3 << EOF
import re, numpy as np
ids_file, raw_file, bim_file = "$TMP_HDA", "$RAW", "$BIM"
out_file = "$OUTDIR/snp_${TISSUE_DIR}.txt"

with open(ids_file) as f: keep_ids = set(line.strip() for line in f if line.strip())
chrpos_names = []
with open(bim_file) as f:
    for line in f:
        p = line.split()
        if len(p) >= 4: chrpos_names.append(f"chr{p[0]}:{p[3]}")

sample_ids, rows = [], []
with open(raw_file) as f:
    f.readline() # Skip header
    n_snps = len(chrpos_names)
    for line in f:
        fields = line.split()
        if len(fields) < 6: continue
        iid = re.sub(r'-b\d+$', '', fields[1])
        if iid in keep_ids:
            sample_ids.append(iid)
            geno_fields = fields[6:6+n_snps]
            rows.append([float(v) if v != "NA" else np.nan for v in geno_fields])

geno = np.array(rows, dtype=float)
af = np.nanmean(geno, axis=0) / 2.0
maf = np.minimum(af, 1 - af)
keep_idx = np.where(np.nan_to_num(maf) >= 0.05)[0]

with open(out_file, 'w') as out:
    out.write('snpid\t' + '\t'.join(sample_ids) + '\n')
    for idx in keep_idx:
        vals = [str(int(v)) if not np.isnan(v) else "NA" for v in geno[:, idx]]
        out.write(chrpos_names[idx] + '\t' + '\t'.join(vals) + '\n')
EOF

##############################################
# STEP 4-5: Alignment & HDA Renaming
##############################################
echo "[4-5] Aligning RNA and DNA IDs..."
python3 << EOF
import pandas as pd
# tmp_matched format: Subj,HRA,RNA_ID,HDA
mapping = pd.read_csv("$TMP_MATCHED", header=None, names=['sub', 'hra', 'rna_id', 'hda'])
hda_to_hra = dict(zip(mapping['hda'], mapping['hra']))

with open("$OUTDIR/snp_${TISSUE_DIR}.txt", 'r') as f:
    snp_hda_samples = f.readline().strip().split('\t')[1:]

target_hra_expr = [hda_to_hra[h].replace('-', '_') for h in snp_hda_samples]
target_hra_cov = [hda_to_hra[h] for h in snp_hda_samples]

df_expr = pd.read_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t', index_col=0)
df_expr = df_expr[target_hra_expr]
df_expr.columns = snp_hda_samples
df_expr.to_csv("$OUTDIR/expression_${TISSUE_DIR}.txt", sep='\t')

df_cov_raw = pd.read_csv("$CLEAN_COV", sep='\t', index_col=0)
df_cov_aligned = df_cov_raw.loc[target_hra_cov]
df_cov_aligned.index = snp_hda_samples
df_cov_aligned.to_csv("$OUTDIR/covariates_${TISSUE_DIR}.txt", sep='\t')
EOF

##############################################
# STEP 6: Encode Covariates (Sex + PCA)
##############################################
echo "[6] Encoding Covariates and Adding PCA..."
python3 << EOF
import pandas as pd
df = pd.read_csv("$OUTDIR/covariates_${TISSUE_DIR}.txt", sep='\t', index_col=0)
pca = pd.read_csv("$PCA", sep='\s+', header=None)
pca = pca.iloc[:, 1:7]
pca.columns = ['ID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
pca.set_index('ID', inplace=True)

df['sex_binary'] = df['sex'].str.lower().map({'male': 1, 'female': 0})
cols_to_keep = ['sex_binary', 'age_at_death', 'rin', 'rna_skew', 'is_als']
cov_matrix = df[cols_to_keep].join(pca, how='inner')
cov_matrix.T.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", sep='\t')
EOF

##############################################
# STEP 7-8: Location Files & Validation
##############################################
echo "[7] Generating SNP location file..."
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk 'BEGIN{OFS="\t"} {print "chr"$1":"$4, "chr"$1, $4}' $BIM >> $OUTDIR/snp_location.txt

echo "[8] Generating gene location file..."
echo -e "geneid\tchr\tleft\tright" > $OUTDIR/gene_location.txt
awk 'BEGIN{OFS="\t"} {
    gene=$4; sub(/___.*$/, "", gene); print gene, $1, $2, $3
}' $GTF_BED6 >> $OUTDIR/gene_location.txt

python3 << EOF
expr_file = "$OUTDIR/expression_${TISSUE_DIR}.txt"
loc_file  = "$OUTDIR/gene_location.txt"

with open(loc_file) as f:
    f.readline()
    loc_genes = set(line.split('\t')[0] for line in f)
with open(expr_file) as f:
    f.readline()
    expr_genes = [line.split('\t')[0] for line in f]

found = sum(1 for g in expr_genes if g in loc_genes)
print(f"  Validation: Found location for {found} / {len(expr_genes)} expression genes.")
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
echo "  Prep complete for  : $TISSUE ($STRAT_SEX)"
echo "  Output directory   : $OUTDIR"
echo ""
echo "  snp_${TISSUE_DIR}.txt                  -> args[1] SNP_file_name"
echo "  expression_${TISSUE_DIR}.txt           -> args[2] expression_file_name"
echo "  covariates_${TISSUE_DIR}_encoded.txt   -> args[3] covariates_file_name"
echo "  gene_location.txt                      -> args[5] gene_location_file_name"
echo "  snp_location.txt                       -> args[6] snp_location_file_name"
echo ""
echo "  Next step: sbatch run_eQTL.sh \"$TISSUE\" \"$STRAT_SEX\""
echo "  $(date)"
echo "============================================"
