#!/bin/bash
#SBATCH --job-name=prep_sQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.err
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=499G

set -euo pipefail
module load R

### Arguments: ###
if [ $# -lt 2 ]; then
    echo "ERROR: Missing arguments."
    echo "Usage: sbatch prep_sQTL.sh \"Cerebellum\" \"Male\""
    exit 1
fi

TISSUE="$1"
STRAT_SEX="$2"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

echo "============================================"
echo "  sQTL Prep for tissue : $TISSUE"
echo "  Stratification       : $STRAT_SEX"
echo "  $(date)"
echo "============================================"

### Paths: ###
BASE=/home/zw529/donglab/data/target_ALS
PLINK=$BASE/QTL/plink
SPLICING_RAW=$BASE/QTL/splicing_matrix.txt 
SPLICING_LOC_SRC=$BASE/QTL/splicing_events_hg38.bed 

RAW=$PLINK/joint_autosomes_matrixEQTL.raw
META=$BASE/QTL/expression_sample_metadata.csv
WGS_MAP=$BASE/QTL/wgs_samples_for_vcf_merge.csv
COV=$BASE/QTL/covariates.tsv
BIM=$PLINK/joint_autosomes_filtered_bed.bim

OUTDIR=$BASE/$TISSUE_DIR/sQTL/$STRAT_SEX
mkdir -p $OUTDIR

# Temp Files
TMP_META=$OUTDIR/tmp_meta.txt
TMP_WGS=$OUTDIR/tmp_wgs_map.txt
TMP_MATCHED=$OUTDIR/tmp_matched.txt
TMP_HRA=$OUTDIR/tmp_HRA_ids.txt
TMP_HRA_UNDER=$OUTDIR/tmp_HRA_ids_underscore.txt
TMP_HDA=$OUTDIR/tmp_HDA_ids.txt
CLEAN_COV=$OUTDIR/covariates_stratified.tsv

##############################################
# STEP 0: Generate BED from Matrix IDs
##############################################
echo "[0] Generating BED file from matrix IDs..."
# Fixed the awk warning by removing the backslash
tail -n +2 "$SPLICING_RAW" | cut -f1 | awk -F'[:-]' 'BEGIN{OFS="\t"} {print $1, $3, $4, $0}' > "$SPLICING_LOC_SRC"

##############################################
# STEP 1: Build sample mapping
##############################################
echo "[1] Building sample mapping..."
python3 << EOF
import pandas as pd
# Added encoding='latin1' to handle non-utf8 characters
df = pd.read_csv("$COV", sep='\t', encoding='latin1')
bad_samples = df[df['rin'].isna() & (df['rna_skew'] > 1.0)]['externalsampleid'].tolist()
df = df[~df['externalsampleid'].isin(bad_samples)]
df = df[df['sex'] == "$STRAT_SEX"]
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

tail -n +2 $META | awk -F',' -v tissue="$TISSUE" '$3 == tissue {print $2","$1}' | grep -Ff $OUTDIR/keep_samples.txt | sort -t',' -k1,1 > $TMP_META
tail -n +2 $WGS_MAP | awk -F',' '{print $1","$2}' | sort -t',' -k1,1 > $TMP_WGS
join -t',' -1 1 -2 1 $TMP_META $TMP_WGS > $TMP_MATCHED
awk -F',' '{print $2}' $TMP_MATCHED > $TMP_HRA
awk -F',' '{print $3}' $TMP_MATCHED > $TMP_HDA
sed 's/-/_/g' $TMP_HRA > $TMP_HRA_UNDER

##############################################
# STEP 2: Subset Splicing Matrix
##############################################
echo "[2] Processing Splicing Matrix..."
Rscript - << EOF
library(data.table)
keep_ids <- fread("$TMP_HRA_UNDER", header=FALSE)\$V1
full_data <- fread("$SPLICING_RAW", header=TRUE)
cols_to_keep <- c(names(full_data)[1], intersect(names(full_data), keep_ids))
splicing <- full_data[, ..cols_to_keep]
rm(full_data); gc()
row_vars <- apply(splicing[, -1, with=FALSE], 1, function(x) var(x, na.rm=TRUE))
splicing <- splicing[which(row_vars > 1e-10), ]
fwrite(splicing, "$OUTDIR/splicing_${TISSUE_DIR}.txt", sep="\t", quote=FALSE)
EOF

##############################################
# STEP 3: Subset + Transpose SNPs
##############################################
echo "[3] Subsetting SNPs..."
python3 << EOF
import re, numpy as np, os
with open("$TMP_HDA") as f: keep_ids = set(line.strip() for line in f if line.strip())
chrpos_names = []
with open("$BIM") as f:
    for line in f:
        p = line.split()
        if len(p) >= 4: chrpos_names.append(f"chr{p[0]}:{p[3]}")

sample_ids, rows = [], []
with open("$RAW") as f:
    header = f.readline().split()
    n_snps = len(chrpos_names)
    for line in f:
        fields = line.split()
        if len(fields) < 6: continue
        iid = re.sub(r'-b\d+$', '', fields[1])
        if iid in keep_ids:
            sample_ids.append(iid)
            geno_fields = fields[6:6+n_snps]
            rows.append([float(v) if v!="NA" else np.nan for v in geno_fields])

geno = np.array(rows, dtype=float)
af = np.nanmean(geno, axis=0) / 2.0
maf = np.minimum(af, 1 - af)
keep_idx = np.where(np.nan_to_num(maf) >= 0.05)[0]

with open("$OUTDIR/snp_${TISSUE_DIR}.txt", 'w') as out:
    out.write('snpid\t' + '\t'.join(sample_ids) + '\n')
    for idx in keep_idx:
        vals = [str(int(v)) if not np.isnan(v) else "NA" for v in geno[:, idx]]
        out.write(chrpos_names[idx] + '\t' + '\t'.join(vals) + '\n')
EOF

##############################################
# STEP 4: Subset Covariates
##############################################
echo "[4] Subsetting Covariates..."
python3 << EOF
import pandas as pd
with open("$TMP_HRA") as f: keep_ids = set(line.strip() for line in f if line.strip())
# Encoding handled here as well
df_cov = pd.read_csv("$CLEAN_COV", sep='\t', encoding='latin1')
df_cov = df_cov[df_cov['externalsampleid'].isin(keep_ids)]
df_cov.to_csv("$OUTDIR/covariates_temp.txt", sep='\t', index=False)
EOF

##############################################
# STEP 5: ID Translation & Final Alignment
##############################################
echo "[5] Final Alignment..."
python3 << EOF
import pandas as pd
mapping = pd.read_csv("$TMP_MATCHED", header=None, names=['sub', 'hra', 'hda'])
hda_to_hra = dict(zip(mapping['hda'], mapping['hra']))
with open("$OUTDIR/snp_${TISSUE_DIR}.txt", 'r') as f:
    snp_hda_samples = f.readline().strip().split('\t')[1:]

target_hra_samples = [hda_to_hra[h].replace('-', '_') for h in snp_hda_samples]
df_splicing = pd.read_csv("$OUTDIR/splicing_${TISSUE_DIR}.txt", sep='\t', index_col=0)
df_splicing_aligned = df_splicing[target_hra_samples]
df_splicing_aligned.columns = snp_hda_samples
df_splicing_aligned.to_csv("$OUTDIR/splicing_${TISSUE_DIR}.txt", sep='\t')

target_hra_cov = [h.replace('_', '-') for h in target_hra_samples]
df_cov = pd.read_csv("$OUTDIR/covariates_temp.txt", sep='\t', index_col=0)
df_cov_aligned = df_cov.loc[target_hra_cov].T
df_cov_aligned.columns = snp_hda_samples
df_cov_aligned.to_csv("$OUTDIR/covariates_${TISSUE_DIR}.txt", sep='\t')
EOF

##############################################
# STEP 6: Encode Covariates
##############################################
echo "[6] Encoding Covariates..."
python3 << EOF
import csv
DROP = {'externalsampleid', 'externalsubjectid', 'tissue', 'subject_group', 'sex'}
NUMERIC = {'age_at_death', 'rin', 'rna_skew', 'is_als'}
with open("$OUTDIR/covariates_${TISSUE_DIR}.txt", encoding='latin1') as f:
    r = csv.reader(f, delimiter='\t'); h = next(r); rows = {row[0]: row[1:] for row in r}
out_rows = []
for name, vals in rows.items():
    if name.lower() in DROP: continue
    if name.lower() in NUMERIC:
        out_rows.append((name, [str(float(v)) if v and v!='NA' else 'NA' for v in vals]))
with open("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt", 'w') as f:
    f.write('covariate\t' + '\t'.join(h[1:]) + '\n')
    for n, v in out_rows: f.write(n + '\t' + '\t'.join(v) + '\n')
EOF

##############################################
# STEP 7: SNP location
##############################################
echo "[7] Generating SNP location file..."
echo -e "snpid\tchr\tpos" > $OUTDIR/snp_location.txt
awk 'BEGIN{OFS="\t"} {print "chr"$1":"$4, "chr"$1, $4}' $BIM >> $OUTDIR/snp_location.txt

##############################################
# STEP 8: Splicing location
##############################################
echo "[8] Generating splicing location file..."
echo -e "geneid\tchr\tleft\tright" > $OUTDIR/splicing_location.txt
awk 'BEGIN{OFS="\t"} {print $4, $1, $2, $3}' "$SPLICING_LOC_SRC" >> $OUTDIR/splicing_location.txt

##############################################
# STEP 9: Final Summary
##############################################
echo "============================================"
echo "  sQTL Prep complete for $TISSUE ($STRAT_SEX)"
echo "  $(date)"
echo "============================================"
