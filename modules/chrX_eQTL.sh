#!/bin/bash
#SBATCH --job-name=chrX_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/chrX_eQTL_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/chrX_eQTL_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

set -euo pipefail

module load R

# ── Arguments & Setup ──────────────────────────────────────────────────
if [ $# -lt 1 ]; then
    echo "ERROR: Missing tissue argument."
    echo "Usage: sbatch prep_chrx.sh \"Cerebellum\""
    exit 1
fi

TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

# ── Paths ─────────────────────────────────────────────────────────────
BASE=/home/zw529/donglab/data/target_ALS
ORIG_DIR=$BASE/$TISSUE_DIR/eQTL

SNP_LOC=$ORIG_DIR/snp_location.txt
GENE_LOC=$ORIG_DIR/gene_location.txt

# Perfectly clean folder matching your standard structure layouts
MALE_DIR=$ORIG_DIR/Male_ChrX
FEMALE_DIR=$ORIG_DIR/Female_ChrX

mkdir -p "$MALE_DIR" "$FEMALE_DIR"

echo "======================================================="
echo "  Prepping Sex-Stratified ChrX Run for: $TISSUE"
echo "  Source Directory: $ORIG_DIR"
echo "  Target Male Workspace: $MALE_DIR"
echo "  Target Female Workspace: $FEMALE_DIR"
echo "  $(date)"
echo "======================================================="

# ── Step 1: Subset SNP Location and Gene Location Files ──────────────
echo "[1] Extracting Chromosome X locations and variants..."

awk 'NR==1 || $2 == "chrX"' "$SNP_LOC" > "$MALE_DIR/snp_location.txt"
cp "$MALE_DIR/snp_location.txt" "$FEMALE_DIR/snp_location.txt"

awk 'NR==1 || $2 == "chrX"' "$GENE_LOC" > "$MALE_DIR/gene_location.txt"
cp "$MALE_DIR/gene_location.txt" "$FEMALE_DIR/gene_location.txt"

awk 'NR>1 {print $1}' "$MALE_DIR/snp_location.txt" > "$ORIG_DIR/chrx_variants_list.txt"

export ORIG_DIR MALE_DIR FEMALE_DIR TISSUE_DIR

# ── Step 2: Identify Samples by Sex via R ────────────────────────────
echo "[2] Parsing sample IDs by sex from covariates..."
Rscript - << 'EOF'
orig_dir <- Sys.getenv("ORIG_DIR")
tissue_dir <- Sys.getenv("TISSUE_DIR")

cov_file <- file.path(orig_dir, paste0("covariates_", tissue_dir, "_encoded.txt"))
cov <- read.table(cov_file, header=TRUE, sep="\t", row.names=1, check.names=FALSE, na.strings=c("", "NA", "NaN"))

sex_row <- as.numeric(cov["Sex", ])
names(sex_row) <- colnames(cov)

males <- names(sex_row)[which(sex_row == 1)]
females <- names(sex_row)[which(sex_row == 2)]

males <- males[!is.na(males) & males != "NA" & males != ""]
females <- females[!is.na(females) & females != "NA" & females != ""]

write.table(males, file.path(orig_dir, "male_ids.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(females, file.path(orig_dir, "female_ids.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
EOF

# ── Step 3: Extract Matrices Using Python ─────────────────────────────
echo "[3] Slicing expression, covariate, and genotype matrices into flat structures..."
python3 - << 'EOF'
import os
import pandas as pd

orig_dir = os.environ["ORIG_DIR"]
m_dir = os.environ["MALE_DIR"]
f_dir = os.environ["FEMALE_DIR"]
tissue_dir = os.environ["TISSUE_DIR"]

# Set up the internal slash subdirectories to handle run_eQTL.sh constraints cleanly
# INDIR evaluates as Male_ChrX/ because of the path workaround
os.makedirs(f"{m_dir}/snp_{tissue_dir}/eQTL", exist_ok=True)
os.makedirs(f"{m_dir}/expression_{tissue_dir}/eQTL", exist_ok=True)
os.makedirs(f"{m_dir}/covariates_{tissue_dir}", exist_ok=True)

os.makedirs(f"{f_dir}/snp_{tissue_dir}/eQTL", exist_ok=True)
os.makedirs(f"{f_dir}/expression_{tissue_dir}/eQTL", exist_ok=True)
os.makedirs(f"{f_dir}/covariates_{tissue_dir}", exist_ok=True)

with open(f"{orig_dir}/male_ids.txt") as f: males = [line.strip() for line in f]
with open(f"{orig_dir}/female_ids.txt") as f: females = [line.strip() for line in f]
with open(f"{orig_dir}/chrx_variants_list.txt") as f: chrx_snps = set(line.strip() for line in f)

# 1. Process Covariates 
cov = pd.read_csv(f"{orig_dir}/covariates_{tissue_dir}_encoded.txt", sep="\t", index_col=0, keep_default_na=False, dtype=str)
cov_clean = cov.drop(index="Sex", errors="ignore")
cov_clean[males].to_csv(f"{m_dir}/covariates_{tissue_dir}/eQTL_Male_ChrX_encoded.txt", sep="\t")
cov_clean[females].to_csv(f"{f_dir}/covariates_{tissue_dir}/eQTL_Female_ChrX_encoded.txt", sep="\t")

# 2. Process Expression
expr = pd.read_csv(f"{orig_dir}/expression_{tissue_dir}.txt", sep="\t", index_col=0, keep_default_na=False, dtype=str)
expr[males].to_csv(f"{m_dir}/expression_{tissue_dir}/eQTL_Male_ChrX.txt", sep="\t")
expr[females].to_csv(f"{f_dir}/expression_{tissue_dir}/eQTL_Female_ChrX.txt", sep="\t")
del expr

# 3. Process SNPs
snp = pd.read_csv(f"{orig_dir}/snp_{tissue_dir}.txt", sep="\t", index_col=0, keep_default_na=False, dtype=str)
chrx_snp = snp.loc[snp.index.isin(chrx_snps)]
chrx_snp[males].to_csv(f"{m_dir}/snp_{tissue_dir}/eQTL_Male_ChrX.txt", sep="\t")
chrx_snp[females].to_csv(f"{f_dir}/snp_{tissue_dir}/eQTL_Female_ChrX.txt", sep="\t")
EOF

# ── Step 4: Submit the QTL Jobs via standard script ──────────────────
echo "[4] Submitting sex-stratified Matrix eQTL jobs..."

# The custom trail path cleanly resolves back up to make the output folders perfectly flat
sbatch run_eQTL.sh "${TISSUE}/eQTL/Male_ChrX/.."
sbatch run_eQTL.sh "${TISSUE}/eQTL/Female_ChrX/.."

echo "======================================================="
echo " Complete! Jobs dispatched flat inside $MALE_DIR and $FEMALE_DIR"
echo "======================================================="
