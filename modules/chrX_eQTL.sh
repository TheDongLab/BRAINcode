#!/bin/bash
#SBATCH --job-name=chrX_eQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/chrX_eQTL_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/chrX_eQTL_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G

set -euo pipefail
module load R

# ── Arguments & Setup ──────────────────────────────────────────────────
if [ $# -lt 1 ]; then
    echo "ERROR: Missing tissue argument."
    echo "Usage: sbatch chrX_eQTL.sh \"Cerebellum\""
    exit 1
fi

TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

# ── Paths ─────────────────────────────────────────────────────────────
BASE=/home/zw529/donglab/data/target_ALS
ORIG_DIR=$BASE/$TISSUE_DIR/eQTL

SNP_LOC=$ORIG_DIR/snp_location.txt
GENE_LOC=$ORIG_DIR/gene_location.txt

# Create perfectly flat folders inside the primary eQTL directory
MALE_DIR=$ORIG_DIR/Male_ChrX
FEMALE_DIR=$ORIG_DIR/Female_ChrX

rm -rf "$MALE_DIR" "$FEMALE_DIR"
mkdir -p "$MALE_DIR" "$FEMALE_DIR"

echo "======================================================="
echo "  Prepping Flat Sex-Stratified ChrX Run for: $TISSUE"
echo "  Source Workspace: $ORIG_DIR"
echo "  Target Male Directory: $MALE_DIR"
echo "  Target Female Directory: $FEMALE_DIR"
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

sex_row <- as.numeric(cov["sex_bin", ])
names(sex_row) <- colnames(cov)

# 0 is Female, 1 is Male based on verification
males <- names(sex_row)[which(sex_row == 1)]
females <- names(sex_row)[which(sex_row == 0)]

males <- males[!is.na(males) & males != "NA" & males != ""]
females <- females[!is.na(females) & females != "NA" & females != ""]

write.table(males, file.path(orig_dir, "male_ids.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(females, file.path(orig_dir, "female_ids.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
EOF

# ── Step 3: Extract Matrices Using Memory-Efficient Streaming ─────────
echo "[3] Slicing matrices via memory-safe line streaming..."
python3 - << 'EOF'
import os
import pandas as pd

orig_dir = os.environ["ORIG_DIR"]
m_dir = os.environ["MALE_DIR"]
f_dir = os.environ["FEMALE_DIR"]
tissue_dir = os.environ["TISSUE_DIR"]

# 1. Parse IDs and variants
with open(f"{orig_dir}/male_ids.txt") as f: 
    males = [line.strip() for line in f if line.strip()]
with open(f"{orig_dir}/female_ids.txt") as f: 
    females = [line.strip() for line in f if line.strip()]
with open(f"{orig_dir}/chrx_variants_list.txt") as f: 
    chrx_snps = set(line.strip() for line in f if line.strip())

# 2. Process Covariates (Small file, fine for memory)
cov = pd.read_csv(f"{orig_dir}/covariates_{tissue_dir}_encoded.txt", sep="\t", index_col=0, keep_default_na=False, dtype=str)
cov.columns = cov.columns.str.strip()

# Corrected to drop 'sex_bin'
cov_clean = cov.drop(index="sex_bin", errors="ignore")

males_present = [m for m in males if m in cov_clean.columns]
females_present = [f for f in females if f in cov_clean.columns]

print(f"Verified columns -> Males: {len(males_present)}, Females: {len(females_present)}")

cov_clean[males_present].to_csv(f"{m_dir}/covariates_Male_ChrX_encoded.txt", sep="\t")
cov_clean[females_present].to_csv(f"{f_dir}/covariates_Female_ChrX_encoded.txt", sep="\t")

# 3. Process Expression
expr = pd.read_csv(f"{orig_dir}/expression_{tissue_dir}.txt", sep="\t", index_col=0, keep_default_na=False, dtype=str)
expr.columns = expr.columns.str.strip()
expr[males_present].to_csv(f"{m_dir}/expression_Male_ChrX.txt", sep="\t")
expr[females_present].to_csv(f"{f_dir}/expression_Female_ChrX.txt", sep="\t")
del expr

# 4. Stream Genotypes Line-by-Line (Memory-Safe)
snp_in_path = f"{orig_dir}/snp_{tissue_dir}.txt"
m_out_path = f"{m_dir}/snp_Male_ChrX.txt"
f_out_path = f"{f_dir}/snp_Female_ChrX.txt"

with open(snp_in_path, "r") as infile, \
     open(m_out_path, "w") as m_out, \
     open(f_out_path, "w") as f_out:
    
    header = infile.readline().strip().split("\t")
    header = [h.strip() for h in header]
    
    m_indices = [header.index(m) for m in males_present if m in header]
    f_indices = [header.index(f) for f in females_present if f in header]
    
    m_out.write("snpid\t" + "\t".join(males_present) + "\n")
    f_out.write("snpid\t" + "\t".join(females_present) + "\n")
    
    for line in infile:
        parts = line.strip().split("\t")
        if not parts:
            continue
        snp_id = parts[0].strip()
        
        if snp_id in chrx_snps:
            m_vals = [parts[idx] for idx in m_indices]
            f_vals = [parts[idx] for idx in f_indices]
            
            m_out.write(f"{snp_id}\t" + "\t".join(m_vals) + "\n")
            f_out.write(f"{snp_id}\t" + "\t".join(f_vals) + "\n")
EOF

# ── Step 4: Submit the QTL Jobs Natively (With Sample Guards) ─────────
echo "[4] Dispatching stratified jobs where cohorts are present..."

M_COUNT=$(awk 'END{print NR}' "$ORIG_DIR/male_ids.txt")
F_COUNT=$(awk 'END{print NR}' "$ORIG_DIR/female_ids.txt")

if [ "$M_COUNT" -gt 0 ] && [ -f "$MALE_DIR/snp_Male_ChrX.txt" ]; then
    echo "-> Submitting Male Run ($M_COUNT samples)..."
    sbatch run_eQTL.sh "$TISSUE" "standard" "Male_ChrX"
else
    echo "-> Skipping Male Run: No samples found for this tissue."
fi

if [ "$F_COUNT" -gt 0 ] && [ -f "$FEMALE_DIR/snp_Female_ChrX.txt" ]; then
    echo "-> Submitting Female Run ($F_COUNT samples)..."
    sbatch run_eQTL.sh "$TISSUE" "standard" "Female_ChrX"
else
    echo "-> Skipping Female Run: No samples found for this tissue."
fi

echo "======================================================="
echo " Complete! Dispatches evaluated."
echo "======================================================="
