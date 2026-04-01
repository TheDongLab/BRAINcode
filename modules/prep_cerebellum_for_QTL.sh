#!/bin/bash
#SBATCH --job-name=prep_cerebellum
#SBATCH --output=/home/zw529/donglab/data/target_ALS/eQTL/prep_cerebellum.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/eQTL/prep_cerebellum.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=499G

###########################################
# prep_cerebellum.sh
# Purpose: Subset all Matrix eQTL inputs to Cerebellum samples only,
#           and prepare gene/SNP location files.
#
# Logic: Keeps subjects who have BOTH:
#   - Cerebellum RNAseq in expression_sample_metadata.csv
#   - WGS data (any tissue) in wgs_samples_for_vcf_merge.csv
#
# Usage: sbatch prep_cerebellum.sh
#        or:   bash prep_cerebellum.sh
# Output: $BASE/cerebellum/
###########################################

set -euo pipefail

BASE=/home/zw529/donglab/data/target_ALS
PLINK=$BASE/eQTL/plink
REFS=/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode

EXPR=$BASE/eQTL/expression_matrix.txt
RAW=$PLINK/joint_autosomes_matrixEQTL.raw
META=$BASE/eQTL/expression_sample_metadata.csv          # HRA <-> subjectid <-> tissue
WGS_MAP=$BASE/eQTL/wgs_samples_for_vcf_merge.csv        # subjectid <-> HDA
COV=$BASE/eQTL/covariates.tsv
BIM=$PLINK/joint_autosomes_filtered_bed.bim
GTF_BED6=$REFS/gencode.v49.annotation.gene.bed6

OUTDIR=$BASE/Cerebellum/eQTL
mkdir -p $OUTDIR

echo "============================================"
echo "  Cerebellum eQTL Prep"
echo "  $(date)"
echo "============================================"

##############################################
# STEP 1: Build sample mapping table
# Join metadata (HRA + subjectid) with WGS map (subjectid + HDA)
# Keeps subjects with BOTH Cerebellum RNAseq AND WGS (any tissue)
##############################################
echo ""
echo "[1] Building Cerebellum sample mapping (HRA <-> subjectid <-> HDA)..."

# Filter metadata to Cerebellum → subjectid, HRA_sampleid
tail -n +2 $META \
    | awk -F',' '$3 == "Cerebellum" {print $2","$1}' \
    | sort -t',' -k1,1 \
    > /tmp/cerebellum_meta.txt

N_META=$(wc -l < /tmp/cerebellum_meta.txt)
echo "  Cerebellum RNA samples in metadata : $N_META"

# WGS map: subjectid -> HDA_sampleid
tail -n +2 $WGS_MAP \
    | awk -F',' '{print $1","$2}' \
    | sort -t',' -k1,1 \
    > /tmp/wgs_map.txt

N_WGS=$(wc -l < /tmp/wgs_map.txt)
echo "  Subjects with WGS data             : $N_WGS"

# Inner join on subjectid → subjectid, HRA_id, HDA_id
join -t',' -1 1 -2 1 \
    /tmp/cerebellum_meta.txt \
    /tmp/wgs_map.txt \
    > /tmp/cerebellum_matched.txt

N_MATCHED=$(wc -l < /tmp/cerebellum_matched.txt)
echo "  Subjects with Cerebellum RNA + WGS : $N_MATCHED"

if [ "$N_MATCHED" -eq 0 ]; then
    echo "ERROR: No matched samples found. Check subjectid formatting between metadata and WGS map."
    exit 1
fi

# Extract ID lists
awk -F',' '{print $2}' /tmp/cerebellum_matched.txt > /tmp/cerebellum_HRA_ids.txt
awk -F',' '{print $3}' /tmp/cerebellum_matched.txt > /tmp/cerebellum_HDA_ids.txt

echo "  Sample HRA IDs (expression): $(head -3 /tmp/cerebellum_HRA_ids.txt | tr '\n' ' ')..."
echo "  Sample HDA IDs (SNP):        $(head -3 /tmp/cerebellum_HDA_ids.txt | tr '\n' ' ')..."

##############################################
# STEP 2: Subset expression matrix
# Expression matrix uses underscores (CGND_HRA_03594)
# Metadata uses hyphens (CGND-HRA-03594) — normalize before matching
##############################################
echo ""
echo "[2] Subsetting expression matrix to Cerebellum samples..."

sed 's/-/_/g' /tmp/cerebellum_HRA_ids.txt > /tmp/cerebellum_HRA_ids_underscore.txt

python3 << 'EOF'
ids_file = "/tmp/cerebellum_HRA_ids_underscore.txt"
expr_file = "/home/zw529/donglab/data/target_ALS/eQTL/expression_matrix.txt"
out_file  = "/home/zw529/donglab/data/target_ALS/eQTL/cerebellum/expression_cerebellum.txt"

with open(ids_file) as f:
    keep_ids = set(line.strip() for line in f if line.strip())

with open(expr_file) as f, open(out_file, 'w') as out:
    header = f.readline().rstrip('\n').split('\t')
    keep_idx = [0] + [i for i, h in enumerate(header) if h in keep_ids]
    found_ids = [header[i] for i in keep_idx[1:]]
    missing = keep_ids - set(found_ids)

    print(f"  Requested samples : {len(keep_ids)}", flush=True)
    print(f"  Found in matrix   : {len(found_ids)}", flush=True)
    if missing:
        print(f"  WARNING - {len(missing)} IDs not found in expression matrix", flush=True)
        print(f"  First 5 missing   : {list(missing)[:5]}", flush=True)

    out.write('\t'.join(header[i] for i in keep_idx) + '\n')
    for line in f:
        fields = line.rstrip('\n').split('\t')
        out.write('\t'.join(fields[i] for i in keep_idx) + '\n')

print(f"  Written to: {out_file}", flush=True)
EOF

echo "  Expression matrix done."

##############################################
# STEP 3: Subset covariates
# Auto-detects wide (samples as cols) vs long (samples as rows) format
##############################################
echo ""
echo "[3] Subsetting covariates to Cerebellum samples..."

python3 << 'EOF'
ids_file = "/tmp/cerebellum_HRA_ids.txt"   # hyphens, as in covariates.tsv
cov_file = "/home/zw529/donglab/data/target_ALS/eQTL/covariates.tsv"
out_file = "/home/zw529/donglab/data/target_ALS/eQTL/cerebellum/covariates_cerebellum.txt"

with open(ids_file) as f:
    keep_ids = set(line.strip() for line in f if line.strip())

with open(cov_file) as f:
    header = f.readline().rstrip('\n').split('\t')

sample_in_header = any(h in keep_ids for h in header)

with open(cov_file) as f, open(out_file, 'w') as out:
    header = f.readline().rstrip('\n').split('\t')
    if sample_in_header:
        # Wide format: covariates as rows, samples as columns
        keep_idx = [0] + [i for i, h in enumerate(header) if h in keep_ids]
        print(f"  Covariate format  : WIDE (samples as columns)", flush=True)
        print(f"  Keeping {len(keep_idx)-1} sample columns", flush=True)
        out.write('\t'.join(header[i] for i in keep_idx) + '\n')
        for line in f:
            fields = line.rstrip('\n').split('\t')
            out.write('\t'.join(fields[i] for i in keep_idx) + '\n')
    else:
        # Long format: samples as rows
        print(f"  Covariate format  : LONG (samples as rows)", flush=True)
        out.write('\t'.join(header) + '\n')
        n = 0
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if fields[0] in keep_ids:
                out.write(line)
                n += 1
        print(f"  Kept {n} sample rows", flush=True)

print(f"  Written to: {out_file}", flush=True)
EOF

echo "  Covariates done."

##############################################
# STEP 4: Subset and transpose SNP matrix
# .raw is space-delimited, samples are ROWS, SNPs are COLUMNS
# Matrix eQTL wants SNPs as ROWS, samples as COLUMNS
# HDA IDs in .raw may have a -b38 suffix — strip it for matching
##############################################
echo ""
echo "[4] Subsetting and transposing SNP matrix (chunked to manage memory)..."

python3 << 'EOF'
import re

ids_file = "/tmp/cerebellum_HDA_ids.txt"
raw_file = "/home/zw529/donglab/data/target_ALS/eQTL/plink/joint_autosomes_matrixEQTL.raw"
out_file = "/home/zw529/donglab/data/target_ALS/eQTL/cerebellum/snp_cerebellum.txt"

with open(ids_file) as f:
    keep_ids = set(line.strip() for line in f if line.strip())

print(f"  HDA IDs to keep: {len(keep_ids)}", flush=True)

def normalize_hda(iid):
    # Strip trailing -b38 / -b37 suffixes added by PLINK
    return re.sub(r'-b\d+$', '', iid)

snp_names = None
rows = []   # list of (IID, [genotype values...])

with open(raw_file) as f:
    header = f.readline().split()
    snp_names = header[6:]

    for line in f:
        fields = line.split()
        iid_raw = fields[1]
        iid_norm = normalize_hda(iid_raw)
        matched_id = iid_norm if iid_norm in keep_ids else (iid_raw if iid_raw in keep_ids else None)
        if matched_id:
            rows.append((matched_id, fields[6:]))

print(f"  Samples found in SNP matrix : {len(rows)}", flush=True)
missing = keep_ids - set(r[0] for r in rows)
if missing:
    print(f"  WARNING: {len(missing)} HDA IDs not found in SNP matrix", flush=True)
    print(f"  First 5 missing: {list(missing)[:5]}", flush=True)

# Transpose and write in chunks
sample_ids = [r[0] for r in rows]
n_samples = len(rows)
n_snps = len(snp_names)
CHUNK = 50000

with open(out_file, 'w') as out:
    out.write('snpid\t' + '\t'.join(sample_ids) + '\n')
    for chunk_start in range(0, n_snps, CHUNK):
        chunk_end = min(chunk_start + CHUNK, n_snps)
        for j in range(chunk_start, chunk_end):
            snp_row = [rows[i][1][j] for i in range(n_samples)]
            out.write(snp_names[j] + '\t' + '\t'.join(snp_row) + '\n')
        print(f"  Written SNPs {chunk_start}–{chunk_end} of {n_snps}", flush=True)

print(f"  Done: {n_snps} SNPs x {n_samples} samples -> {out_file}", flush=True)
EOF

echo "  SNP matrix done."

##############################################
# STEP 5: SNP location file from .bim
# .bim: chr, snpid, cm, pos, a1, a2
# Matrix eQTL needs: snpid, chr, pos
##############################################
echo ""
echo "[5] Generating SNP location file from .bim..."

SNP_LOC=$OUTDIR/snp_location.txt
echo -e "snpid\tchr\tpos" > $SNP_LOC
awk 'BEGIN{OFS="\t"} {print $2, "chr"$1, $4}' $BIM >> $SNP_LOC

N_SNPS=$(tail -n +2 $SNP_LOC | wc -l)
echo "  SNPs in location file : $N_SNPS"
echo "  Written to: $SNP_LOC"

##############################################
# STEP 6: Gene location file from gencode v49 BED6
# BED6: chr, start, end, gene_id, score, strand
# Matrix eQTL needs: geneid, chr, left, right
# Strip version suffix from gene IDs (ENSG00000123.7 -> ENSG00000123)
##############################################
echo ""
echo "[6] Generating gene location file from gencode v49 BED6..."

GENE_LOC=$OUTDIR/gene_location.txt
echo -e "geneid\tchr\tleft\tright" > $GENE_LOC
awk 'BEGIN{OFS="\t"} {
    gene=$4; sub(/\.[0-9]+$/, "", gene)
    print gene, $1, $2, $3
}' $GTF_BED6 >> $GENE_LOC

N_GENES=$(tail -n +2 $GENE_LOC | wc -l)
echo "  Genes in location file : $N_GENES"
echo "  Written to: $GENE_LOC"

# Sanity check: gene ID overlap with expression matrix
python3 << 'EOF'
expr_file = "/home/zw529/donglab/data/target_ALS/eQTL/cerebellum/expression_cerebellum.txt"
loc_file  = "/home/zw529/donglab/data/target_ALS/eQTL/cerebellum/gene_location.txt"

with open(loc_file) as f:
    f.readline()
    loc_genes = set(line.split('\t')[0] for line in f)

with open(expr_file) as f:
    f.readline()
    expr_genes = [line.split('\t')[0] for line in f]

found = sum(1 for g in expr_genes if g in loc_genes)
missing = [g for g in expr_genes if g not in loc_genes]
print(f"  Expression genes              : {len(expr_genes)}", flush=True)
print(f"  Found in gene location file   : {found}", flush=True)
print(f"  Missing from location file    : {len(missing)}", flush=True)
if missing:
    print(f"  First 5 missing gene IDs      : {missing[:5]}", flush=True)
EOF

##############################################
# STEP 7: Final summary
##############################################
echo ""
echo "============================================"
echo "  Prep complete. Output files in $OUTDIR:"
echo ""
echo "  snp_cerebellum.txt         -> args[1] SNP_file_name"
echo "  expression_cerebellum.txt  -> args[2] expression_file_name"
echo "  covariates_cerebellum.txt  -> args[3] covariates_file_name"
echo "  gene_location.txt          -> args[5] gene_location_file_name"
echo "  snp_location.txt           -> args[6] snp_location_file_name"
echo ""
echo "  Next step: sbatch run_eQTL_cerebellum.sh"
echo "  $(date)"
echo "============================================"
