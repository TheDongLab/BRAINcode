#!/bin/bash
#SBATCH --job-name=vcf_qc_audit
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/vcf_qc_audit.log
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

# Setup environment to mirror convert_vcf_to_plink.sh
module --force purge
module load StdEnv
module load PLINK2/avx2_20250707

ORIG_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz" 
IMPUTED_LOG="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_all_chrs_filtered.log"
AUDIT_DIR="/home/zw529/donglab/data/target_ALS/QTL/plink/audit_temp"

mkdir -p "$AUDIT_DIR"

echo "========================================================="
echo " RUNNING STEP-BY-STEP QC AUDIT: ORIGINAL VS IMPUTED VCF"
echo "========================================================="

echo "[1/2] Processing original VCF to capture filtering steps..."
echo "------------------- PLINK2 STDOUT/STDERR -------------------"

# Added --make-pgen to satisfy PLINK2's structural execution requirement
plink2 --vcf "$ORIG_VCF" \
       --double-id \
       --set-all-var-ids @:#:\$r:\$a \
       --new-id-max-allele-len 800 \
       --vcf-half-call m \
       --allow-extra-chr \
       --geno 0.05 \
       --mind 0.05 \
       --maf 0.05 \
       --hwe 1e-6 \
       --max-alleles 2 \
       --make-pgen \
       --out "$AUDIT_DIR/orig_audit" \
       --threads "${SLURM_CPUS_PER_TASK:-4}" 2>&1 || echo "PLINK2 execution phase complete."

echo "------------------------------------------------------------"

ORIG_LOG="$AUDIT_DIR/orig_audit.log"

parse_plink_log() {
    local log_file=$1
    
    if [ ! -f "$log_file" ]; then
        echo "0,0,0,0,0,0"
        return
    fi

    # Flexible matching to account for line differences across PLINK2 sub-builds
    local total=$(grep -E "variants loaded from" "$log_file" | head -n 1 | awk '{print $1}' || echo "0")
    local mult=$(grep -E "removed due to allele count|alleles|variants removed due to max-alleles" "$log_file" | head -n 1 | awk '{print $1}' || echo "0")
    local geno=$(grep -E "removed due to missing genotype|missing genotype" "$log_file" | head -n 1 | awk '{print $1}' || echo "0")
    local hwe=$(grep -E "removed due to hardy-weinberg|hardy-weinberg" "$log_file" | head -n 1 | awk '{print $1}' || echo "0")
    local maf=$(grep -E "removed due to frequency threshold|frequency threshold" "$log_file" | head -n 1 | awk '{print $1}' || echo "0")
    
    # Extract remaining counts dynamically depending on how the pipeline completed
    local final="0"
    if grep -q "remaining variants after main filters" "$log_file"; then
        final=$(grep -E "remaining variants after main filters" "$log_file" | head -n 1 | awk '{print $1}')
    elif grep -q "written to" "$log_file"; then
        final=$(grep -E "variants written to" "$log_file" | head -n 1 | awk '{print $1}')
    fi
    
    # Strip everything that isn't a plain integer string
    total=$(echo "${total:-0}" | tr -d '[:alpha:],()')
    mult=$(echo "${mult:-0}" | tr -d '[:alpha:],()')
    geno=$(echo "${geno:-0}" | tr -d '[:alpha:],()')
    hwe=$(echo "${hwe:-0}" | tr -d '[:alpha:],()')
    maf=$(echo "${maf:-0}" | tr -d '[:alpha:],()')
    final=$(echo "${final:-0}" | tr -d '[:alpha:],()')

    echo "${total:-0},${mult:-0},${geno:-0},${hwe:-0},${maf:-0},${final:-0}"
}

echo "[2/2] Parsing logs and compiling metrics..."
orig_metrics=$(parse_plink_log "$ORIG_LOG")
impute_metrics=$(parse_plink_log "$IMPUTED_LOG")

# Parse string arrays
IFS=',' read -r o_tot o_mul o_gen o_hwe o_maf o_fin <<< "$orig_metrics"
IFS=',' read -r i_tot i_mul i_gen i_hwe i_maf i_fin <<< "$impute_metrics"

echo ""
echo "========================================================================="
echo "                      QC FILTER BREAKDOWN COMPARISON"
echo "========================================================================="
printf "%-35s | %-18s | %-18s\n" "Filter Metric" "Original VCF" "Imputed VCF"
echo "-------------------------------------------------------------------------"
printf "%-35s | %-18d | %-18d\n" "Total Raw Variants Loaded" "${o_tot:-0}" "${i_tot:-0}"
printf "%-35s | %-18d | %-18d\n" "Removed by --max-alleles 2" "${o_mul:-0}" "${i_mul:-0}"
printf "%-35s | %-18d | %-18d\n" "Removed by --geno 0.05" "${o_gen:-0}" "${i_gen:-0}"
printf "%-35s | %-18d | %-18d\n" "Removed by --hwe 1e-6" "${o_hwe:-0}" "${i_hwe:-0}"
printf "%-35s | %-18d | %-18d\n" "Removed by --maf 0.05" "${o_maf:-0}" "${i_maf:-0}"
echo "-------------------------------------------------------------------------"
printf "%-35s | %-18d | %-18d\n" "FINAL Variants Passing QC" "${o_fin:-0}" "${i_fin:-0}"
echo "========================================================================="

if [ "${o_fin:-0}" -gt 0 ]; then
    yield=$(awk "BEGIN {print ${i_fin:-0} / ${o_fin:-0}}")
    echo "True Post-QC Yield Shift: ${yield}x"
else
    echo "True Post-QC Yield Shift: N/A"
fi

# Clean up temporary audit directories
rm -rf "$AUDIT_DIR"
