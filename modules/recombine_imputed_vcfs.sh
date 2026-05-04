#!/bin/bash
#SBATCH --job-name=recombine_imputed_vcfs
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/imputed_vcfs/recombine_imputed_vcfs.log
#SBATCH --mem=80G
#SBATCH --cpus-per-task=4
#SBATCH --time=20:00:00

###############################
# Run this script AFTER performing imputation on individual chromosome vcfs using the TOPMed Imputation Server
###############################

export LC_ALL=C
module load BCFtools/1.21

WORKDIR="/home/zw529/donglab/data/target_ALS/QTL/imputed_vcfs"
OUTDIR="/home/zw529/donglab/data/target_ALS/QTL"
PASSWORD="l_hpXa9mH*#J3Z4YtR"   ### edit your password here

cd $WORKDIR

# --- STEP 0: DOWNLOAD FROM TOPMED SERVER ---
# while read -r url; do
#      echo "Downloading: $url"
#      wget -q --show-progress "$url"
# done < download_list.txt

# --- STEP 1: EXTRACTION ---
echo "Unzipping files..."
for zipfile in chr_*.zip; do
    unzip -P "$PASSWORD" -n "$zipfile"
done

# --- STEP 2: INDIVIDUAL CHROMOSOME FILTERING ---
# Filtering here (R2 > 0.5 and MAF > 0.05) reduces the size of the final master file.
echo "Filtering variants (R2 > 0.5 & MAF > 0.05)..."
mkdir -p filtered_tmp
for vcf in *.dose.vcf.gz; do
    out_name="filtered_tmp/filtered_${vcf}"
    echo "Processing $vcf..."
    # INFO/R2 for Michigan/TOPMed; MAF calculated on the fly
    bcftools filter -i 'R2>=0.5 && MAF>=0.05' -Oz -o "$out_name" "$vcf"
    bcftools index -t "$out_name"
done

# --- STEP 3: PREP CONCATENATION ---
# Sort numerically: 1..22 then X (if present)
ls filtered_tmp/filtered_*.dose.vcf.gz | sort -V > vcf_list.txt

# --- STEP 4: CONCATENATE ---
echo "Concatenating filtered VCFs into master file..."
bcftools concat -f vcf_list.txt -na -Oz -o ${OUTDIR}/target_ALS_imputed_filtered_joint.vcf.gz

# --- STEP 5: INDEXING ---
echo "Indexing final filtered master VCF..."
bcftools index -f -t ${OUTDIR}/target_ALS_imputed_filtered_joint.vcf.gz

# --- STEP 6: CLEANUP (Optional) ---
# Uncomment the line below to delete unzipped raw VCFs and save space
# rm *.dose.vcf.gz filtered_tmp/filtered_*.dose.vcf.gz

echo "Workflow complete."
echo "Master Filtered VCF: ${OUTDIR}/target_ALS_imputed_filtered_joint.vcf.gz"
