#!/bin/bash
#SBATCH --job-name=recombine_imputed_vcfs
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/imputed_vcfs/recombine_imputed_vcfs.log
#SBATCH --mem=80G
#SBATCH --cpus-per-task=4
#SBATCH --time=20:00:00

export LC_ALL=C
module load BCFtools/1.21

WORKDIR="/home/zw529/donglab/data/target_ALS/QTL/imputed_vcfs"
OUTDIR="/home/zw529/donglab/data/target_ALS/QTL"
MALES="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs/males.txt"
ORIGINAL_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
IMPUTED_VCF="${OUTDIR}/target_ALS_imputed_filtered_joint.vcf.gz"
PASSWORD="4EDpN9tIAkPy@&6;db"   ### edit your password here

cd $WORKDIR

# --- STEP 0: DOWNLOAD FROM TOPMED SERVER ---
if [ -f "download_list.txt" ]; then
    echo "Starting download of imputed chromosome zip files..."
    while read -r url; do
        if [ ! -z "$url" ]; then
            echo "Downloading: $url"
            wget -q --show-progress "$url"
        fi
    done < download_list.txt
else
    echo "ERROR: download_list.txt not found in $WORKDIR. Skipping download step."
fi

# --- STEP 1: EXTRACTION ---
echo "Unzipping files..."
for zipfile in chr_*.zip; do
    unzip -P "$PASSWORD" -n "$zipfile"
done

# --- STEP 2: INDIVIDUAL CHROMOSOME FILTERING ---
echo "Filtering variants (R2 > 0.5 & MAF > 0.05)..."
mkdir -p filtered_tmp
for vcf in *.dose.vcf.gz; do
    out_name="filtered_tmp/filtered_${vcf}"
    echo "Processing $vcf..."
    bcftools filter -i 'R2>=0.5 && MAF>=0.05' -Oz -o "$out_name" "$vcf"
    bcftools index -t "$out_name"
done

# --- STEP 3: HAPLOID REVERSION FOR CHRX MALES ---
# We apply the fix to the filtered chrX file before concatenation
CHRX_FILTERED="filtered_tmp/filtered_chrX.dose.vcf.gz"
CHRX_FIXED="filtered_tmp/filtered_chrX_fixed.dose.vcf.gz"

if [ -f "$CHRX_FILTERED" ]; then
    echo "Applying PAR-aware haploid fix to males on chrX..."
    bcftools view "$CHRX_FILTERED" | \
    awk -v males_file="$MALES" -F'\t' 'BEGIN {
        OFS="\t";
        PAR1_START=10001; PAR1_END=2781479;
        PAR2_START=155701383; PAR2_END=156030895;
        while ((getline < males_file) > 0) {
            gsub(/[[:space:]]/, "", $1);
            if ($1 != "") male_list[$1] = 1;
        }
    }
    {
        if (/^#/) {
            if (/^#CHROM/) {
                for (i=10; i<=NF; i++) {
                    sample_name = $i;
                    gsub(/[[:space:]]/, "", sample_name);
                    if (sample_name in male_list) is_male[i] = 1;
                }
            }
            print $0;
        } else {
            POS = $2;
            is_par = ((POS >= PAR1_START && POS <= PAR1_END) || (POS >= PAR2_START && POS <= PAR2_END));
            for (i=10; i<=NF; i++) {
                if (is_male[i] && !is_par) {
                    split($i, a, ":");
                    split(a[1], gt, /[\/|]/);
                    a[1] = gt[1]; 
                    new_field = a[1];
                    for (j=2; j<=length(a); j++) new_field = new_field ":" a[j];
                    $i = new_field;
                }
            }
            print $0;
        }
    }' | bcftools view -Oz -o "$CHRX_FIXED"
    bcftools index -t "$CHRX_FIXED"
    mv "$CHRX_FIXED" "$CHRX_FILTERED"
    mv "${CHRX_FIXED}.tbi" "${CHRX_FILTERED}.tbi"
fi

# --- STEP 4: PREP CONCATENATION ---
# Sort numerically: 1..22 then X
ls filtered_tmp/filtered_*.dose.vcf.gz | sort -V > vcf_list.txt

# --- STEP 5: CONCATENATE ---
echo "Concatenating filtered VCFs into master file..."
bcftools concat -f vcf_list.txt -Oz -o "$IMPUTED_VCF"

# --- STEP 6: INDEXING ---
echo "Indexing final filtered master VCF..."
bcftools index -f -t "$IMPUTED_VCF"

# --- STEP 7: STATISTICAL COMPARISON ---
echo "-------------------------------------------------------"
echo "Generating Variant Count Comparison..."

# Count biallelic SNPs in Original (excluding MAF/Quality if desired, but here we just count total biallelic SNPs)
orig_count=$(bcftools view -m2 -M2 -v snps "$ORIGINAL_VCF" | grep -v "^#" | wc -l)

# Count biallelic SNPs in Imputed (Already filtered for R2 > 0.5 and MAF > 0.05 in Step 2)
impute_count=$(bcftools view -m2 -M2 -v snps "$IMPUTED_VCF" | grep -v "^#" | wc -l)

echo "Biallelic SNPs in Original VCF: $orig_count"
echo "Biallelic SNPs in Imputed VCF (R2>0.5, MAF>0.05): $impute_count"

if [ "$impute_count" -gt 0 ]; then
    fold_increase=$(awk "BEGIN {print $impute_count / $orig_count}")
    echo "Yield increase: ${fold_increase}x"
fi
echo "-------------------------------------------------------"

echo "Workflow complete."
echo "Master Filtered VCF: $IMPUTED_VCF"
