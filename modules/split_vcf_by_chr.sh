#!/bin/bash
#SBATCH --job-name=split_vcf_QC_filtered
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/split_vcf.log
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

#########################################################################
# This script:
# 1. Extracts the "Pass QC" samples from the PLINK2 filtered dataset.
# 2. Creates a clean male list for those samples.
# 3. Splits the joint VCF by chromosome (1-22, X, Y).
# 4. Forces haploid genotypes for males on ChrX to satisfy TOPMed/eQTL requirements.
#########################################################################

module load BCFtools

# Define paths
INPUT_VCF="/home/zw529/donglab/data/target_ALS/QTL/joint_genotyped_GQ.vcf.gz"
OUTPUT_DIR="/home/zw529/donglab/data/target_ALS/QTL/chromosome_joint_vcfs"
METADATA="/home/zw529/donglab/data/target_ALS/targetALS_rnaseq_metadata.csv"
PSAM_FILE="/home/zw529/donglab/data/target_ALS/QTL/plink/joint_autosomes_filtered.psam"

mkdir -p $OUTPUT_DIR

# --- STEP 1: CREATE SAMPLE LISTS ---
PASS_SAMPLES="${OUTPUT_DIR}/pass_samples.txt"
awk 'NR>1 {print $1}' $PSAM_FILE > $PASS_SAMPLES

FEMALE_LIST="${OUTPUT_DIR}/females.txt"
MALE_LIST="${OUTPUT_DIR}/males.txt"

# Clean metadata and create sex map
tr -d '\r' < $METADATA | awk -F',' 'NR>1 && $2 != "" {print $2, tolower($5)}' | sort -u > ${OUTPUT_DIR}/subject_sex_map.txt

grep -Fwf $PASS_SAMPLES ${OUTPUT_DIR}/subject_sex_map.txt | awk '$2=="female" {print $1}' > $FEMALE_LIST
grep -Fwf $PASS_SAMPLES ${OUTPUT_DIR}/subject_sex_map.txt | awk '$2=="male" {print $1}' > $MALE_LIST

echo "Females: $(wc -l < $FEMALE_LIST)"
echo "Males: $(wc -l < $MALE_LIST)"

# --- STEP 2: PROCESS CHRX ---
echo "Processing chrX..."

# 2a. Extract females (keep as-is, diploid)
bcftools view -r chrX -S $FEMALE_LIST --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/tmp_X_females.vcf.gz
bcftools index -t ${OUTPUT_DIR}/tmp_X_females.vcf.gz

# 2b. Extract males as diploid (uncompressed intermediate)
bcftools view -r chrX -S $MALE_LIST --force-samples $INPUT_VCF -Ou -o ${OUTPUT_DIR}/tmp_X_males_diploid.vcf

# 2c. Convert males to haploid using bcftools +setGT
# -t q -i 'GT!="."': apply to all non-missing genotypes
# -n c:1: convert to ploidy 1 (haploid)
bcftools +setGT -Oz -o ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz ${OUTPUT_DIR}/tmp_X_males_diploid.vcf -- -t q -i 'GT!="."' -n c:1

bcftools index -t ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz

# Verify the conversion worked
echo "Female ploidy check (should be diploid like 0/0, 0/1, 1/1):"
bcftools query -f '[%GT\n]' ${OUTPUT_DIR}/tmp_X_females.vcf.gz | head -5

echo "Male ploidy check (should be haploid like 0, 1):"
bcftools query -f '[%GT\n]' ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz | head -5

# 2d. Merge females and haploid males using embedded Python script
# This preserves per-sample ploidy (females diploid, males haploid)
python3 << "PYTHON_SCRIPT"
import sys
import gzip

def read_vcf_header(filename):
    """Read VCF header and return metadata lines and sample list."""
    header_lines = []
    samples = []
    
    with gzip.open(filename, 'rt') if filename.endswith('.gz') else open(filename, 'r') as f:
        for line in f:
            if line.startswith('##'):
                header_lines.append(line.rstrip())
            elif line.startswith('#CHROM'):
                fields = line.rstrip().split('\t')
                samples = fields[9:]  # Sample names start at column 9
                header_lines.append(line.rstrip())
                break
    
    return header_lines, samples

def read_vcf_variants(filename):
    """Read VCF variants (non-header lines) and return list of tuples."""
    variants = []
    
    with gzip.open(filename, 'rt') if filename.endswith('.gz') else open(filename, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                variants.append(line.rstrip().split('\t'))
    
    return variants

def merge_vcfs(female_vcf, male_vcf, output_vcf):
    """Merge female and male VCFs, preserving per-sample ploidy."""
    
    # Read female VCF
    female_header, female_samples = read_vcf_header(female_vcf)
    female_variants = read_vcf_variants(female_vcf)
    
    # Read male VCF
    male_header, male_samples = read_vcf_header(male_vcf)
    male_variants = read_vcf_variants(male_vcf)
    
    print(f"Female samples: {len(female_samples)}")
    print(f"Male samples: {len(male_samples)}")
    print(f"Female variants: {len(female_variants)}")
    print(f"Male variants: {len(male_variants)}")
    
    # Create combined sample list (females first, then males)
    combined_samples = female_samples + male_samples
    
    # Index variants by position for quick lookup
    female_by_pos = {(v[0], v[1], v[3], v[4]): v for v in female_variants}
    male_by_pos = {(v[0], v[1], v[3], v[4]): v for v in male_variants}
    
    # Get all unique positions
    all_positions = set(female_by_pos.keys()) | set(male_by_pos.keys())
    
    # Sort by chromosome and position (assuming chrX)
    all_positions = sorted(all_positions, key=lambda x: int(x[1]))
    
    print(f"Total unique positions: {len(all_positions)}")
    
    # Write output VCF
    with gzip.open(output_vcf, 'wt') as out:
        # Write header
        for line in female_header[:-1]:  # Exclude old #CHROM line
            out.write(line + '\n')
        
        # Write new #CHROM line with combined samples
        chrom_line = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(combined_samples)
        out.write(chrom_line + '\n')
        
        # Write merged variants
        for pos in all_positions:
            female_var = female_by_pos.get(pos)
            male_var = male_by_pos.get(pos)
            
            # Use female variant as template (they should be identical for shared positions)
            if female_var:
                template = female_var
            else:
                template = male_var
            
            # Build genotype columns
            genotypes = []
            
            # Female genotypes
            if female_var:
                # Extract genotypes from columns 9 onward
                for gt in female_var[9:]:
                    genotypes.append(gt)
            else:
                # Missing for this variant
                for _ in female_samples:
                    genotypes.append('./.')
            
            # Male genotypes
            if male_var:
                # Extract genotypes from columns 9 onward
                for gt in male_var[9:]:
                    genotypes.append(gt)
            else:
                # Missing for this variant
                for _ in male_samples:
                    genotypes.append('.')
            
            # Build output line
            out_line = '\t'.join(template[0:9]) + '\t' + '\t'.join(genotypes)
            out.write(out_line + '\n')
    
    print(f"Merged VCF written to: {output_vcf}")

female_vcf = sys.argv[1]
male_vcf = sys.argv[2]
output_vcf = sys.argv[3]

merge_vcfs(female_vcf, male_vcf, output_vcf)
PYTHON_SCRIPT "${OUTPUT_DIR}/tmp_X_females.vcf.gz" "${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz" "${OUTPUT_DIR}/target_ALS_chrX.vcf.gz"
PYTHON_SCRIPT

bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz

# --- STEP 3: ALL OTHER CHROMOSOMES ---
for chr in {1..22} Y; do
    echo "Processing chr${chr}..."
    bcftools view -r chr${chr} -S $PASS_SAMPLES --force-samples $INPUT_VCF -Oz -o ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
    bcftools index -f -t ${OUTPUT_DIR}/target_ALS_chr${chr}.vcf.gz
done

# --- STEP 4: CLEANUP ---
rm -f ${OUTPUT_DIR}/tmp_X_females.vcf.gz*
rm -f ${OUTPUT_DIR}/tmp_X_males_diploid.vcf*
rm -f ${OUTPUT_DIR}/tmp_X_males_hap.vcf.gz*

echo "Successfully created mixed-ploidy chrX VCF!"
echo "chrX output: ${OUTPUT_DIR}/target_ALS_chrX.vcf.gz"
