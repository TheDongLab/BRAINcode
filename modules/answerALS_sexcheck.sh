#!/usr/bin/env bash
#SBATCH --job-name=answerALS_chrY
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/iPSC_Motor_Neuron/eQTL/chrY_ploidy_%j.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/iPSC_Motor_Neuron/eQTL/chrY_ploidy_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

module load BCFtools/1.21-GCC-13.3.0
module load Python/3.12.3-GCCcore-13.3.0

BASE="/home/zw529/donglab/data/answer_ALS"
EQTL="$BASE/iPSC_Motor_Neuron/eQTL"
VCF="$BASE/QTL/joint_genotyped_GQ.vcf.gz"
METADATA="$BASE/aals_dataportal_datatable.csv"
MISMATCH="$EQTL/potential_sex_mismatches.txt"
WORK="$EQTL/chrY_ploidy_check"
OUTPUT="$EQTL/sex_metadata_nonPAR_chrY_ploidy.tsv"

mkdir -p "$WORK"

for file in "$VCF" "$METADATA" "$MISMATCH"; do
    [[ -s "$file" ]] || {
        echo "[ERROR] Missing or empty input: $file" >&2
        exit 1
    }
done

if [[ ! -s "$VCF.tbi" && ! -s "$VCF.csi" ]]; then
    bcftools index -t "$VCF"
fi

# GRCh38 non-PAR chrY intervals in BED coordinates.
# These exclude chrY PAR1 (10,001-2,781,479) and PAR2
# (56,887,903-57,217,415).
REGIONS="$WORK/chrY_nonPAR_hg38.bed"
cat > "$REGIONS" <<'EOF'
chrY	0	10000
chrY	2781479	56887902
chrY	57217415	57227415
EOF

python3 - "$VCF" "$METADATA" "$MISMATCH" "$REGIONS" "$OUTPUT" <<'PY'
import csv
import subprocess
import sys

vcf, metadata_file, mismatch_file, regions_file, output_file = sys.argv[1:]

with open(metadata_file, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    id_col = reader.fieldnames[0]
    sex_col = next(
        (c for c in reader.fieldnames if c.strip().lower() == "sex"),
        None
    )
    if sex_col is None:
        raise SystemExit("No Sex column found in AnswerALS metadata.")

    metadata_sex = {}
    for row in reader:
        sample = (row.get(id_col) or "").strip()
        if sample:
            value = (row.get(sex_col) or "").strip().lower()
            metadata_sex[sample] = (
                "Male" if value == "male"
                else "Female" if value == "female"
                else "Unknown"
            )

with open(mismatch_file) as f:
    mismatch_ids = {
        line.strip() for line in f if line.strip()
    }

samples = subprocess.check_output(
    ["bcftools", "query", "-l", vcf],
    text=True
).splitlines()

counts = {
    sample: {"haploid": 0, "diploid": 0, "partial": 0, "nonref": 0}
    for sample in samples
}

query = subprocess.Popen(
    [
        "bcftools", "query",
        "-R", regions_file,
        "-f", "%POS[\\t%GT]\\n",
        vcf
    ],
    stdout=subprocess.PIPE,
    text=True
)

for line in query.stdout:
    fields = line.rstrip("\n").split("\t")
    genotypes = fields[1:]

    for sample, gt in zip(samples, genotypes):
        alleles = gt.replace("|", "/").split("/")

        if all(a == "." for a in alleles):
            continue
        if any(a == "." for a in alleles):
            counts[sample]["partial"] += 1
            continue

        if len(alleles) == 1:
            counts[sample]["haploid"] += 1
        elif len(alleles) == 2:
            counts[sample]["diploid"] += 1

        if any(a != "0" for a in alleles):
            counts[sample]["nonref"] += 1

if query.wait() != 0:
    raise SystemExit("bcftools query failed.")

with open(output_file, "w") as out:
    out.write(
        "sample_id\tmetadata_sex\tnonPAR_chrY_haploid_GT_sites\t"
        "nonPAR_chrY_diploid_GT_sites\tpartial_GT_sites\t"
        "nonPAR_chrY_nonref_GT_sites\tploidy_check\t"
        "in_potential_sex_mismatches\n"
    )

    for sample in samples:
        sex = metadata_sex.get(sample, "Unknown")
        c = counts[sample]

        if sex == "Male":
            if c["haploid"] > 0 and c["diploid"] == 0:
                check = "HAPLOID_CALLS_AS_EXPECTED"
            elif c["diploid"] > 0:
                check = "DIPLOID_CALLS_REVIEW"
            else:
                check = "NO_NONPAR_Y_CALLS_REVIEW"
        elif sex == "Female":
            if c["haploid"] == 0 and c["diploid"] == 0:
                check = "NO_NONPAR_Y_CALLS_AS_EXPECTED"
            else:
                check = "NONPAR_Y_CALLS_REVIEW"
        else:
            check = "UNKNOWN_METADATA_SEX"

        out.write(
            f"{sample}\t{sex}\t{c['haploid']}\t{c['diploid']}\t"
            f"{c['partial']}\t{c['nonref']}\t{check}\t"
            f"{'YES' if sample in mismatch_ids else 'NO'}\n"
        )

print(f"Samples checked: {len(samples)}")
print(f"Mismatch-list IDs loaded: {len(mismatch_ids)}")
print(f"Results: {output_file}")
PY
