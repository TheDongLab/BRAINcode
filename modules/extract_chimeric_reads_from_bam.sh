#!/bin/bash
#SBATCH --job-name=remap_unmapped_chimeric
#SBATCH --array=1-5%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

###########################################
# Reference paths
###########################################
GENOME_BASE="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38"
STAR_GENOME="${GENOME_BASE}/Sequence/STAR"
GENOME_FA="${GENOME_BASE}/Sequence/WholeGenomeFasta/genome.fa"
GENOME_FAI="${GENOME_BASE}/Sequence/WholeGenomeFasta/genome.fa.fai"
GTF="${GENOME_BASE}/Annotation/gencode/gencode.v49.annotation.gtf"
REFFLAT="${GENOME_BASE}/Annotation/gencode/refFlat_cIRS7.txt"
ANNOT_BEDS="${STAR_GENOME}"

###########################################
# Environment
###########################################
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate RNAseq

CPU=${SLURM_CPUS_PER_TASK:-4}

###########################################
# Resolve sample for this array task
###########################################
ALS_ROOT="/home/zw529/donglab/data/target_ALS"

mapfile -t SAMPLE_DIRS < <(find "${ALS_ROOT}" \
    -mindepth 4 -maxdepth 4 \
    -type d \
    -path "*/RNAseq/Processed/*" | sort)

SAMPLE_DIR="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"

if [[ -z "${SAMPLE_DIR}" ]]; then
    echo "No sample for array index ${SLURM_ARRAY_TASK_ID}. Exiting."
    exit 0
fi

samplename=$(basename "${SAMPLE_DIR}")

exec > "${SAMPLE_DIR}/output_remap.out" 2> "${SAMPLE_DIR}/output_remap.err"

echo "[$(date)] Sample     : ${samplename}"
echo "[$(date)] Sample dir : ${SAMPLE_DIR}"

###########################################
# Locate the original sorted BAM — never touch it
# Exclude any *remap* files in case this script is re-run
###########################################
ORIG_BAM=$(find "${SAMPLE_DIR}" -maxdepth 1 -name "*.sorted.bam" ! -name "*remap*" | head -n 1)

if [[ -z "${ORIG_BAM}" ]]; then
    echo "[ERROR] No *.sorted.bam found in ${SAMPLE_DIR}. Exiting." >&2
    exit 1
fi

echo "[$(date)] Input BAM  : ${ORIG_BAM}  (will not be modified)"

# All remap outputs live in a dedicated subdirectory — fully isolated
# from the original BAM and its index
REMAP_DIR="${SAMPLE_DIR}/remap_chimeric"
mkdir -p "${REMAP_DIR}"

###########################################
# STEP 1: Extract unmapped reads
###########################################
UNMAPPED_BAM="${REMAP_DIR}/${samplename}.unmapped.bam"
UNMAPPED_FQ="${REMAP_DIR}/${samplename}.unmapped.fastq"

if [ ! -f "${REMAP_DIR}/.status.unmapped_extracted" ]; then
    echo "[STEP 1] Extracting unmapped reads with samtools view -f 4 ..."

    samtools view -@ "${CPU}" -f 4 -b "${ORIG_BAM}" -o "${UNMAPPED_BAM}"
    samtools fastq -@ "${CPU}" "${UNMAPPED_BAM}" > "${UNMAPPED_FQ}"

    touch "${REMAP_DIR}/.status.unmapped_extracted"
    echo "[STEP 1] Unmapped reads written to: ${UNMAPPED_FQ}"
else
    echo "[STEP 1] Unmapped reads already extracted, skipping."
fi

###########################################
# STEP 2: STAR remap (chimeric)
#
# Prefix: remap_chimeric/<samplename>.remap.
# Key outputs:
#   <samplename>.remap.Aligned.sortedByCoord.out.bam   <- new BAM, original untouched
#   <samplename>.remap.Chimeric.out.junction            <- input for CIRCexplorer2
###########################################
STAR_PREFIX="${REMAP_DIR}/${samplename}.remap."
REMAP_BAM="${STAR_PREFIX}Aligned.sortedByCoord.out.bam"
CHIM_JXN="${STAR_PREFIX}Chimeric.out.junction"

if [ ! -f "${REMAP_DIR}/.status.star_chimeric" ]; then
    echo "[STEP 2] STAR chimeric remapping starting..."

    STAR_TMP_DIR="${REMAP_DIR}/STARtmp_${SLURM_JOB_ID:-NA}_${SLURM_ARRAY_TASK_ID:-NA}"

    STAR --runThreadN "${CPU}" \
        --genomeDir "${STAR_GENOME}" \
        --readFilesIn "${UNMAPPED_FQ}" \
        --outFileNamePrefix "${STAR_PREFIX}" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattrRGline ID:${samplename} SM:${samplename} LB:${samplename} PL:ILLUMINA PU:${samplename} \
        --outFilterMultimapNmax 10 \
        --alignEndsType Local \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10 \
        --chimOutType Junctions \
        --outSAMstrandField intronMotif \
        --quantMode GeneCounts \
        --outTmpDir "${STAR_TMP_DIR}"

    touch "${REMAP_DIR}/.status.star_chimeric"
    rm -rf "${STAR_TMP_DIR}"
    echo "[STEP 2] STAR chimeric remapping completed."
else
    echo "[STEP 2] STAR remap already complete, skipping."
fi

###########################################
# STEP 3: Symlink key outputs to sample root
# for easy access (IGV, CIRCexplorer2, etc.)
# Symlinks are clearly named *remap* so they
# can never be confused with the original BAM.
###########################################
REMAP_BAM_LINK="${SAMPLE_DIR}/${samplename}.remap.sorted.bam"
if [[ ! -L "${REMAP_BAM_LINK}" ]]; then
    ln -s "${REMAP_BAM}" "${REMAP_BAM_LINK}"
    echo "[STEP 3] Symlink: ${REMAP_BAM_LINK} -> ${REMAP_BAM}"
fi

# Named STAR.Chimeric.out.junction so CIRCexplorer2 parse can be called
# from SAMPLE_DIR without specifying the full path into remap_chimeric/
CHIM_JXN_LINK="${SAMPLE_DIR}/STAR.Chimeric.out.junction"
if [[ ! -L "${CHIM_JXN_LINK}" ]]; then
    ln -s "${CHIM_JXN}" "${CHIM_JXN_LINK}"
    echo "[STEP 3] Symlink: ${CHIM_JXN_LINK} -> ${CHIM_JXN}"
fi

###########################################
# STEP 4: circRNA calling with CIRCexplorer2
###########################################
if [ ! -f "${SAMPLE_DIR}/.status.RNAseq.circRNA" ]; then
    echo "[STEP 4] circRNA calling starting..."

    # Check file exists (may be 0 lines if no chimeric reads — still valid to proceed)
    NJXN=0
    if [[ -f "${CHIM_JXN}" ]]; then
        NJXN=$(wc -l < "${CHIM_JXN}")
    fi
    echo "[STEP 4] Chimeric.out.junction: ${NJXN} lines (${CHIM_JXN})"

    if [[ ! -f "${CHIM_JXN}" ]]; then
        echo "[STEP 4] WARNING: Chimeric.out.junction not found — skipping circRNA calling." >&2
        touch "${SAMPLE_DIR}/.status.RNAseq.circRNA"
    else
        cd "${SAMPLE_DIR}" || { echo "[STEP 4] ERROR: cannot cd to ${SAMPLE_DIR}" >&2; exit 1; }

        echo "[STEP 4] Running CIRCexplorer2 parse..."
        CIRCexplorer2 parse \
            -t STAR \
            STAR.Chimeric.out.junction \
            > back_spliced_junction.txt
        if [[ $? -ne 0 ]]; then
            echo "[STEP 4] ERROR: CIRCexplorer2 parse failed." >&2; exit 1
        fi

        echo "[STEP 4] Running CIRCexplorer2 annotate..."
        CIRCexplorer2 annotate \
            -r "${REFFLAT}" \
            -g "${GENOME_FA}" \
            -b back_spliced_junction.bed \
            -o circularRNA_known.txt \
            --low-confidence
        if [[ $? -ne 0 ]]; then
            echo "[STEP 4] ERROR: CIRCexplorer2 annotate failed." >&2; exit 1
        fi

        if [ -s circularRNA_known.txt ]; then

            [ -f "${REMAP_BAM}.bai" ] || samtools index "${REMAP_BAM}"

            python3 ~/donglab/pipelines/scripts/rnaseq/circ_percent_calculation.py \
                "${SAMPLE_DIR}/circularRNA_known.txt"

            HEADER="chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\texonCount\texonSizes\texonOffsets\treadNumber\tcircType\tgeneName\tisoformName\tindex\tflankIntron"

            { echo -e "${HEADER}"; cat "${SAMPLE_DIR}/circularRNA_known.txt"; } \
                > "${SAMPLE_DIR}/circularRNA_known.tmp" \
                && mv "${SAMPLE_DIR}/circularRNA_known.tmp" "${SAMPLE_DIR}/circularRNA_known.txt"

            { echo -e "${HEADER}"; cat "${SAMPLE_DIR}/low_conf_circularRNA_known.txt"; } \
                > "${SAMPLE_DIR}/low_conf_circularRNA_known.tmp" \
                && mv "${SAMPLE_DIR}/low_conf_circularRNA_known.tmp" "${SAMPLE_DIR}/low_conf_circularRNA_known.txt"

            touch "${SAMPLE_DIR}/.status.RNAseq.circRNA"
            echo "[STEP 4] SUCCESS: Found $(wc -l < "${SAMPLE_DIR}/circularRNA_known.txt") circRNAs."
            echo "[STEP 4] circRNA calling completed successfully."

        else
            echo "[STEP 4] ERROR: circularRNA_known.txt is empty — no circRNAs annotated." >&2
            exit 1
        fi
    fi
else
    echo "[STEP 4] circRNA calling already complete, skipping."
fi

###########################################
# STEP 5: Cleanup intermediate files
# Runs after circRNA calling so nothing is
# removed prematurely on re-runs
###########################################
echo "[STEP 5] Cleaning up intermediate files..."

if [[ -f "${UNMAPPED_FQ}" ]]; then
    rm -f "${UNMAPPED_FQ}"
    echo "[STEP 5] Removed: ${UNMAPPED_FQ}"
fi

echo "[$(date)] Done: ${samplename}"
