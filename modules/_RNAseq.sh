#!/bin/bash
###########################################
# RNAseq pipeline (paired-end, stranded)
# Author: Xianjun Dong & Zachery Wolfe (Zachery updated)
# Date: 3/18/2026
# Version: 5.2 (Added preliminary step to handle existing .bam output. Changed BAM references in further steps so it is not automatically hard-coded to STAR .bam)
###########################################

set -euo pipefail
 
export xml_catalog_files_libxml2="${xml_catalog_files_libxml2:-}"
export CONDA_BACKUP_CXX="${CONDA_BACKUP_CXX:-}"
export CONDA_BACKUP_GXX="${CONDA_BACKUP_GXX:-}"
 
###########################################
# Conda activation
###########################################
set +u
source ~/donglab/pipelines/modules/miniconda3/etc/profile.d/conda.sh
conda activate RNAseq
set -u
 
###########################################
# Clean PATH & Kent_tools
###########################################
export PATH="$HOME/.conda/envs/RNAseq/bin:/apps/software/2022b/software/Kent_tools/468-GCC-13.3.0/bin:/opt/slurm/current/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin"
 
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
chmod +x "$SCRIPT_DIR/bam2bigwig.sh" "$SCRIPT_DIR/_bam2annotation.r" "$SCRIPT_DIR/normalize_counts.py"
 
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
# Inputs
###########################################
R1="$1"
R2="$2"
outputdir="$3"
 
mkdir -p "$outputdir"
 
CPU=${CPU:-4}
samplename=$(basename "$outputdir")
SAMPLE_DIR="$outputdir"
 
###########################################
# Dynamic BAM resolution
# Use STAR output if present, otherwise find existing BAM in sample dir
###########################################
STAR_BAM="$SAMPLE_DIR/STAR.Aligned.sortedByCoord.out.bam"
if [ -f "$STAR_BAM" ]; then
    BAM="$STAR_BAM"
else
    BAM=$(ls "$SAMPLE_DIR"/*.sorted.bam "$SAMPLE_DIR"/*.final.bam 2>/dev/null | head -1 || true)
    if [ -z "$BAM" ]; then
        echo "[INFO] No BAM file found yet — will be created by STAR in step 4."
        BAM="$STAR_BAM"
    else
        echo "[INFO] Using existing BAM: $(basename "$BAM")"
        # Symlink to expected name so all downstream steps reference consistently
        ln -sf "$BAM" "$STAR_BAM"
        [ -f "${BAM}.bai" ] && ln -sf "${BAM}.bai" "${STAR_BAM}.bai"
        BAM="$STAR_BAM"
    fi
fi
 
###########################################
# STEP 1: Quality check (FastQC)
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.fastqc" ]; then
    echo "[STEP 1] FastQC quality check starting..."
    fastqc -t $CPU --nogroup -o "$SAMPLE_DIR" "$R1" "$R2" && \
    touch "$SAMPLE_DIR/.status.RNAseq.fastqc" && \
    echo "[STEP 1] FastQC completed successfully."
fi
 
###########################################
# STEP 2: Trimming (Trim Galore)
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.trim" ]; then
    echo "[STEP 2] Read trimming (Trim Galore) starting..."
    trim_galore --paired --cores $CPU --output_dir "$SAMPLE_DIR" "$R1" "$R2" && \
    touch "$SAMPLE_DIR/.status.RNAseq.trim" && \
    echo "[STEP 2] Trimming completed successfully."
fi
 
R1="$SAMPLE_DIR/$(basename "$R1" .fastq.gz)_val_1.fq.gz"
R2="$SAMPLE_DIR/$(basename "$R2" .fastq.gz)_val_2.fq.gz"
 
###########################################
# STEP 3: kPAL contamination check
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.kpal" ]; then
    echo "[STEP 3] kPAL contamination check starting..."
    KPAL_DIR="$SAMPLE_DIR/kpal_results"
    KPAL_OUT="$KPAL_DIR/counts.k9.txt"
    mkdir -p "$KPAL_DIR" && \
    gunzip -c "$R1" > "$KPAL_DIR/temp_R1.fastq" && \
    gunzip -c "$R2" > "$KPAL_DIR/temp_R2.fastq" && \
    "$HOME/.conda/envs/RNAseq/bin/kpal" count \
        "$KPAL_DIR/temp_R1.fastq" \
        "$KPAL_DIR/temp_R2.fastq" \
        "$KPAL_OUT" && \
    rm -f "$KPAL_DIR"/temp_R{1,2}.fastq && \
    touch "$SAMPLE_DIR/.status.RNAseq.kpal" && \
    echo "[STEP 3] kPAL contamination check completed successfully."
fi
 
###########################################
# STEP 4: STAR mapping
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.mapping" ]; then
    echo "[STEP 4] STAR mapping starting..."
    STAR_TMP_DIR="$SAMPLE_DIR/STARtmp_${SLURM_JOB_ID:-NA}_${SLURM_ARRAY_TASK_ID:-NA}"
    STAR --runThreadN $CPU \
        --genomeDir "$STAR_GENOME" \
        --readFilesIn "$R1" "$R2" --readFilesCommand zcat \
        --outFileNamePrefix "$SAMPLE_DIR/STAR." \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattrRGline ID:$samplename SM:$samplename LB:$samplename PL:ILLUMINA PU:$samplename \
        --outFilterMultimapNmax 10 \
        --alignEndsType Local \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10 \
        --chimOutType Junctions \
        --outSAMstrandField intronMotif \
        --quantMode GeneCounts \
        --outTmpDir "$STAR_TMP_DIR" && \
    touch "$SAMPLE_DIR/.status.RNAseq.mapping" && \
    rm -rf "$STAR_TMP_DIR" && \
    # Update BAM to the newly created STAR output
    BAM="$STAR_BAM" && \
    echo "[STEP 4] STAR mapping completed successfully."
fi
 
###########################################
# STEP 5: circRNA calling
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.circRNA" ]; then
    echo "[STEP 5] circRNA calling starting..."
 
    # Skip if no chimeric junction file (e.g. pre-processed BAM-only samples)
    if [ ! -f "$SAMPLE_DIR/STAR.Chimeric.out.junction" ]; then
        echo "[STEP 5] WARNING: No STAR.Chimeric.out.junction found — skipping circRNA calling."
        touch "$SAMPLE_DIR/.status.RNAseq.circRNA"
    else
        cd "$SAMPLE_DIR" && \
 
        # Parse STAR chimeric junctions (Local alignment; default CIRCexplorer2 behavior)
        CIRCexplorer2 parse \
            -t STAR \
            STAR.Chimeric.out.junction \
            > back_spliced_junction.txt && \
 
        # Annotate with default autofix
        CIRCexplorer2 annotate \
            -r "$REFFLAT" \
            -g "$GENOME_FA" \
            -b back_spliced_junction.bed \
            -o circularRNA_known.txt \
            --low-confidence && \
 
        if [ -s circularRNA_known.txt ]; then
 
            # Ensure BAM index exists (required for circ percentage calculation)
            [ -f "${BAM}.bai" ] || samtools index "$BAM"
 
            # Circ percentage calculation
            python3 ~/donglab/pipelines/scripts/rnaseq/circ_percent_calculation.py \
                "$SAMPLE_DIR/circularRNA_known.txt" && \
 
            # Prepend official CIRCexplorer2 column headers to circularRNA_known.txt
            { echo -e "chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\texonCount\texonSizes\texonOffsets\treadNumber\tcircType\tgeneName\tisoformName\tindex\tflankIntron"; \
            cat "$SAMPLE_DIR/circularRNA_known.txt"; } > "$SAMPLE_DIR/circularRNA_known.tmp" && \
            mv "$SAMPLE_DIR/circularRNA_known.tmp" "$SAMPLE_DIR/circularRNA_known.txt" && \
 
            # Prepend official CIRCexplorer2 column headers to low_conf_circularRNA_known.txt
            { echo -e "chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\texonCount\texonSizes\texonOffsets\treadNumber\tcircType\tgeneName\tisoformName\tindex\tflankIntron"; \
            cat "$SAMPLE_DIR/low_conf_circularRNA_known.txt"; } > "$SAMPLE_DIR/low_conf_circularRNA_known.tmp" && \
            mv "$SAMPLE_DIR/low_conf_circularRNA_known.tmp" "$SAMPLE_DIR/low_conf_circularRNA_known.txt" && \
 
            touch "$SAMPLE_DIR/.status.RNAseq.circRNA" && \
            echo "[STEP 5] SUCCESS: Found $(wc -l < "$SAMPLE_DIR/circularRNA_known.txt") circRNAs." && \
            echo "[STEP 5] circRNA calling completed successfully."
 
        else
            echo "[STEP 5] ERROR: No circRNAs annotated." && exit 1
        fi
    fi
fi
 
###########################################
# STEP 6: BAM post-processing and annotation (summary table + BEDs + PDF)
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.bam2annotation" ]; then
    echo "[STEP 6] BAM post-processing and annotation starting..."
 
    SUMMARY="$SAMPLE_DIR/bam.summary.txt"
    echo -e "total_non_rRNA_mt\texons\tintrons\trRNA\tmtRNA\tLINE\tSINE\tintergenic\tintergenic_near_genes\tintergenic_not_near_genes\ttotal" > $SUMMARY
 
    count() { samtools view -c -F 0x100 -L "$1" "$BAM"; }
 
    total_non_rRNA_mt=$(count "${ANNOT_BEDS}/non_rRNA_mt.sorted.bed")
    exons=$(count "${ANNOT_BEDS}/exons.sorted.bed")
    introns=$(count "${ANNOT_BEDS}/introns.sorted.bed")
    rRNA=$(count "${ANNOT_BEDS}/rRNA.sorted.bed")
    mtRNA=$(count "${ANNOT_BEDS}/mtRNA.sorted.bed")
    LINE=$(count "${ANNOT_BEDS}/LINE.sorted.bed")
    SINE=$(count "${ANNOT_BEDS}/SINE.sorted.bed")
    intergenic=$(count "${ANNOT_BEDS}/intergenic.sorted.bed")
    intergenic_not_near_genes=$(count "${ANNOT_BEDS}/intergenic_notneargene.sorted.bed")
    intergenic_near_genes=$((intergenic - intergenic_not_near_genes))
    total=$(samtools view -c -F 0x100 "$BAM")
 
    echo -e "${total_non_rRNA_mt}\t${exons}\t${introns}\t${rRNA}\t${mtRNA}\t${LINE}\t${SINE}\t${intergenic}\t${intergenic_near_genes}\t${intergenic_not_near_genes}\t${total}" >> $SUMMARY
 
    samtools view -cF 0x100 "$BAM" > "$SAMPLE_DIR/bam.stat"
    samtools flagstat "$BAM" >> "$SAMPLE_DIR/bam.stat"
 
    Rscript "$SCRIPT_DIR/_bam2annotation.r" "$SAMPLE_DIR/bam.summary.txt" "$SAMPLE_DIR/bam2annotation.pdf"
 
    touch "$SAMPLE_DIR/.status.RNAseq.bam2annotation"
    echo "[STEP 6] BAM annotation completed successfully."
fi
 
###########################################
# STEP 7: Gene counting (HTSeq)
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.htseqcount" ]; then
    echo "[STEP 7] Gene counting (HTSeq) starting..."
 
    htseq-count -m intersection-strict -t exon -i gene_id -s yes -q -f bam -r pos \
        "$BAM" "$GTF" \
        > "$SAMPLE_DIR/htseqcount.tab" 2> "$SAMPLE_DIR/htseqcount.stderr" && \
 
    { echo -e "gene_id\tcount"; cat "$SAMPLE_DIR/htseqcount.tab"; } > "$SAMPLE_DIR/htseqcount.tab.tmp" && \
    mv "$SAMPLE_DIR/htseqcount.tab.tmp" "$SAMPLE_DIR/htseqcount.tab" && \
 
    touch "$SAMPLE_DIR/.status.RNAseq.htseqcount" && \
    echo "[STEP 7] Gene counting completed successfully."
fi
 
###########################################
# STEP 8: Read normalization (TPM, RPKM, FPKM)
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.normalization" ]; then
    echo "[STEP 8] Read normalization (TPM/RPKM/FPKM) starting..."
 
    python3 "$SCRIPT_DIR/normalize_counts.py" \
        "$SAMPLE_DIR/htseqcount.tab" \
        "$GTF" \
        "$SAMPLE_DIR" \
        --library-type PE && \
 
    touch "$SAMPLE_DIR/.status.RNAseq.normalization" && \
    echo "[STEP 8] Normalization completed successfully."
fi
 
###########################################
# STEP 9: LeafCutter junction quantification, clustering, and PSI
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.leafcutter" ]; then
    echo "[STEP 9] LeafCutter junction quantification starting..."
    LEAFCUTTER_BASE="/home/zw529/donglab/pipelines/modules/rnaseq/bin/leafcutter"
    LEAFCUTTER_OUT="$SAMPLE_DIR/leafcutter"
    JUNC_DIR="$LEAFCUTTER_OUT/juncs"
    CLUSTER_DIR="$LEAFCUTTER_OUT/clusters"
    PSI_DIR="$LEAFCUTTER_OUT/psi"
    export PATH="$LEAFCUTTER_BASE/scripts:$LEAFCUTTER_BASE/clustering:$PATH" && \
    mkdir -p "$JUNC_DIR" "$CLUSTER_DIR" "$PSI_DIR" && \
    "$LEAFCUTTER_BASE/scripts/bam2junc.sh" \
        "$BAM" \
        "$JUNC_DIR/${samplename}.junc" && \
    python "$LEAFCUTTER_BASE/clustering/leafcutter_quant_only.py" \
        -j "$JUNC_DIR"/*.junc \
        -o "$LEAFCUTTER_OUT/${samplename}" && \
    JUNC_LIST="$JUNC_DIR/juncfile_list.txt" && \
    ls "$JUNC_DIR"/*.junc > "$JUNC_LIST" && \
    python "$LEAFCUTTER_BASE/clustering/leafcutter_cluster.py" \
        -j "$JUNC_LIST" \
        -o leafcutter_clusters \
        -r "$CLUSTER_DIR" && \
    LC_PSI="$PSI_DIR/${samplename}.leafcutter.PSI.tsv" && \
    python "$LEAFCUTTER_BASE/scripts/leafcutter_psi.py" \
        "$CLUSTER_DIR/leafcutter_clusters_refined" \
        "$JUNC_DIR/${samplename}.junc" \
        "$LC_PSI" && \
    { head -n 1 "$LC_PSI" && tail -n +2 "$LC_PSI" | sort -k7,7nr; } > "$PSI_DIR/${samplename}.leafcutter.PSI.sorted.tsv" && \
    touch "$SAMPLE_DIR/.status.RNAseq.leafcutter" && \
    echo "[STEP 9] LeafCutter processing completed successfully."
fi
 
###########################################
# STEP 10: bigWig generation (+ / - split)
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.bigwig" ]; then
    echo "[STEP 10] bigWig generation starting..."
 
    # Ensure BAM index exists (required by genomecov)
    [ -f "${BAM}.bai" ] || samtools index "$BAM"
 
    "$SCRIPT_DIR/bam2bigwig.sh" \
        "$BAM" \
        -split && \
    touch "$SAMPLE_DIR/.status.RNAseq.bigwig" && \
    echo "[STEP 10] bigWig generation completed successfully."
fi
 
###########################################
# STEP 11: Cleanup
###########################################
if [ -f "$SAMPLE_DIR/.status.RNAseq.bigwig" ]; then
    echo "[STEP 11] Cleaning up trimmed FASTQ files, bedGraph files, and bed files..."
    rm -f "$SAMPLE_DIR"/*_val_1.fq.gz "$SAMPLE_DIR"/*_val_2.fq.gz && \
    rm -f "$SAMPLE_DIR"/*.bedGraph && \
    rm -f "$SAMPLE_DIR"/back_spliced* && \
    echo "[STEP 11] Trimmed FASTQs, bedGraphs, and circRNA back_spliced intermediate files removed."
fi
 
echo "[$(date)] RNAseq pipeline finished successfully and cleanup completed."
