#!/bin/bash
###########################################
# RNAseq pipeline (paired-end)
# Author: Xianjun Dong & Zachery Wolfe (Zachery updated)
# Date: 2/9/2026
# Version: 4.0 (removed cufflinks step, added bigWig generation step, added && to ensure that an error will stop the entire script)
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
chmod +x "$SCRIPT_DIR"/{bam2bigwig.sh,_bam2annotation.r}

###########################################
# Reference paths
###########################################
GENOME_BASE="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38"
STAR_GENOME="${GENOME_BASE}/Sequence/STAR"
GENOME_FA="${GENOME_BASE}/Sequence/WholeGenomeFasta/genome.fa"
GENOME_FAI="${GENOME_BASE}/Sequence/WholeGenomeFasta/genome.fa.fai"
GTF="${GENOME_BASE}/Annotation/gencode/gencode.v49.annotation.gtf"
REFFLAT="${GENOME_BASE}/Annotation/gencode/refFlat.txt"
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
# STEP 4: STAR mapping — modified for standard junctions
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.mapping" ]; then
    echo "[STEP 4] STAR mapping starting..."
    STAR_TMP_DIR="$SAMPLE_DIR/STARtmp_${SLURM_JOB_ID:-NA}_${SLURM_ARRAY_TASK_ID:-NA}"
    STAR --runThreadN $CPU \
        --genomeDir "$STAR_GENOME" \
        --readFilesIn "$R1" "$R2" --readFilesCommand zcat \
        --outFileNamePrefix "$SAMPLE_DIR/" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattrRGline ID:$samplename SM:$samplename LB:$samplename PL:ILLUMINA PU:$samplename \
        --outFilterMultimapNmax 10 \
        --alignEndsType EndToEnd \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10 \
        --chimOutType Junctions \
        --outSAMstrandField intronMotif \
        --outTmpDir "$STAR_TMP_DIR" && \
    touch "$SAMPLE_DIR/.status.RNAseq.mapping" && \
    rm -rf "$STAR_TMP_DIR" && \
    echo "[STEP 4] STAR mapping completed successfully."
fi

###########################################
# STEP 5: circRNA calling
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.circRNA" ]; then
    echo "[STEP 5] circRNA calling starting..."
    cd "$SAMPLE_DIR" && \
    awk 'BEGIN{OFS="\t"} ($1==$4 && ($3=="+" || $3=="-") && $1!~/^#/){s=($2<$5)?$2:$5; e=($2<$5)?$5:$2; if((e-s)<=1000000) print $1,s-1,e,"JUNC/"NR"/1",0,$3}' Chimeric.out.junction | sort -k1,1 -k2,2n > bubble.junction.raw && \
    python3 ~/donglab/pipelines/scripts/rnaseq/snap_junctions_to_exons.py \
        "$REFFLAT" \
        bubble.junction.raw \
        bubble.junction && \
    CIRCexplorer2 annotate \
        -r "$REFFLAT" \
        -g "$GENOME_FA" \
        -b bubble.junction \
        -o circularRNA_known.txt \
        --low-confidence \
        --no-fix && \
    if [ -s circularRNA_known.txt ]; then \
        touch "$SAMPLE_DIR/.status.RNAseq.circRNA" && \
        echo "[STEP 5] SUCCESS: Found $(wc -l < circularRNA_known.txt) circRNAs." && \
        if [ ! -f "$SAMPLE_DIR/Aligned.sortedByCoord.out.bam.bai" ]; then \
            samtools index "$SAMPLE_DIR/Aligned.sortedByCoord.out.bam"; \
        fi && \
        python3 ~/donglab/pipelines/scripts/rnaseq/circ_percent_calculation.py \
            "$SAMPLE_DIR/Aligned.sortedByCoord.out.bam" \
            "$SAMPLE_DIR/circularRNA_known.txt" \
            "$SAMPLE_DIR/bubble.junction.raw" && \
        echo "[STEP 5] circRNA calling completed successfully."; \
    else \
        echo "[STEP 5] ERROR: No circRNAs annotated." && exit 1; \
    fi
fi

###########################################
# STEP 6: BAM post-processing and annotation (summary table + BEDs + PDF)
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.bam2annotation" ]; then
    echo "[STEP 6] BAM post-processing and annotation starting..."
    BEDS_DIR="$SAMPLE_DIR/Aligned.sortedByCoord.out_beds"
    mkdir -p "$BEDS_DIR" && \
    for bed in exons genes genic intergenic intergenic_notneargene introns mtRNA non_rRNA_mt rRNA LINE SINE; do \
        src="${ANNOT_BEDS}/${bed}.sorted.bed"; \
        dest="$BEDS_DIR/${bed}.bed"; \
        if [ -f "$src" ]; then cp "$src" "$dest"; else echo "[STEP 6] WARNING: $src missing, creating empty file"; touch "$dest"; fi; \
    done && \
    echo "$(samtools view -cF 0x100 Aligned.sortedByCoord.out.bam) primary alignments" > Aligned.sortedByCoord.out.bam.stat && \
    samtools flagstat Aligned.sortedByCoord.out.bam >> Aligned.sortedByCoord.out.bam.stat && \
    SUMMARY=Aligned.sortedByCoord.out.summary.txt && \
    echo -e "total_non_rRNA_mt\texons\tintrons\trRNA\tmtRNA\tLINE\tSINE\tintergenic\tintergenic_near_genes\tintergenic_not_near_genes\ttotal" > $SUMMARY && \
    count() { samtools view -c -F 0x100 -L "$1" Aligned.sortedByCoord.out.bam; } && \
    total_non_rRNA_mt=$(count "$BEDS_DIR/non_rRNA_mt.bed") && \
    exons=$(count "$BEDS_DIR/exons.bed") && \
    introns=$(count "$BEDS_DIR/introns.bed") && \
    rRNA=$(count "$BEDS_DIR/rRNA.bed") && \
    mtRNA=$(count "$BEDS_DIR/mtRNA.bed") && \
    LINE=$(count "$BEDS_DIR/LINE.bed") && \
    SINE=$(count "$BEDS_DIR/SINE.bed") && \
    intergenic=$(count "$BEDS_DIR/intergenic.bed") && \
    intergenic_not_near_genes=$(count "$BEDS_DIR/intergenic_notneargene.bed") && \
    intergenic_near_genes=$((intergenic - intergenic_not_near_genes)) && \
    total=$(samtools view -c -F 0x100 Aligned.sortedByCoord.out.bam) && \
    echo -e "${total_non_rRNA_mt}\t${exons}\t${introns}\t${rRNA}\t${mtRNA}\t${LINE}\t${SINE}\t${intergenic}\t${intergenic_near_genes}\t${intergenic_not_near_genes}\t${total}" >> $SUMMARY && \
    Rscript "$SCRIPT_DIR/_bam2annotation.r" "$SUMMARY" "Aligned.sortedByCoord.out.bam2annotation.pdf" && \
    touch "$SAMPLE_DIR/.status.RNAseq.bam2annotation" && \
    echo "[STEP 6] BAM annotation completed successfully."
fi

###########################################
# STEP 7: gene counting (htseq-count)
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.htseqcount" ]; then
    echo "[STEP 7] Gene counting (HTSeq) starting..."
    samtools sort -n -o Aligned.sortedByName.bam Aligned.sortedByCoord.out.bam && \
    htseq-count -m intersection-strict -t exon -i gene_id -s yes -q -f bam \
        Aligned.sortedByName.bam "$GTF" \
        > hgseqcount.by.gene.tab 2> hgseqcount.by.gene.tab.stderr && \
    touch "$SAMPLE_DIR/.status.RNAseq.htseqcount" && \
    echo "[STEP 7] Gene counting completed successfully."
fi

###########################################
# STEP 8: LeafCutter junction quantification, clustering, and PSI
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.leafcutter" ]; then
    echo "[STEP 8] LeafCutter junction quantification starting..."
    LEAFCUTTER_BASE="/home/zw529/donglab/pipelines/modules/rnaseq/bin/leafcutter"
    LEAFCUTTER_OUT="$SAMPLE_DIR/leafcutter"
    JUNC_DIR="$LEAFCUTTER_OUT/juncs"
    CLUSTER_DIR="$LEAFCUTTER_OUT/clusters"
    PSI_DIR="$LEAFCUTTER_OUT/psi"
    export PATH="$LEAFCUTTER_BASE/scripts:$LEAFCUTTER_BASE/clustering:$PATH" && \
    mkdir -p "$JUNC_DIR" "$CLUSTER_DIR" "$PSI_DIR" && \
    "$LEAFCUTTER_BASE/scripts/bam2junc.sh" \
        Aligned.sortedByCoord.out.bam \
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
    SS_PSI="$PSI_DIR/${samplename}.splice_site.PSI.tsv" && \
    python "$LEAFCUTTER_BASE/scripts/splice_site_psi.py" \
        "$JUNC_DIR/${samplename}.junc" \
        "$SS_PSI" && \
    { head -n 1 "$SS_PSI" && tail -n +2 "$SS_PSI" | sort -k9,9nr; } > "$PSI_DIR/${samplename}.splice_site.PSI.sorted.tsv" && \
    touch "$SAMPLE_DIR/.status.RNAseq.leafcutter" && \
    echo "[STEP 8] LeafCutter processing completed successfully."
fi

###########################################
# STEP 9: bigWig generation
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.bigwig" ]; then
    echo "[STEP 9] bigWig generation starting..."
    BED="$SAMPLE_DIR/Aligned.sortedByCoord.out.bed"
    if [ ! -s "$BED" ]; then
        bedtools bamtobed -bed12 -i Aligned.sortedByCoord.out.bam > "$BED"
    fi && \
    PRIMARY_READS=$(samtools view -c -F 0x100 Aligned.sortedByCoord.out.bam) && \
    "$SCRIPT_DIR/bam2bigwig.sh" \
        "$BED" \
        -split \
        "$PRIMARY_READS" && \
    touch "$SAMPLE_DIR/.status.RNAseq.bigwig" && \
    echo "[STEP 9] bigWig generation completed successfully."
fi

###########################################
# STEP 10: Cleanup (_val.gz fastqs)
###########################################
if [ -f "$SAMPLE_DIR/.status.RNAseq.bigwig" ]; then
    echo "[STEP 10] Cleaning up trimmed FASTQ files..."
    rm -f "$SAMPLE_DIR"/*_val_1.fq.gz "$SAMPLE_DIR"/*_val_2.fq.gz && \
    echo "[STEP 10] Trimmed FASTQs removed."
fi

echo "[$(date)] RNAseq pipeline finished successfully and cleanup completed."
