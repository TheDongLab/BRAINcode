#!/bin/bash
###########################################
# RNAseq pipeline (paired-end)
# Author: Xianjun Dong & Zachery Wolfe (Zachery updated)
# Date: 1/30/2026
# Version: 3.14 (perfected CIRCexplorer2 steps and PSI steps in leafcutter, including sorting by highest PSI)
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
chmod +x "$SCRIPT_DIR"/{_callSNP.sh,bam2bigwig.sh,_bam2annotation.sh,_bam2annotation.r,_sam2variation.awk}

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
    fastqc -t $CPU --nogroup -o "$SAMPLE_DIR" "$R1" "$R2"
    touch "$SAMPLE_DIR/.status.RNAseq.fastqc"
fi

###########################################
# STEP 2: Trimming (Trim Galore)
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.trim" ]; then
    trim_galore --paired --cores $CPU --output_dir "$SAMPLE_DIR" "$R1" "$R2"
    touch "$SAMPLE_DIR/.status.RNAseq.trim"
fi

R1="$SAMPLE_DIR/$(basename "$R1" .fastq.gz)_val_1.fq.gz"
R2="$SAMPLE_DIR/$(basename "$R2" .fastq.gz)_val_2.fq.gz"

###########################################
# STEP 3: (REMOVED)
# No manual FASTQ rewriting.
# Trust Trim Galore + aligner invariants.
###########################################

###########################################
# STEP 4: kPAL contamination check
###########################################
KPAL_DIR="$SAMPLE_DIR/kpal_results"
KPAL_OUT="$KPAL_DIR/counts.txt"
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.kpal" ]; then
    mkdir -p "$KPAL_DIR"
    gunzip -c "$R1" > "$KPAL_DIR/temp_R1.fastq"
    gunzip -c "$R2" > "$KPAL_DIR/temp_R2.fastq"
    "$HOME/.conda/envs/RNAseq/bin/kpal" count "$KPAL_DIR/temp_R1.fastq" "$KPAL_DIR/temp_R2.fastq" "$KPAL_OUT"
    rm -f "$KPAL_DIR/temp_R1.fastq" "$KPAL_DIR/temp_R2.fastq"
    touch "$SAMPLE_DIR/.status.RNAseq.kpal"
fi

###########################################
# STEP 5: STAR mapping — modified for standard junctions
###########################################
STAR_TMP_DIR="$SAMPLE_DIR/STARtmp_${SLURM_JOB_ID:-NA}_${SLURM_ARRAY_TASK_ID:-NA}"
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.mapping" ]; then
    STAR --runThreadN $CPU \
         --genomeDir /home/zw529/donglab/pipelines/genome \
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
         --outTmpDir "$STAR_TMP_DIR"
    touch "$SAMPLE_DIR/.status.RNAseq.mapping"
    rm -rf "$STAR_TMP_DIR"
fi

###########################################
# STEP 6: circRNA calling
###########################################
cd "$SAMPLE_DIR"
if [ ! -f .status.RNAseq.circRNA ]; then
    echo "[STEP 6] circRNA calling start"

    # 6.1 Create the compatible bubble.junction
    # This version correctly handles coordinate logic and naming requirements
    awk 'BEGIN{OFS="\t"} ($1==$4 && ($3=="+" || $3=="-") && $1!~/^#/){
        s = ($2 < $5) ? $2 : $5;
        e = ($2 < $5) ? $5 : $2;
        if ((e - s) <= 1000000) {
            print $1, s-1, e, "JUNC/"NR"/1", 0, $3
        }
    }' Chimeric.out.junction | sort -k1,1 -k2,2n > bubble.junction

    echo "[STEP 6] Manual BSJ count = $(wc -l < bubble.junction)"

    # 6.2 Run annotation using the proven refFlat reference
    CIRCexplorer2 annotate \
        -r /home/zw529/donglab/pipelines/genome/refFlat.txt \
        -g /home/zw529/donglab/pipelines/genome/hg38.fa.with_chrM \
        -b bubble.junction \
        -o circularRNA_known.txt \
        --low-confidence

    # 6.3 Validation
    if [ -s circularRNA_known.txt ]; then
        touch .status.RNAseq.circRNA
        echo "[STEP 6] SUCCESS: Annotated $(wc -l < circularRNA_known.txt) circRNAs."
        echo "--- TOP RESULTS ---"
        head -n 5 circularRNA_known.txt
    else
        echo "[STEP 6] ERROR: Final validation failed. circularRNA_known.txt is empty."
        exit 1
    fi
fi

###########################################
# STEP 7: BAM post-processing and annotation (summary table + BEDs + PDF)
###########################################
if [ ! -f "$SAMPLE_DIR/.status.RNAseq.bam2annotation" ]; then
    ANNOT_DIR=~/donglab/pipelines/genome/annotation_safe
    BEDS_DIR="$SAMPLE_DIR/Aligned.sortedByCoord.out_beds"
    mkdir -p "$BEDS_DIR"

    # Copy all annotation-safe BEDs
    for bed in exons genes genic intergenic intergenic_notneargene introns mtRNA non_rRNA_mt rRNA LINE SINE; do
        src="$ANNOT_DIR/${bed}.sorted.bed"
        dest="$BEDS_DIR/${bed}.bed"
        if [ -f "$src" ]; then
            cp "$src" "$dest"
        else
            echo "[STEP 7] WARNING: $src missing, creating empty file"
            touch "$dest"
        fi
    done

    # Index BAM
    samtools index Aligned.sortedByCoord.out.bam

    # Generate BAM stats
    echo "$(samtools view -cF 0x100 Aligned.sortedByCoord.out.bam) primary alignments" > Aligned.sortedByCoord.out.bam.stat
    samtools flagstat Aligned.sortedByCoord.out.bam >> Aligned.sortedByCoord.out.bam.stat

    # Generate summary table (overlaps)
    SUMMARY=Aligned.sortedByCoord.out.summary.txt
    echo -e "total_non_rRNA_mt\texons\tintrons\trRNA\tmtRNA\tLINE\tSINE\tintergenic\tintergenic_near_genes\tintergenic_not_near_genes\ttotal" > $SUMMARY

    # Count reads overlapping each BED
    count() { samtools view -c -F 0x100 -L "$1" Aligned.sortedByCoord.out.bam; }

    total_non_rRNA_mt=$(count "$BEDS_DIR/non_rRNA_mt.bed")
    exons=$(count "$BEDS_DIR/exons.bed")
    introns=$(count "$BEDS_DIR/introns.bed")
    rRNA=$(count "$BEDS_DIR/rRNA.bed")
    mtRNA=$(count "$BEDS_DIR/mtRNA.bed")
    LINE=$(count "$BEDS_DIR/LINE.bed")
    SINE=$(count "$BEDS_DIR/SINE.bed")
    intergenic=$(count "$BEDS_DIR/intergenic.bed")
    intergenic_not_near_genes=$(count "$BEDS_DIR/intergenic_notneargene.bed")
    intergenic_near_genes=$((intergenic - intergenic_not_near_genes))
    total=$(samtools view -c -F 0x100 Aligned.sortedByCoord.out.bam)

    echo -e "${total_non_rRNA_mt}\t${exons}\t${introns}\t${rRNA}\t${mtRNA}\t${LINE}\t${SINE}\t${intergenic}\t${intergenic_near_genes}\t${intergenic_not_near_genes}\t${total}" >> $SUMMARY

    # Generate PDF with _bam2annotation.r
    Rscript "$SCRIPT_DIR/_bam2annotation.r" "$SUMMARY" "Aligned.sortedByCoord.out.bam2annotation.pdf"

    # Mark step complete
    touch "$SAMPLE_DIR/.status.RNAseq.bam2annotation"
fi

###########################################
# STEP 8: transcript assembly (cufflinks)
###########################################
strandoption="--library-type fr-firststrand"
if [ ! -f .status.RNAseq.cufflinks ]; then
    cufflinks --no-update-check $strandoption -o "$SAMPLE_DIR" -p $CPU \
        -G /home/zw529/donglab/pipelines/genome/hg38.primary.with_chr.gtf \
        -M "" --compatible-hits-norm --multi-read-correct \
        Aligned.sortedByCoord.out.bam
    touch .status.RNAseq.cufflinks
fi

###########################################
# STEP 9: gene counting (htseq-count)
###########################################
if [ ! -f .status.RNAseq.htseqcount ]; then
    samtools sort -n -o Aligned.sortedByName.bam Aligned.sortedByCoord.out.bam
    htseq-count -m intersection-strict -t exon -i gene_id -s yes -q -f bam \
        Aligned.sortedByName.bam /home/zw529/donglab/pipelines/genome/hg38.primary.with_chr.gtf \
        > hgseqcount.by.gene.tab 2> hgseqcount.by.gene.tab.stderr

    touch .status.RNAseq.htseqcount
fi

###########################################
# STEP 10: LeafCutter junction quantification, clustering, and PSI
###########################################

LEAFCUTTER_BASE="/home/zw529/donglab/pipelines/modules/rnaseq/bin/leafcutter"
LEAFCUTTER_OUT="$SAMPLE_DIR/leafcutter"
JUNC_DIR="$LEAFCUTTER_OUT/juncs"
CLUSTER_DIR="$LEAFCUTTER_OUT/clusters"
PSI_DIR="$LEAFCUTTER_OUT/psi"

export PATH="$LEAFCUTTER_BASE/scripts:$LEAFCUTTER_BASE/clustering:$PATH"

if [ ! -f "$SAMPLE_DIR/.status.RNAseq.leafcutter" ]; then
    echo "[`date`] STEP 10. LeafCutter junction quantification, clustering, and PSI"

    mkdir -p "$JUNC_DIR" "$CLUSTER_DIR" "$PSI_DIR"

    # 10.1 Generate .junc from BAM:
    "$LEAFCUTTER_BASE/scripts/bam2junc.sh" \
        Aligned.sortedByCoord.out.bam \
        "$JUNC_DIR/${samplename}.junc"

    # 10.2 Quantify junctions (per-sample counts):
    python "$LEAFCUTTER_BASE/clustering/leafcutter_quant_only.py" \
        -j "$JUNC_DIR"/*.junc \
        -o "$LEAFCUTTER_OUT/${samplename}"

    # 10.3 Prepare junction list for clustering:
    JUNC_LIST="$JUNC_DIR/juncfile_list.txt"
    ls "$JUNC_DIR"/*.junc > "$JUNC_LIST"

    # 10.4 Cluster junctions (strand-aware):
    RUNDIR="$CLUSTER_DIR"
    PREFIX="leafcutter_clusters"

    python "$LEAFCUTTER_BASE/clustering/leafcutter_cluster.py" \
        -j "$JUNC_LIST" \
        -o "$PREFIX" \
        -r "$RUNDIR"

    # 10.5 Compute PSI per junction (strand from .junc):
    # uses refined clusters by default
    python "$LEAFCUTTER_BASE/scripts/leafcutter_psi.py" \
        "$CLUSTER_DIR/leafcutter_clusters_refined" \
        "$JUNC_DIR/${samplename}.junc" \
        "$PSI_DIR/${samplename}.leafcutter.PSI.tsv"
        
    # 10.6 Sort PSI highest → lowest (keep header)
    PSI_IN="$PSI_DIR/${samplename}.leafcutter.PSI.tsv"
    PSI_OUT="$PSI_DIR/${samplename}.leafcutter.PSI.sorted.tsv"

    { head -n 1 "$PSI_IN" && tail -n +2 "$PSI_IN" | sort -k7,7nr; } > "$PSI_OUT"

    touch "$SAMPLE_DIR/.status.RNAseq.leafcutter"
fi

echo "[`date`] RNAseq pipeline finished successfully."
