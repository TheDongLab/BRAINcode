#!/bin/bash
###########################################
# RNAseq pipeline (paired-end)
# Author: Xianjun Dong & Zachery Wolfe (Zachery updated)
# Date: 1/28/2026
# Version: 3.12 (added leafcutter junction-counting and clustering scripts)
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
# STEP 3: Paired-end filtering
###########################################
R1_FIXED="$SAMPLE_DIR/$(basename "$R1" .fq.gz)_fixed.fq.gz"
R2_FIXED="$SAMPLE_DIR/$(basename "$R2" .fq.gz)_fixed.fq.gz"

if [ ! -f "$SAMPLE_DIR/.status.RNAseq.paired_filter" ]; then
    TMP_R1="$SAMPLE_DIR/tmp_R1.fq"
    TMP_R2="$SAMPLE_DIR/tmp_R2.fq"

    paste <(zcat "$R1") <(zcat "$R2") | \
    awk -F'\t' -v r1="$TMP_R1" -v r2="$TMP_R2" '
        BEGIN {passed=0; failed=0}
        NR%4==1 {h1=$1; h2=$2; gsub(/[^ -~]/,"",h1); gsub(/[^ -~]/,"",h2)}
        NR%4==2 {s1=$1; s2=$2; gsub(/[^ACGTNacgtn]/,"",s1); gsub(/[^ACGTNacgtn]/,"",s2)}
        NR%4==3 {p1=$1; p2=$2}
        NR%4==0 {
            q1=$1; q2=$2; gsub(/[^ -~]/,"",q1); gsub(/[^ -~]/,"",q2)
            if(length(s1)==length(q1) && length(s2)==length(q2) && length(s1)>0) {
                print h1"\n"s1"\n"p1"\n"q1 >> r1
                print h2"\n"s2"\n"p2"\n"q2 >> r2
                passed++
            } else {failed++}
        }
        END {print "Passed:", passed > "/dev/stderr"; print "Failed:", failed > "/dev/stderr"}'

    gzip -c "$TMP_R1" > "$R1_FIXED"
    gzip -c "$TMP_R2" > "$R2_FIXED"
    rm -f "$TMP_R1" "$TMP_R2"
    touch "$SAMPLE_DIR/.status.RNAseq.paired_filter"
fi

R1="$R1_FIXED"
R2="$R2_FIXED"

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
# STEP 6: circRNA calling (CIRCexplorer2)
###########################################
cd "$SAMPLE_DIR"
if [ ! -f .status.RNAseq.circRNA ]; then
    echo "[STEP 6] circRNA calling start"

    echo "[STEP 6] running CIRCexplorer2 parse"
    CIRCexplorer2 parse -t STAR Chimeric.out.junction > circ.parse.txt

    echo "[STEP 6] building BSJ BED (CE2-compatible)"
    awk '$1==$4 && $2>$5 && ($2-$5)<=1e6 {
        OFS="\t";
        print $1,$5-1,$2,"BSJ_"NR"/1",0,$3
    }' Chimeric.out.junction > back_spliced_junction.bed

    echo "[STEP 6] BSJ count = $(wc -l < back_spliced_junction.bed)"

    echo "[STEP 6] running CIRCexplorer2 annotate"
    CIRCexplorer2 annotate \
        -r /home/zw529/donglab/pipelines/genome/hg38.primary.with_chr.ce2.clean.txt \
        -g /home/zw529/donglab/pipelines/genome/hg38.fa.with_chrM \
        -b back_spliced_junction.bed \
        -o circularRNA_known.txt \
        --low-confidence \
        > .CIRCexplorer2_annotate.log 2>&1

    touch .status.RNAseq.circRNA
    echo "[STEP 6] circRNA calling complete"
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
# STEP 10: LeafCutter junction quantification & clustering
###########################################

LEAFCUTTER_BASE="/home/zw529/donglab/pipelines/modules/rnaseq/bin/leafcutter"
LEAFCUTTER_OUT="$SAMPLE_DIR/leafcutter"
JUNC_DIR="$LEAFCUTTER_OUT/juncs"
CLUSTER_DIR="$LEAFCUTTER_OUT/clusters"

export PATH="$LEAFCUTTER_BASE/scripts:$LEAFCUTTER_BASE/clustering:$PATH"

if [ ! -f "$SAMPLE_DIR/.status.RNAseq.leafcutter" ]; then
    echo "[`date`] STEP 10. LeafCutter junction quantification & clustering"

    mkdir -p "$JUNC_DIR" "$CLUSTER_DIR"

    # Generate .junc file from BAM
    "$LEAFCUTTER_BASE/scripts/bam2junc.sh" \
        Aligned.sortedByCoord.out.bam \
        "$JUNC_DIR/${samplename}.junc"

    # Quantify junctions per sample
    python "$LEAFCUTTER_BASE/clustering/leafcutter_quant_only.py" \
        -j "$JUNC_DIR"/*.junc \
        -o "$LEAFCUTTER_OUT/${samplename}"

    # Prepare a text file listing all .junc files for clustering
    JUNC_LIST="$JUNC_DIR/juncfile_list.txt"
    ls "$JUNC_DIR"/*.junc > "$JUNC_LIST"

    # Cluster junctions across shared splice sites
    RUNDIR="$CLUSTER_DIR"
    PREFIX="leafcutter_clusters"

    python "$LEAFCUTTER_BASE/clustering/leafcutter_cluster.py" \
        -j "$JUNC_LIST" \
        -o "$PREFIX" \
        -r "$RUNDIR"

    touch "$SAMPLE_DIR/.status.RNAseq.leafcutter"
fi

echo "[`date`] RNAseq pipeline finished successfully."
