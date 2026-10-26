#!/bin/bash
###########################################
# RNA-seq pipeline (paired-end, stranded)
# Supports TargetALS/FASTQ and AnswerALS/original-BAM inputs.
###########################################

set -euo pipefail

export xml_catalog_files_libxml2="${xml_catalog_files_libxml2:-}"
export CONDA_BACKUP_CXX="${CONDA_BACKUP_CXX:-}"
export CONDA_BACKUP_GXX="${CONDA_BACKUP_GXX:-}"

###########################################
# Inputs
###########################################
R1="${1:-}"
R2="${2:-}"
outputdir="${3:?Usage: _RNAseq.sh R1 R2 outputdir [target_user] [answerALS_original_bam]}"
TARGET_USER="${4:-$USER}"

mkdir -p "$outputdir"

CPU="${CPU:-4}"
samplename=$(basename "$outputdir")
SAMPLE_DIR="$outputdir"

###########################################
# Detect AnswerALS from its output directory.
# TargetALS/FASTQ behavior remains the default.
###########################################
ANSWERALS_MODE=false
case "$(dirname "$(readlink -f "$SAMPLE_DIR")")" in
    */answer_ALS/iPSC_Motor_Neuron/RNAseq/Processed)
        ANSWERALS_MODE=true
        ;;
esac

###########################################
# Conda
###########################################
set +u
source "/home/${TARGET_USER}/donglab/pipelines/modules/miniconda3/etc/profile.d/conda.sh"
conda activate RNAseq
set -u

###########################################
# PATH and script locations
###########################################
export PATH="/home/${TARGET_USER}/.conda/envs/bulkRNAseq/bin:/apps/software/2022b/software/Kent_tools/468-GCC-13.3.0/bin:/opt/slurm/current/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
chmod +x "$SCRIPT_DIR/bam2bigwig.sh" "$SCRIPT_DIR/_bam2annotation.r" "$SCRIPT_DIR/normalize_counts.py"

###########################################
# Reference paths
###########################################
GENOME_BASE="/home/${TARGET_USER}/donglab/references/genome/Homo_sapiens/UCSC/hg38"
GENCODE_DIR="${GENOME_BASE}/Annotation/gencode"
STAR_GENOME="${GENOME_BASE}/Sequence/STAR"
GENOME_FA="${GENOME_BASE}/Sequence/WholeGenomeFasta/genome.fa"
GENOME_FAI="${GENOME_BASE}/Sequence/WholeGenomeFasta/genome.fa.fai"
GTF="${GENCODE_DIR}/gencode.v49.annotation.gtf"
REFFLAT="${GENCODE_DIR}/refFlat_cIRS7.txt"
ANNOT_BEDS="$STAR_GENOME"
CHROMINFO="${STAR_GENOME}/genome.sizes"

# AnswerALS reference files contain only chromosomes 1–22, X, Y, and MT.
ANSWERALS_ANNOT_DIR="${GENCODE_DIR}/answerALS_nochr"
ANSWERALS_CHROMINFO="${STAR_GENOME}/genome.sizes.nochr"

if "$ANSWERALS_MODE"; then
    GTF="${ANSWERALS_ANNOT_DIR}/gencode.v49.annotation.nochr.gtf"
    ANNOT_BEDS="$ANSWERALS_ANNOT_DIR"
    CHROMINFO="$ANSWERALS_CHROMINFO"

    for required in \
        "$GTF" \
        "$CHROMINFO" \
        "$ANSWERALS_ANNOT_DIR/ncbi_to_ucsc_no_chr.txt" \
        "$ANNOT_BEDS/non_rRNA_mt.sorted.bed" \
        "$ANNOT_BEDS/exons.sorted.bed" \
        "$ANNOT_BEDS/introns.sorted.bed" \
        "$ANNOT_BEDS/rRNA.sorted.bed" \
        "$ANNOT_BEDS/mtRNA.sorted.bed" \
        "$ANNOT_BEDS/LINE.sorted.bed" \
        "$ANNOT_BEDS/SINE.sorted.bed" \
        "$ANNOT_BEDS/intergenic.sorted.bed" \
        "$ANNOT_BEDS/intergenic_notneargene.sorted.bed"
    do
        [[ -s "$required" ]] || {
            echo "[ERROR] Required AnswerALS reference file is missing or empty: $required" >&2
            exit 1
        }
    done
fi

###########################################
# Select or create the input BAM
###########################################
STAR_BAM="$SAMPLE_DIR/STAR.Aligned.sortedByCoord.out.bam"
STAR_CREATED_BAM=false
EXISTING_BAM=""

if "$ANSWERALS_MODE"; then
    [[ -n "${5:-}" && -f "${5:-}" ]] || {
        echo "[ERROR] AnswerALS requires the selected original BAM as argument 5." >&2
        exit 1
    }

    EXISTING_BAM="$(readlink -f "$5")"

    [[ "$(dirname "$EXISTING_BAM")" == "$(readlink -f "$SAMPLE_DIR")" ]] || {
        echo "[ERROR] Selected AnswerALS BAM must be in this donor directory." >&2
        exit 1
    }

    case "$(basename "$EXISTING_BAM")" in
        *release5*|*remap*|STAR.*)
            echo "[ERROR] Use the selected original AnswerALS BAM, not a release5, remap, or STAR BAM." >&2
            exit 1
            ;;
        *.sorted.bam)
            ;;
        *)
            echo "[ERROR] Expected an original AnswerALS .sorted.bam input." >&2
            exit 1
            ;;
    esac

    [[ ! -e "$STAR_BAM" || -L "$STAR_BAM" ]] || {
        echo "[ERROR] Refusing to replace a real file at $STAR_BAM." >&2
        exit 1
    }
else
    # Existing TargetALS behavior: discover a pre-existing BAM and skip Steps 1–4.
    while IFS= read -r candidate; do
        if [[ -f "$candidate" && ! "$candidate" =~ remap ]]; then
            EXISTING_BAM="$candidate"
            break
        fi
    done < <(find "$SAMPLE_DIR" -maxdepth 1 -type f -name "*.bam" 2>/dev/null)
fi

if [[ -n "$EXISTING_BAM" ]]; then
    echo "[INFO] Pre-existing BAM found: $(basename "$EXISTING_BAM") — skipping Steps 1–4."
    ln -sf "$EXISTING_BAM" "$STAR_BAM"
    [[ ! -f "${EXISTING_BAM}.bai" ]] || ln -sf "${EXISTING_BAM}.bai" "${STAR_BAM}.bai"
    BAM="$STAR_BAM"
else
    ###########################################
    # STEP 1: FastQC
    ###########################################
    if [[ ! -f "$SAMPLE_DIR/.status.RNAseq.fastqc" ]]; then
        echo "[STEP 1] FastQC quality check starting..."
        fastqc -t "$CPU" --nogroup -o "$SAMPLE_DIR" "$R1" "$R2"
        touch "$SAMPLE_DIR/.status.RNAseq.fastqc"
        echo "[STEP 1] FastQC completed successfully."
    fi

    ###########################################
    # STEP 2: Trim Galore
    ###########################################
    if [[ ! -f "$SAMPLE_DIR/.status.RNAseq.trim" ]]; then
        echo "[STEP 2] Read trimming starting..."
        trim_galore --paired --cores "$CPU" --output_dir "$SAMPLE_DIR" "$R1" "$R2"
        touch "$SAMPLE_DIR/.status.RNAseq.trim"
        echo "[STEP 2] Trimming completed successfully."
    fi

    R1="$SAMPLE_DIR/$(basename "$R1" .fastq.gz)_val_1.fq.gz"
    R2="$SAMPLE_DIR/$(basename "$R2" .fastq.gz)_val_2.fq.gz"

    ###########################################
    # STEP 3: kPAL
    ###########################################
    if [[ ! -f "$SAMPLE_DIR/.status.RNAseq.kpal" ]]; then
        echo "[STEP 3] kPAL contamination check starting..."
        KPAL_DIR="$SAMPLE_DIR/kpal_results"
        mkdir -p "$KPAL_DIR"

        gunzip -c "$R1" > "$KPAL_DIR/temp_R1.fastq"
        gunzip -c "$R2" > "$KPAL_DIR/temp_R2.fastq"

        "/home/${TARGET_USER}/.conda/envs/RNAseq/bin/kpal" count \
            "$KPAL_DIR/temp_R1.fastq" \
            "$KPAL_DIR/temp_R2.fastq" \
            "$KPAL_DIR/counts.k9.txt"

        rm -f "$KPAL_DIR"/temp_R{1,2}.fastq
        touch "$SAMPLE_DIR/.status.RNAseq.kpal"
        echo "[STEP 3] kPAL contamination check completed successfully."
    fi

    ###########################################
    # STEP 4: STAR mapping
    ###########################################
    if [[ ! -f "$SAMPLE_DIR/.status.RNAseq.mapping" ]]; then
        echo "[STEP 4] STAR mapping starting..."
        STAR_TMP_DIR="$SAMPLE_DIR/STARtmp_${SLURM_JOB_ID:-NA}_${SLURM_ARRAY_TASK_ID:-NA}"

        STAR --runThreadN "$CPU" \
            --genomeDir "$STAR_GENOME" \
            --readFilesIn "$R1" "$R2" \
            --readFilesCommand zcat \
            --outFileNamePrefix "$SAMPLE_DIR/STAR." \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattrRGline ID:"$samplename" SM:"$samplename" LB:"$samplename" PL:ILLUMINA PU:"$samplename" \
            --outFilterMultimapNmax 10 \
            --alignEndsType Local \
            --chimSegmentMin 10 \
            --chimJunctionOverhangMin 10 \
            --chimOutType Junctions \
            --outSAMstrandField intronMotif \
            --quantMode GeneCounts \
            --outTmpDir "$STAR_TMP_DIR"

        touch "$SAMPLE_DIR/.status.RNAseq.mapping"
        rm -rf "$STAR_TMP_DIR"
        echo "[STEP 4] STAR mapping completed successfully."
    fi

    [[ -f "$STAR_BAM" ]] || {
        echo "[ERROR] STAR BAM not found after Step 4 — mapping may have failed." >&2
        exit 1
    }

    BAM="$STAR_BAM"
    STAR_CREATED_BAM=true
fi

###########################################
# STEP 5: circRNA calling
###########################################
if [[ ! -f "$SAMPLE_DIR/.status.RNAseq.circRNA" ]]; then
    echo "[STEP 5] circRNA calling starting..."

    if [[ ! -f "$SAMPLE_DIR/STAR.Chimeric.out.junction" ]]; then
        echo "[STEP 5] WARNING: No STAR.Chimeric.out.junction found — skipping circRNA calling."
        touch "$SAMPLE_DIR/.status.RNAseq.circRNA"
    else
        cd "$SAMPLE_DIR"

        CIRCexplorer2 parse -t STAR STAR.Chimeric.out.junction > back_spliced_junction.txt

        CIRCexplorer2 annotate \
            -r "$REFFLAT" \
            -g "$GENOME_FA" \
            -b back_spliced_junction.bed \
            -o circularRNA_known.txt \
            --low-confidence

        if [[ ! -s circularRNA_known.txt ]]; then
            echo "[STEP 5] ERROR: No circRNAs annotated." >&2
            exit 1
        fi

        [[ -f "${STAR_BAM}.bai" ]] || samtools index "$STAR_BAM"

        echo "[STEP 5a] Extracting splice-boundary coordinates..."
        if "$ANSWERALS_MODE"; then
            python3 - "$STAR_BAM" "$GENOME_FAI" <<'ANSWERALS_BOUNDARIES'
import subprocess
import sys

bam, reference_fai = sys.argv[1:]
header = subprocess.check_output(
    ["samtools", "view", "-H", bam], universal_newlines=True
)
sequences = []
for line in header.splitlines():
    if line.startswith("@SQ\t"):
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:])
        sequences.append((fields["SN"], int(fields["LN"])))

lengths = dict(sequences)
order = {name: i for i, (name, _) in enumerate(sequences)}
if not sequences or len(lengths) != len(sequences):
    raise ValueError("Missing or duplicate BAM reference sequences")

with open(reference_fai) as handle:
    reference_lengths = {
        fields[0]: int(fields[1])
        for fields in (line.split() for line in handle)
    }

boundaries, seen = [], set()
with open("circularRNA_known.txt") as handle:
    for line in handle:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 18:
            raise ValueError("Expected headerless CIRCexplorer annotations with 18 fields")

        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        bam_chrom = (
            chrom if chrom in lengths
            else chrom[3:] if chrom.startswith("chr") and chrom[3:] in lengths
            else None
        )
        if bam_chrom is None:
            raise ValueError("No matching BAM chromosome for " + chrom)
        if reference_lengths.get(chrom) != lengths[bam_chrom]:
            raise ValueError("Reference/BAM chromosome length mismatch for " + chrom)
        if not 0 <= start < end <= lengths[bam_chrom]:
            raise ValueError("Invalid circRNA coordinates")

        circ_id = chrom + "_" + str(start) + "_" + str(end)
        if circ_id in seen:
            raise ValueError("Duplicate circRNA coordinate ID: " + circ_id)
        seen.add(circ_id)

        for left, right in ((start, start + 1), (end - 1, end)):
            boundaries.append((bam_chrom, left, right, circ_id, "0", fields[5]))

boundaries.sort(key=lambda row: (order[row[0]], row[1], row[2], row[3]))

with open("temp.ends.genome", "w") as handle:
    for name, length in sequences:
        handle.write("%s\t%d\n" % (name, length))

with open("temp.ends.bed", "w") as handle:
    for row in boundaries:
        handle.write("\t".join(map(str, row)) + "\n")
ANSWERALS_BOUNDARIES
        else
            awk 'BEGIN {FS="\t"; OFS="\t"} {
                circ_id=$1"_" $2"_" $3;
                print $1, $2, $2+1, circ_id, 0, $6;
                print $1, $3-1, $3, circ_id, 0, $6;
            }' circularRNA_known.txt |
                sortBed -i stdin -faidx "$GENOME_FAI" > temp.ends.bed
        fi

        echo "[STEP 5b] Counting linear reads at circRNA boundaries..."
        if "$ANSWERALS_MODE"; then
            bedtools coverage -split \
                -b "$STAR_BAM" -a temp.ends.bed -sorted -counts \
                -g temp.ends.genome > temp.ends.counts
        else
            bedtools coverage -split \
                -b "$STAR_BAM" -a temp.ends.bed -sorted -counts \
                -g "$GENOME_FAI" > temp.ends.counts
        fi

        echo "[STEP 5c] Calculating circularization ratios..."
        python3 "/home/${TARGET_USER}/donglab/pipelines/scripts/rnaseq/circ_percent_calculation.py" \
            circularRNA_known.txt \
            temp.ends.counts \
            circularRNA_known_circ_percentage.txt

        rm -f temp.ends.bed temp.ends.counts
        "$ANSWERALS_MODE" && rm -f temp.ends.genome

        {
            echo -e "chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\texonCount\texonSizes\texonOffsets\treadNumber\tcircType\tgeneName\tisoformName\tindex\tflankIntron\tcirc_reads\tlinear_reads\tcirc_percent"
            cat circularRNA_known_circ_percentage.txt
        } > circularRNA_known_circ_percentage.tmp
        mv circularRNA_known_circ_percentage.tmp circularRNA_known_circ_percentage.txt

        {
            echo -e "chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\texonCount\texonSizes\texonOffsets\treadNumber\tcircType\tgeneName\tisoformName\tindex\tflankIntron"
            cat low_conf_circularRNA_known.txt
        } > low_conf_circularRNA_known.tmp
        mv low_conf_circularRNA_known.tmp low_conf_circularRNA_known.txt

        touch "$SAMPLE_DIR/.status.RNAseq.circRNA"
        echo "[STEP 5] SUCCESS: Found $(tail -n +2 circularRNA_known_circ_percentage.txt | wc -l) circRNAs."
        echo "[STEP 5] circRNA calling completed successfully."
    fi
fi

###########################################
# STEP 6: BAM post-processing and annotation
###########################################
if [[ ! -f "$SAMPLE_DIR/.status.RNAseq.bam2annotation" ]]; then
    echo "[STEP 6] BAM post-processing and annotation starting..."

    SUMMARY="$SAMPLE_DIR/bam.summary.txt"
    echo -e "total_non_rRNA_mt\texons\tintrons\trRNA\tmtRNA\tLINE\tSINE\tintergenic\tintergenic_near_genes\tintergenic_not_near_genes\ttotal" > "$SUMMARY"

    count() { samtools view -c -F 0x100 -L "$1" "$BAM"; }

    total_non_rRNA_mt=$(count "$ANNOT_BEDS/non_rRNA_mt.sorted.bed")
    exons=$(count "$ANNOT_BEDS/exons.sorted.bed")
    introns=$(count "$ANNOT_BEDS/introns.sorted.bed")
    rRNA=$(count "$ANNOT_BEDS/rRNA.sorted.bed")
    mtRNA=$(count "$ANNOT_BEDS/mtRNA.sorted.bed")
    LINE=$(count "$ANNOT_BEDS/LINE.sorted.bed")
    SINE=$(count "$ANNOT_BEDS/SINE.sorted.bed")
    intergenic=$(count "$ANNOT_BEDS/intergenic.sorted.bed")
    intergenic_not_near_genes=$(count "$ANNOT_BEDS/intergenic_notneargene.sorted.bed")
    intergenic_near_genes=$((intergenic - intergenic_not_near_genes))
    total=$(samtools view -c -F 0x100 "$BAM")

    echo -e "${total_non_rRNA_mt}\t${exons}\t${introns}\t${rRNA}\t${mtRNA}\t${LINE}\t${SINE}\t${intergenic}\t${intergenic_near_genes}\t${intergenic_not_near_genes}\t${total}" >> "$SUMMARY"

    samtools view -cF 0x100 "$BAM" > "$SAMPLE_DIR/bam.stat"
    samtools flagstat "$BAM" >> "$SAMPLE_DIR/bam.stat"

    Rscript "$SCRIPT_DIR/_bam2annotation.r" \
        "$SAMPLE_DIR/bam.summary.txt" "$SAMPLE_DIR/bam2annotation.pdf"

    touch "$SAMPLE_DIR/.status.RNAseq.bam2annotation"
    echo "[STEP 6] BAM annotation completed successfully."
fi

###########################################
# STEP 7: Gene counting (HTSeq)
###########################################
if [[ ! -f "$SAMPLE_DIR/.status.RNAseq.htseqcount" ]]; then
    echo "[STEP 7] Gene counting (HTSeq) starting..."

    htseq-count -m intersection-strict -t exon -i gene_id -s reverse -q -f bam -r pos \
        "$BAM" "$GTF" \
        > "$SAMPLE_DIR/htseqcount.tab" 2> "$SAMPLE_DIR/htseqcount.stderr"

    { echo -e "gene_id\tcount"; cat "$SAMPLE_DIR/htseqcount.tab"; } > "$SAMPLE_DIR/htseqcount.tab.tmp"
    mv "$SAMPLE_DIR/htseqcount.tab.tmp" "$SAMPLE_DIR/htseqcount.tab"

    touch "$SAMPLE_DIR/.status.RNAseq.htseqcount"
    echo "[STEP 7] Gene counting completed successfully."
fi

###########################################
# STEP 8: Read normalization
###########################################
if [[ ! -f "$SAMPLE_DIR/.status.RNAseq.normalization" ]]; then
    echo "[STEP 8] Read normalization (TPM/RPKM/FPKM) starting..."

    python3 "$SCRIPT_DIR/normalize_counts.py" \
        "$SAMPLE_DIR/htseqcount.tab" \
        "$GTF" \
        "$SAMPLE_DIR" \
        --library-type PE

    touch "$SAMPLE_DIR/.status.RNAseq.normalization"
    echo "[STEP 8] Normalization completed successfully."
fi

###########################################
# STEP 9: LeafCutter junctions
# AnswerALS workers create per-sample junction files only.
# Cohort clustering and PSI should run after the full sample array finishes.
###########################################
if "$ANSWERALS_MODE"; then
    if [[ ! -f "$SAMPLE_DIR/.status.RNAseq.leafcutter_junctions" ]]; then
        echo "[STEP 9] Extracting AnswerALS LeafCutter junctions..."
        LEAFCUTTER_BASE="/home/${TARGET_USER}/donglab/pipelines/modules/rnaseq/bin/leafcutter"
        LEAFCUTTER_OUT="$SAMPLE_DIR/leafcutter"
        mkdir -p "$LEAFCUTTER_OUT/juncs" "$LEAFCUTTER_OUT/inputs"

        # bam2junc.sh creates temporary files next to its BAM argument.
        LC_BAM="$LEAFCUTTER_OUT/inputs/${samplename}.bam"
        ln -sf "$BAM" "$LC_BAM"
        [[ ! -f "${BAM}.bai" ]] || ln -sf "${BAM}.bai" "${LC_BAM}.bai"

        LC_JUNC="$LEAFCUTTER_OUT/juncs/${samplename}.junc"
        LC_JUNC_RAW="$LEAFCUTTER_OUT/juncs/${samplename}.all_contigs.junc"
        LC_JUNC_FILTERED="${LC_JUNC}.tmp.${SLURM_JOB_ID:-$$}"

        "$LEAFCUTTER_BASE/scripts/bam2junc.sh" \
            "$LC_BAM" \
            "$LC_JUNC_RAW"

        # Keep only junctions on contigs listed in the canonical-only CHROMINFO.
        # AnswerALS BAMs use unprefixed names, matching genome.sizes.nochr.
        awk '
            NR == FNR { allowed[$1] = 1; next }
            NF && ($1 in allowed)
        ' "$CHROMINFO" "$LC_JUNC_RAW" > "$LC_JUNC_FILTERED"

        [[ -s "$LC_JUNC_FILTERED" ]] || {
            echo "[ERROR] LeafCutter produced no junctions on canonical contigs." >&2
            echo "[ERROR] Raw junction file retained at: $LC_JUNC_RAW" >&2
            exit 1
        }

        mv "$LC_JUNC_FILTERED" "$LC_JUNC"
        rm -f "$LC_JUNC_RAW"

        touch "$SAMPLE_DIR/.status.RNAseq.leafcutter_junctions"
        echo "[STEP 9] AnswerALS canonical-contig junction extraction completed."
    fi
else
    # Keep the existing TargetALS per-sample LeafCutter workflow.
    if [[ ! -f "$SAMPLE_DIR/.status.RNAseq.leafcutter" ]]; then
        echo "[STEP 9] LeafCutter junction quantification starting..."
        LEAFCUTTER_BASE="/home/${TARGET_USER}/donglab/pipelines/modules/rnaseq/bin/leafcutter"
        LEAFCUTTER_OUT="$SAMPLE_DIR/leafcutter"
        JUNC_DIR="$LEAFCUTTER_OUT/juncs"
        CLUSTER_DIR="$LEAFCUTTER_OUT/clusters"
        PSI_DIR="$LEAFCUTTER_OUT/psi"

        export PATH="$LEAFCUTTER_BASE/scripts:$LEAFCUTTER_BASE/clustering:$PATH"
        mkdir -p "$JUNC_DIR" "$CLUSTER_DIR" "$PSI_DIR"

        "$LEAFCUTTER_BASE/scripts/bam2junc.sh" \
            "$BAM" "$JUNC_DIR/${samplename}.junc"

        JUNC_LIST="$JUNC_DIR/juncfile_list.txt"
        ls "$JUNC_DIR"/*.junc > "$JUNC_LIST"

        python "$LEAFCUTTER_BASE/clustering/leafcutter_cluster.py" \
            -j "$JUNC_LIST" \
            -o leafcutter_clusters \
            -r "$CLUSTER_DIR"

        python "$LEAFCUTTER_BASE/clustering/leafcutter_quant_only.py" \
            -j "$JUNC_LIST" \
            -o leafcutter_clusters \
            -r "$CLUSTER_DIR"

        LC_PSI="$PSI_DIR/${samplename}.leafcutter.PSI.tsv"
        python "$LEAFCUTTER_BASE/scripts/leafcutter_psi.py" \
            "$CLUSTER_DIR/leafcutter_clusters_refined" \
            "$JUNC_DIR/${samplename}.junc" \
            "$LC_PSI"

        { head -n 1 "$LC_PSI" && tail -n +2 "$LC_PSI" | sort -k7,7nr; } \
            > "$PSI_DIR/${samplename}.leafcutter.PSI.sorted.tsv"

        touch "$SAMPLE_DIR/.status.RNAseq.leafcutter"
        echo "[STEP 9] LeafCutter processing completed successfully."
    fi
fi

###########################################
# STEP 10: bigWig generation
###########################################
if [[ ! -f "$SAMPLE_DIR/.status.RNAseq.bigwig" ]]; then
    echo "[STEP 10] bigWig generation starting..."

    [[ -f "${BAM}.bai" ]] || samtools index "$BAM"

    if "$ANSWERALS_MODE"; then
        "$SCRIPT_DIR/bam2bigwig.sh" "$BAM" -split "$CHROMINFO"
    else
        # Existing TargetALS call remains unchanged.
        "$SCRIPT_DIR/bam2bigwig.sh" "$BAM" -split
    fi

    touch "$SAMPLE_DIR/.status.RNAseq.bigwig"
    echo "[STEP 10] bigWig generation completed successfully."
fi

###########################################
# STEP 11: Cleanup
###########################################
if [[ -f "$SAMPLE_DIR/.status.RNAseq.bigwig" ]]; then
    echo "[STEP 11] Cleanup starting..."

    rm -f "$SAMPLE_DIR"/*_val_1.fq.gz "$SAMPLE_DIR"/*_val_2.fq.gz
    rm -f "$SAMPLE_DIR"/*.bedGraph
    rm -f "$SAMPLE_DIR"/back_spliced*
    echo "[STEP 11] Trimmed FASTQs, bedGraphs, and circRNA intermediates removed."

    rm -f "$SAMPLE_DIR/circularRNA_known.txt"
    rm -f "$SAMPLE_DIR/htseqcount.tab"
    echo "[STEP 11] Redundant intermediate files removed."

    STAR_LOGS=$(ls "$SAMPLE_DIR"/STAR.Log* 2>/dev/null || true)
    if [[ -n "$STAR_LOGS" ]]; then
        STAR_LOGS_DIR="$SAMPLE_DIR/STAR_logs"
        mkdir -p "$STAR_LOGS_DIR"
        mv "$SAMPLE_DIR"/STAR.Log* "$STAR_LOGS_DIR/"
        zip -j "$STAR_LOGS_DIR/STAR_logs.zip" "$STAR_LOGS_DIR"/STAR.Log*
        rm -f "$STAR_LOGS_DIR"/STAR.Log*
        echo "[STEP 11] STAR logs archived."
    else
        echo "[STEP 11] No STAR logs found — skipping log archiving."
    fi

    if [[ "$STAR_CREATED_BAM" == true && -f "$STAR_BAM" ]]; then
        echo "[STEP 11] Converting STAR BAM to CRAM..."
        STAR_CRAM="${STAR_BAM%.bam}.cram"
        samtools view -@ "$CPU" -C -T "$GENOME_FA" -o "$STAR_CRAM" "$STAR_BAM"
        samtools index "$STAR_CRAM"
        rm -f "$STAR_BAM" "${STAR_BAM}.bai"
        echo "[STEP 11] CRAM conversion complete: $(basename "$STAR_CRAM")"
    fi

    echo "[STEP 11] Cleanup completed."
fi

echo "[$(date)] RNA-seq pipeline finished successfully."
