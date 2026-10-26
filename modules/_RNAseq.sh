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
        "$GENOME_FA" \
        "$GENOME_FAI" \
        "$REFFLAT" \
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

    # The AnswerALS whitelist must contain exactly the 25 primary contigs.
    awk '
        BEGIN { FS="[ \t]+" }
        NF {
            n++
            if (NF < 2 || $2 !~ /^[0-9]+$/ ||
                $1 !~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|MT)$/) {
                bad=1
            }
        }
        END {
            if (bad || n != 25) exit 1
        }
    ' "$CHROMINFO" || {
        echo "[ERROR] AnswerALS CHROMINFO must list exactly chromosomes 1–22, X, Y, and MT: $CHROMINFO" >&2
        exit 1
    }
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
    [[ ! -f "${EXISTING_BAM}.bai" ]] || ln -sf "$EXISTING_BAM.bai" "$STAR_BAM.bai"
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
# AnswerALS canonical-contig setup
###########################################
WORK_TMPDIR=""
CANONICAL_BED=""

if "$ANSWERALS_MODE"; then
    WORK_TMPDIR=$(mktemp -d "${TMPDIR:-/tmp}/answerALS_rnaseq_${SLURM_JOB_ID:-$$}_XXXXXX")
    trap '[[ -z "${WORK_TMPDIR:-}" ]] || rm -rf "$WORK_TMPDIR"' EXIT

    CANONICAL_BED="$WORK_TMPDIR/canonical_contigs.bed"
    awk 'BEGIN { FS="[ \t]+"; OFS="\t" } NF >= 2 { print $1, 0, $2 }' \
        "$CHROMINFO" > "$CANONICAL_BED"

    BAM_REFS="$WORK_TMPDIR/bam.references.tsv"
    samtools view -H "$BAM" |
        awk -F'\t' '$1 == "@SQ" {
            sn=""; ln="";
            for (i=2; i<=NF; i++) {
                if ($i ~ /^SN:/) sn=substr($i,4);
                if ($i ~ /^LN:/) ln=substr($i,4);
            }
            if (sn != "" && ln != "") print sn "\t" ln;
        }' > "$BAM_REFS"

    while IFS=$'\t ' read -r chrom length rest; do
        [[ -n "$chrom" && -n "$length" ]] || continue

        bam_length=$(awk -F'\t' -v chrom="$chrom" '$1 == chrom { print $2; exit }' "$BAM_REFS")
        [[ -n "$bam_length" ]] || {
            echo "[ERROR] Canonical contig '$chrom' is missing from BAM header: $BAM" >&2
            exit 1
        }
        [[ "$bam_length" == "$length" ]] || {
            echo "[ERROR] BAM/reference length mismatch for $chrom: BAM=$bam_length reference=$length" >&2
            exit 1
        }
    done < "$CHROMINFO"

    echo "[INFO] AnswerALS analysis is restricted to the 25 contigs listed in $CHROMINFO."
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

        CIRC_CHIMERIC="$SAMPLE_DIR/STAR.Chimeric.out.junction"

        if "$ANSWERALS_MODE"; then
            # CIRCexplorer's reference uses UCSC-style names (chr1, chrM).
            # Keep only canonical AnswerALS endpoints, then translate their
            # names for CIRCexplorer. The original junction file is untouched.
            CIRC_CHIMERIC="$WORK_TMPDIR/STAR.Chimeric.canonical.ucsc.out.junction"

            awk -F'\t' -v OFS='\t' '
                function to_answerals(chrom, base) {
                    base=chrom
                    sub(/^chr/, "", base)
                    if (base == "M") base="MT"
                    return base
                }
                function to_ucsc(chrom) {
                    if (chrom == "MT") return "chrM"
                    return "chr" chrom
                }
                NR == FNR {
                    allowed[$1]=1
                    next
                }
                NF >= 6 {
                    left=to_answerals($1)
                    right=to_answerals($4)
                    if (!(left in allowed) || !(right in allowed)) next
                    $1=to_ucsc(left)
                    $4=to_ucsc(right)
                    print
                }
            ' "$CHROMINFO" "$SAMPLE_DIR/STAR.Chimeric.out.junction" > "$CIRC_CHIMERIC"

            [[ -s "$CIRC_CHIMERIC" ]] || {
                echo "[ERROR] No canonical-contig chimeric junctions found for circRNA analysis." >&2
                exit 1
            }
        fi

        CIRCexplorer2 parse -t STAR "$CIRC_CHIMERIC" > back_spliced_junction.txt

        CIRCexplorer2 annotate \
            -r "$REFFLAT" \
            -g "$GENOME_FA" \
            -b back_spliced_junction.bed \
            -o circularRNA_known.txt \
            --low-confidence

        if "$ANSWERALS_MODE"; then
            # Convert CIRCexplorer's chr-prefixed annotation records back to
            # AnswerALS BAM names and discard anything outside the whitelist.
            normalize_circ_annotations() {
                local annotation_file="$1"
                local normalized_file="${annotation_file}.canonical.tmp"

                [[ -f "$annotation_file" ]] || : > "$annotation_file"

                awk -F'\t' -v OFS='\t' '
                    function to_answerals(chrom, base) {
                        base=chrom
                        sub(/^chr/, "", base)
                        if (base == "M") base="MT"
                        return base
                    }
                    NR == FNR {
                        allowed[$1]=1
                        next
                    }
                    NF {
                        chrom=to_answerals($1)
                        if (chrom in allowed) {
                            $1=chrom
                            print
                        }
                    }
                ' "$CHROMINFO" "$annotation_file" > "$normalized_file"

                mv "$normalized_file" "$annotation_file"
            }

            normalize_circ_annotations "circularRNA_known.txt"
            normalize_circ_annotations "low_conf_circularRNA_known.txt"
        fi

        if [[ ! -s circularRNA_known.txt ]]; then
            echo "[STEP 5] ERROR: No canonical circRNAs annotated." >&2
            exit 1
        fi

        [[ -f "${STAR_BAM}.bai" ]] || samtools index "$STAR_BAM"

        echo "[STEP 5a] Extracting canonical splice-boundary coordinates..."
        if "$ANSWERALS_MODE"; then
            python3 - "$STAR_BAM" "$GENOME_FAI" "$CHROMINFO" <<'ANSWERALS_BOUNDARIES'
import subprocess
import sys

bam, reference_fai, chrominfo = sys.argv[1:]

allowed = {}
with open(chrominfo) as handle:
    for line in handle:
        fields = line.split()
        if len(fields) >= 2:
            allowed[fields[0]] = int(fields[1])

header = subprocess.check_output(
    ["samtools", "view", "-H", bam], universal_newlines=True
)

bam_sequences = []
for line in header.splitlines():
    if line.startswith("@SQ\t"):
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:])
        name = fields["SN"]
        length = int(fields["LN"])
        if name in allowed:
            if allowed[name] != length:
                raise ValueError("BAM/reference chromosome length mismatch for " + name)
            bam_sequences.append((name, length))

bam_names = [name for name, _ in bam_sequences]
if set(bam_names) != set(allowed):
    missing = sorted(set(allowed) - set(bam_names))
    raise ValueError("Canonical contigs missing from BAM header: " + ", ".join(missing))

lengths = dict(bam_sequences)
order = {name: i for i, (name, _) in enumerate(bam_sequences)}

with open(reference_fai) as handle:
    reference_lengths = {
        fields[0]: int(fields[1])
        for fields in (line.split() for line in handle)
    }

def ucsc_name(name):
    return "chrM" if name == "MT" else "chr" + name

boundaries, seen = [], set()
with open("circularRNA_known.txt") as handle:
    for line in handle:
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 18:
            raise ValueError("Expected headerless CIRCexplorer annotations with 18 fields")

        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        if chrom not in allowed:
            raise ValueError("Noncanonical circRNA contig reached boundary stage: " + chrom)

        if reference_lengths.get(ucsc_name(chrom)) != lengths[chrom]:
            raise ValueError("FASTA/BAM chromosome length mismatch for " + chrom)

        if not 0 <= start < end <= lengths[chrom]:
            raise ValueError("Invalid circRNA coordinates for " + chrom)

        circ_id = chrom + "_" + str(start) + "_" + str(end)
        if circ_id in seen:
            raise ValueError("Duplicate circRNA coordinate ID: " + circ_id)
        seen.add(circ_id)

        for left, right in ((start, start + 1), (end - 1, end)):
            boundaries.append((chrom, left, right, circ_id, "0", fields[5]))

boundaries.sort(key=lambda row: (order[row[0]], row[1], row[2], row[3]))

with open("temp.ends.genome", "w") as handle:
    for name, length in bam_sequences:
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

        echo "[STEP 5b] Counting linear reads at canonical circRNA boundaries..."
        if "$ANSWERALS_MODE"; then
            samtools view -bh -L "$CANONICAL_BED" "$STAR_BAM" |
                bedtools coverage -split \
                    -b stdin -a temp.ends.bed -sorted -counts \
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
        echo "[STEP 5] SUCCESS: Found $(tail -n +2 circularRNA_known_circ_percentage.txt | wc -l) canonical circRNAs."
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

    if "$ANSWERALS_MODE"; then
        total=$(samtools view -c -F 0x100 -L "$CANONICAL_BED" "$BAM")
    else
        total=$(samtools view -c -F 0x100 "$BAM")
    fi

    echo -e "${total_non_rRNA_mt}\t${exons}\t${introns}\t${rRNA}\t${mtRNA}\t${LINE}\t${SINE}\t${intergenic}\t${intergenic_near_genes}\t${intergenic_not_near_genes}\t${total}" >> "$SUMMARY"

    if "$ANSWERALS_MODE"; then
        samtools view -c -F 0x100 -L "$CANONICAL_BED" "$BAM" > "$SAMPLE_DIR/bam.stat"
        samtools view -bh -L "$CANONICAL_BED" "$BAM" |
            samtools flagstat - >> "$SAMPLE_DIR/bam.stat"
    else
        samtools view -cF 0x100 "$BAM" > "$SAMPLE_DIR/bam.stat"
        samtools flagstat "$BAM" >> "$SAMPLE_DIR/bam.stat"
    fi

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

    if "$ANSWERALS_MODE"; then
        samtools view -bh -L "$CANONICAL_BED" "$BAM" |
            htseq-count -m intersection-strict -t exon -i gene_id \
                -s reverse -q -f bam -r pos - "$GTF" \
                > "$SAMPLE_DIR/htseqcount.tab" \
                2> "$SAMPLE_DIR/htseqcount.stderr"
    else
        htseq-count -m intersection-strict -t exon -i gene_id -s reverse -q -f bam -r pos \
            "$BAM" "$GTF" \
            > "$SAMPLE_DIR/htseqcount.tab" 2> "$SAMPLE_DIR/htseqcount.stderr"
    fi

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
# AnswerALS workers create per-sample canonical junction files only.
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
        LC_JUNC_FILTERED="$WORK_TMPDIR/${samplename}.canonical.junc"

        "$LEAFCUTTER_BASE/scripts/bam2junc.sh" \
            "$LC_BAM" \
            "$LC_JUNC_RAW"

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
    # Keep existing TargetALS per-sample LeafCutter workflow.
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
