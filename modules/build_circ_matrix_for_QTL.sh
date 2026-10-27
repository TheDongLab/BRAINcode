#!/bin/bash
#SBATCH --job-name=build_and_plot_circ
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/build_and_plot_circ.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/build_and_plot_circ.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=160G

# Set to target_ALS or answer_ALS before running
EXPERIMENT="answer_ALS"
export EXPERIMENT

# Outputs in the selected experiment's QTL directory:
#   circ_counts_matrix.txt: CIRCexplorer2 readNumber, one row per BSJ.
#   circ_rpm_matrix.txt: readNumber / primary mapped reads from bam.stat * 1e6.
#   circ_rpm_sample_metadata.tsv: denominators, mapping rates, source files, audit.
# Every selected sample uses bam.stat: primary mapped QC-pass + QC-fail reads.
# No STAR-log, total-read, or read-pair fallback; do not divide the denominator by 2.
# Optional STAR.Chimeric.out.junction/back_spliced_junction.bed enable a parser consistency audit; absent audit inputs are recorded, not claimed as verified.
# Original circ_matrix.txt remains a fraction matrix; no INT is applied to RPM.
set -euo pipefail

unset PYTHONPATH

module purge
module load Python/3.12.3-GCCcore-13.3.0
module load Python-bundle-PyPI/2024.06-GCCcore-13.3.0

python3 - <<'EOF'
import os
import pandas as pd
import numpy as np
import re
from pathlib import Path
from collections import defaultdict

pd.set_option("future.no_silent_downcasting", True)

EXPERIMENT = os.environ["EXPERIMENT"]
BASE = Path("/home/zw529/donglab/data") / EXPERIMENT
OUT = BASE / "QTL"
OUT.mkdir(exist_ok=True)

GENE_BED = Path(
    "/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/"
    "hg38/Annotation/gencode/gencode.v49.annotation.gene.bed6"
)

CANONICAL_CHROMS = {
    *(str(i) for i in range(1, 23)),
    "X", "Y", "MT",
}

def norm_id(x):
    return str(x).strip().replace("-", "_")

def answer_chrom(chrom):
    """Normalize AnswerALS canonical contig names to no-chr form."""
    chrom = str(chrom).strip()
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        chrom = "MT"
    return chrom if chrom in CANONICAL_CHROMS else None

# Use the same original-library statistic for every sample. The standalone
# first line and generic "mapped" field in bam.stat are NOT this denominator.
def mapped_read_stats(sample_dir):
    path = sample_dir / "bam.stat"
    if not path.is_file():
        raise ValueError(f"Missing {path}; no alternate RPM denominator is permitted")
    text = path.read_text()
    matches = re.findall(
        r"^\s*(\d+)\s*\+\s*(\d+)\s+primary mapped\s*(?:\(|$)",
        text,
        re.M,
    )
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one primary mapped field in {path}")
    passed, failed = map(int, matches[0])
    mapped = passed + failed
    if mapped <= 0:
        raise ValueError(f"Zero primary mapped reads in {path}")
    totals = re.findall(r"^\s*(\d+)\s*\+\s*(\d+)\s+primary\s*$", text, re.M)
    if len(totals) != 1:
        raise ValueError(f"Expected exactly one primary total field in {path}")
    total = sum(map(int, totals[0]))
    if total < mapped:
        raise ValueError(f"Primary mapped reads exceed primary total in {path}")
    return dict(
        mapping_stats_file=str(path),
        primary_mapped_qc_pass=passed,
        primary_mapped_qc_fail=failed,
        primary_mapped_reads=mapped,
        primary_reads=total,
        primary_mapped_fraction=mapped / total,
        rpm_multiplier=1e6 / mapped,
        denominator_definition="bam.stat_primary_mapped_QCpass_plus_QCfail_reads",
    )

def audit_star_junctions(sample_dir):
    # Compare pre-annotation counts: annotate can realign/change boundaries.
    path = sample_dir / "STAR.Chimeric.out.junction"
    bed = sample_dir / "back_spliced_junction.bed"
    if not path.is_file():
        return dict(junction_audit="unavailable_no_chimeric_file")
    counts = defaultdict(int)
    seen = set()
    duplicates = 0
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.split()
            if len(c) < 14:
                raise ValueError(f"Malformed STAR junction row in {path}")
            if int(c[6]) < 0 or c[0] != c[3] or c[2] != c[5]:
                continue
            start, end = (
                (int(c[4]), int(c[1]) - 1)
                if c[2] == "+"
                else (int(c[1]), int(c[4]) - 1)
            )
            if start > end:
                continue
            key = (c[0], start, end)
            support = (key, c[9])
            duplicates += support in seen
            seen.add(support)
            counts[key] += 1

    # In the unmapped-read single-FASTQ workflow, mates can share a read name.
    # Report repeated names without deduplicating or reinterpreting readNumber.
    status = "junction_records_bed_unavailable"
    if bed.is_file():
        observed = {}
        for line in bed.read_text().splitlines():
            if not line.strip():
                continue
            c = line.split()
            key = (c[0], int(c[1]), int(c[2]))
            if key in observed:
                raise ValueError(f"Duplicate junction in {bed}: {key}")
            observed[key] = int(c[3].rsplit("/", 1)[1])
        if dict(counts) != observed:
            raise ValueError(f"STAR support counts disagree with {bed}; check parser/version")
        status = "junction_records_match_parsed_bed"
    return dict(
        junction_audit=status,
        qualifying_bsj_records=sum(counts.values()),
        repeated_read_name_junction_records=duplicates,
    )

# =========================================================================
# STEP 1: LOAD GENOMIC ANNOTATIONS (BED6)
# =========================================================================
print("Loading GENCODE v49 gene annotations...")
genes_by_chrom = defaultdict(list)
if GENE_BED.exists():
    with open(GENE_BED, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                chrom = parts[0]
                if EXPERIMENT == "answer_ALS":
                    chrom = answer_chrom(chrom)
                    if chrom is None:
                        continue
                start = int(parts[1])
                end = int(parts[2])
                gene_name = parts[3]
                genes_by_chrom[chrom].append((start, end, gene_name))
    print(f"Loaded annotations for {len(genes_by_chrom)} chromosomes.")
else:
    print(
        f"Warning: Gene annotation BED6 file not found at {GENE_BED}. "
        "Table will use 'NA' for gene names."
    )

def find_overlapping_genes(circ_id):
    try:
        parts = circ_id.split(":")
        chrom = parts[0]
        if EXPERIMENT == "answer_ALS":
            chrom = answer_chrom(chrom)
            if chrom is None:
                return "NA"
        coord_parts = parts[1].split("-")
        c_start = int(coord_parts[0])
        c_end = int(coord_parts[1])
    except Exception:
        return "NA"

    if chrom not in genes_by_chrom:
        return "NA"

    overlapping = []
    for g_start, g_end, g_name in genes_by_chrom[chrom]:
        if max(c_start, g_start) < min(c_end, g_end):
            overlapping.append(g_name)

    if overlapping:
        return ",".join(sorted(set(overlapping)))
    return "NA"

# =========================================================================
# STEP 2: MATRIX GENERATION & VARIANCE FILTERING
# =========================================================================
if EXPERIMENT == "target_ALS":
    # Keep existing TargetALS metadata matching and first-match behavior.
    rna = pd.read_csv(BASE / "targetALS_rnaseq_metadata.csv")
    rna.columns = rna.columns.str.strip().str.lower()

    all_circ = list(
        BASE.glob("**/RNAseq/Processed/*/circularRNA_known_circ_percentage.txt")
    )

    def find_circ(sample_id):
        sid_norm = norm_id(sample_id)
        for circ_path in all_circ:
            if norm_id(circ_path.parent.name) == sid_norm:
                return circ_path
        return None

    rna["circ_path"] = rna["externalsampleid"].apply(find_circ)
    rna_found = rna[rna["circ_path"].notna()].copy()
    rna_found = rna_found.drop_duplicates(
        subset=["externalsampleid"],
        keep="first",
    )
    print(f"Cleaned unique RNA samples for matrix processing: {len(rna_found)}")

elif EXPERIMENT == "answer_ALS":
    rna_processed = BASE / "iPSC_Motor_Neuron/RNAseq/Processed"
    all_circ = sorted(
        rna_processed.glob("*/circularRNA_known_circ_percentage.txt")
    )
    rna_found = pd.DataFrame([
        {
            "matrix_sample_id": circ_path.parent.name,
            "donor_id": circ_path.parent.name,
            "tissue": "iPSC_Motor_Neuron",
            "circ_path": circ_path,
        }
        for circ_path in all_circ
    ])
    print(f"AnswerALS donors with circRNA percentage files: {len(rna_found)}")

else:
    raise ValueError("EXPERIMENT must be target_ALS or answer_ALS")

expr = None
reads_list = []
rpm_metadata = []
seen_samples = set()

for _, row in rna_found.iterrows():
    if EXPERIMENT == "target_ALS":
        sample = norm_id(row["externalsampleid"])
    else:
        # Preserve CASE-... IDs to match the AnswerALS joint VCF.
        sample = str(row["matrix_sample_id"]).strip()

    circ_file = row["circ_path"]

    if sample in seen_samples:
        raise ValueError(f"Duplicate sample ID: {sample}")
    seen_samples.add(sample)

    circ = pd.read_csv(circ_file, sep="\t")
    required = {
        "chrom", "start", "end", "strand", "circ_percent", "readNumber"
    }
    if not required.issubset(circ.columns):
        raise ValueError(f"Missing columns in {circ_file}: {required - set(circ.columns)}")

    if EXPERIMENT == "answer_ALS":
        # Keep only the 25 canonical AnswerALS contigs and normalize to no-chr
        circ["chrom"] = circ["chrom"].map(answer_chrom)
        circ = circ[circ["chrom"].notna()].copy()

    counts = pd.to_numeric(circ["readNumber"], errors="raise")
    if not (np.isfinite(counts) & (counts >= 0) & (counts % 1 == 0)).all():
        raise ValueError(f"Invalid readNumber counts in {circ_file}")
    circ["readNumber"] = counts.astype("int64")

    if "circ_reads" in circ.columns:
        if not np.array_equal(
            pd.to_numeric(circ["circ_reads"], errors="raise"),
            counts,
        ):
            raise ValueError(f"circ_reads disagrees with readNumber in {circ_file}")

    stats = mapped_read_stats(circ_file.parent)
    audit = audit_star_junctions(circ_file.parent)
    print(
        f"{sample}: {stats['primary_mapped_reads']:,} primary mapped reads; "
        f"{audit['junction_audit']}"
    )
    rpm_metadata.append(
        dict(sample_id=sample, circ_file=str(circ_file), **stats, **audit)
    )

    circ["circ_id"] = (
        circ["chrom"].astype(str) + ":" +
        circ["start"].astype(str) + "-" +
        circ["end"].astype(str) + ":" +
        circ["strand"].astype(str)
    )

    if circ["circ_id"].duplicated().any():
        raise ValueError(
            f"Duplicate circ junction IDs in {circ_file}; refusing to double count"
        )

    pct_df = circ[["circ_id", "circ_percent"]].copy()
    pct_df.columns = ["circ_id", sample]

    if expr is None:
        expr = pct_df
    else:
        expr = expr.merge(pct_df, on="circ_id", how="outer")

    read_df = circ[["circ_id", "readNumber"]].copy()
    read_df.columns = ["circ_id", sample]
    reads_list.append(read_df)

if expr is None:
    raise ValueError("No usable circRNA samples found")

expr = expr.fillna(0).infer_objects(copy=False)

print("Aggregating cohort-wide back-spliced read counts...")
reads_master = reads_list[0]
for next_df in reads_list[1:]:
    reads_master = reads_master.merge(next_df, on="circ_id", how="outer")

reads_master = reads_master.fillna(0)
read_columns = reads_master.columns.drop("circ_id")
reads_master[read_columns] = reads_master[read_columns].astype("int64")

# Align factors by sample name, never by incidental file order.
rpm_meta = pd.DataFrame(rpm_metadata).set_index("sample_id").loc[read_columns]
rpm_master = (
    reads_master.set_index("circ_id")
    .div(rpm_meta["primary_mapped_reads"], axis=1) * 1e6
)
rpm_master.to_csv(OUT / "circ_rpm_matrix.txt", sep="\t", index=True)
rpm_meta.to_csv(OUT / "circ_rpm_sample_metadata.tsv", sep="\t", index=True)
print(f"Saved BSJ RPM per primary mapped read: {OUT / 'circ_rpm_matrix.txt'}")

# Export the full BSJ count matrix before the circ-fraction SD filter.
reads_master.to_csv(OUT / "circ_counts_matrix.txt", sep="\t", index=False)
print(f"Saved unfiltered BSJ counts: {OUT / 'circ_counts_matrix.txt'}")

circ_ids = expr["circ_id"]
numeric_data = expr.drop("circ_id", axis=1)
row_sds = numeric_data.std(axis=1)

min_sd_threshold = 0.005
variant_mask = row_sds > min_sd_threshold

filtered_numeric = numeric_data[variant_mask].copy()
filtered_circ_ids = circ_ids[variant_mask]
expr_final = pd.concat(
    [
        filtered_circ_ids.reset_index(drop=True),
        filtered_numeric.reset_index(drop=True),
    ],
    axis=1,
)

# =========================================================================
# STEP 3: CALCULATE DETAILED ABUNDANCE METRICS
# =========================================================================
print("Calculating read metrics and mapping host genes...")
reads_filtered = reads_master[
    reads_master["circ_id"].isin(filtered_circ_ids)
].copy()
numeric_reads_filtered = reads_filtered.drop(columns=["circ_id"])

total_reads_filtered = numeric_reads_filtered.sum(axis=1)
mean_reads_filtered = numeric_reads_filtered.mean(axis=1)
mean_reads_nonzero_only = (
    numeric_reads_filtered.replace(0, np.nan).mean(axis=1).fillna(0)
)

master_metrics = pd.DataFrame({
    "circ_id": reads_filtered["circ_id"],
    "total_reads": total_reads_filtered,
    "avg_reads": mean_reads_filtered,
    "avg_reads_nonzero_only": mean_reads_nonzero_only,
})

master_metrics["gene_name"] = master_metrics["circ_id"].apply(
    find_overlapping_genes
)

col_order = [
    "circ_id", "gene_name", "total_reads", "avg_reads",
    "avg_reads_nonzero_only",
]
master_metrics = master_metrics[col_order]
master_metrics_filtered = master_metrics[
    master_metrics["total_reads"] > 0
].copy()

expr_final.to_csv(OUT / "circ_matrix.txt", sep="\t", index=False)

if EXPERIMENT == "target_ALS":
    rna_found[
        ["externalsampleid", "externalsubjectid", "tissue"]
    ].to_csv(OUT / "circ_sample_metadata.csv", index=False)
else:
    rna_found[
        ["matrix_sample_id", "donor_id", "tissue"]
    ].to_csv(OUT / "circ_sample_metadata.csv", index=False)

master_metrics_filtered.to_csv(
    OUT / "circ_abundance_distribution_table.txt",
    sep="\t",
    index=False,
)

reads_summary = pd.DataFrame({
    "circ_id": reads_master["circ_id"],
    "total_reads": reads_master.drop(columns=["circ_id"]).sum(axis=1),
})
reads_summary[
    reads_summary["circ_id"].isin(filtered_circ_ids)
].to_csv(OUT / "filtered_reads_summary.tmp", sep="\t", index=False)

pd.DataFrame({
    "circ_id": reads_filtered["circ_id"],
    "avg_reads": mean_reads_filtered,
}).to_csv(OUT / "filtered_averages_summary.tmp", sep="\t", index=False)

num_circs, num_samples = filtered_numeric.shape
total_cells = filtered_numeric.size
nonzero_cells = (filtered_numeric > 0).sum().sum()

print("\n" + "=" * 40)
print("        === MATRIX STATS ===")
print("=" * 40)
print(f"Total circRNAs: {num_circs:,}")
print(f"Total Samples: {num_samples:,}")
if total_cells > 0:
    print(f"Global Avg: {filtered_numeric.values.mean():.6f}")
if nonzero_cells > 0:
    print(
        f"Detected Avg (>0): "
        f"{filtered_numeric.values[filtered_numeric.values > 0].mean():.6f}"
    )
if total_cells > 0:
    print(f"Sparsity: {(1 - (nonzero_cells / total_cells)) * 100:.2f}%")

active_counts = (filtered_numeric > 0).sum(axis=0)
if len(active_counts):
    top_sample = active_counts.idxmax()
    print(f"Top Subject (>0%): {top_sample} ({active_counts[top_sample]:,} circs)")

print("\n" + "=" * 40)
print("        === CHROMOSOMES ===")
print("=" * 40)
chrs_series = filtered_circ_ids.reset_index(drop=True).str.split(":").str[0]
chr_counts = chrs_series.value_counts()
if EXPERIMENT == "target_ALS":
    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
else:
    chr_order = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
for chrom in chr_order:
    if chrom in chr_counts:
        print(f"{chrom}: {chr_counts[chrom]:,}")
EOF

module --force purge
module load R

Rscript - <<'EOF'
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(svglite)
  library(scales)
})

experiment <- Sys.getenv("EXPERIMENT")
out_dir <- paste0("/home/zw529/donglab/data/", experiment, "/QTL/")
table_path <- paste0(out_dir, "circ_abundance_distribution_table.txt")

print("Loading unified metrics table...")
df_metrics <- read.delim(table_path, header=TRUE, sep="\t")

y_ticks <- 10^(0:5)

als_genes <- c("ATXN1", "ATXN2", "HOMER1", "C9orf72", "SOD1", "FUS", "STMN2", "TARDBP", "TBK1", "UNC13A")
syn_genes <- c("RIMS1", "RIMS2")
target_genes <- c(als_genes, syn_genes)

df_metrics$cat <- ifelse(grepl(paste(als_genes, collapse="|"), df_metrics$gene_name), "ALS",
                  ifelse(grepl(paste(syn_genes, collapse="|"), df_metrics$gene_name), "SYN", NA))

df_targets <- df_metrics[!is.na(df_metrics$cat) & df_metrics$total_reads >= 3, ]
cat_colors <- c("ALS" = "#E41A1C", "SYN" = "#377EB8")

print(paste("Found", nrow(df_targets), "target circRNAs with total reads >= 3 for annotation."))

save_ggplot <- function(plot_obj, basename, width_in = 11.33, height_in = 5.33, dpi = 300) {
    ggsave(paste0(out_dir, basename, ".png"), plot=plot_obj, width=width_in, height=height_in, dpi=dpi, units="in")
    ggsave(paste0(out_dir, basename, ".svg"), plot=plot_obj, width=width_in, height=height_in, units="in", device=svglite::svglite)
}

counts_table_tot <- table(df_metrics$total_reads)
plot_data_tot <- data.frame(
    reads = as.numeric(names(counts_table_tot)),
    circ_count = as.numeric(counts_table_tot)
)

df_tot_sorted <- df_metrics[order(-df_metrics$total_reads), ]
top_outliers_tot <- head(df_tot_sorted, 5)
outliers_p2_tot <- top_outliers_tot[top_outliers_tot$total_reads >= 5000 & top_outliers_tot$total_reads <= 75000, ]
if (nrow(outliers_p2_tot) > 0) {
    outliers_p2_tot$label_text <- paste0(outliers_p2_tot$circ_id, " (", outliers_p2_tot$gene_name, ")")
}

build_panel1_totals <- function(bar_color) {
    ggplot(plot_data_tot, aes(x=reads, y=circ_count)) +
        geom_segment(aes(xend=reads, y=0, yend=circ_count), color=bar_color, linewidth=1.1, lineend="square") +
        scale_y_log10(limits=c(1, 1e5), breaks=y_ticks, labels=comma) +
        coord_cartesian(xlim=c(0, 250)) +
        labs(x="Number of back-spliced reads", y="Number of circular RNAs",
             title="Distribution of circRNA Expression by Back-spliced Read Support (Totals)") +
        theme_classic(base_size=13) +
        theme(plot.title=element_text(size=12, face="bold"))
}

build_panel2_totals <- function(bar_color) {
    plot_data_tot_p2 <- plot_data_tot[plot_data_tot$reads >= 5000, ]

    p <- ggplot(plot_data_tot_p2, aes(x=reads, y=circ_count)) +
        geom_segment(aes(xend=reads, y=0, yend=circ_count), color=bar_color, linewidth=1.1, lineend="square") +
        scale_y_log10(limits=c(1, 1e5)) +
        scale_x_continuous(breaks=c(5000, 35000, 70000), labels=c("5k", "35k", "70k")) +
        coord_cartesian(xlim=c(5000, 75000), clip="off") +
        labs(x=NULL, y=NULL) +
        theme_classic(base_size=13) +
        theme(axis.line.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
              plot.margin=margin(t=5.5, r=35, b=5.5, l=5.5))

    if (nrow(outliers_p2_tot) > 0) {
        p <- p +
            geom_segment(data=outliers_p2_tot, aes(x=total_reads, xend=total_reads, y=3, yend=25000),
                         inherit.aes=FALSE, linewidth=0.5, linetype="dashed", color="black") +
            geom_text(data=outliers_p2_tot, aes(x=total_reads, y=30000, label=label_text),
                      inherit.aes=FALSE, angle=90, hjust=0, size=1.8, fontface="bold", color="black")
    }
    p
}

panel1_targets <- df_targets[df_targets$total_reads <= 250, ]
panel2_targets <- df_targets[df_targets$total_reads >= 5000 & df_targets$total_reads <= 75000, ]

print("Generating unannotated cohort total split-axis distribution plot...")
plot_totals_unannotated <- build_panel1_totals("red4") + build_panel2_totals("red4") + plot_layout(widths=c(0.76, 0.24))
save_ggplot(plot_totals_unannotated, "circ_abundance_distribution_unannotated", width_in=12.5)

print("Generating annotated cohort total split-axis distribution plot...")

panel1_targets <- df_targets[df_targets$total_reads <= 250, ]
panel1_targets <- panel1_targets[order(panel1_targets$total_reads), ]
if (nrow(panel1_targets) > 0) {
    panel1_targets$y_base <- as.numeric(counts_table_tot[as.character(panel1_targets$total_reads)])
    idx <- seq_len(nrow(panel1_targets))
    stagger_factor <- ((idx - 1) %% 4) + 1
    panel1_targets$y_label_pos <- panel1_targets$y_base * (1.8 + (stagger_factor * 1.5))
    panel1_targets$rotation_angle <- c(45, 60, 75)[((idx - 1) %% 3) + 1]
    panel1_targets$y_arrow_end <- panel1_targets$y_base * 1.8
    panel1_targets$y_arrow_start <- panel1_targets$y_label_pos * 0.9
    panel1_targets$label_text <- paste0(
        "circ", regmatches(panel1_targets$gene_name, regexpr(paste(target_genes, collapse="|"), panel1_targets$gene_name)),
        " (", panel1_targets$total_reads, ")"
    )
}

panel2_targets <- df_targets[df_targets$total_reads >= 5000 & df_targets$total_reads <= 80000, ]
panel2_targets <- panel2_targets[order(panel2_targets$total_reads), ]
if (nrow(panel2_targets) > 0) {
    panel2_targets$y_base <- as.numeric(counts_table_tot[as.character(panel2_targets$total_reads)])
    panel2_targets$y_base[is.na(panel2_targets$y_base)] <- 1
    idx <- seq_len(nrow(panel2_targets))
    stagger_factor <- ((idx - 1) %% 4) + 1
    panel2_targets$y_label_pos <- 400 * (1.5 ^ stagger_factor)
    panel2_targets$rotation_angle <- c(35, 50, 65)[((idx - 1) %% 3) + 1]
    panel2_targets$y_arrow_end <- panel2_targets$y_base * 2.5
    panel2_targets$y_arrow_start <- panel2_targets$y_label_pos * 0.85
    panel2_targets$label_text <- paste0(
        "circ", regmatches(panel2_targets$gene_name, regexpr(paste(target_genes, collapse="|"), panel2_targets$gene_name)),
        " (", format(panel2_targets$total_reads, big.mark=","), ")"
    )
}

cat_labels <- c("ALS" = "ALS-associated", "SYN" = "Synapse-associated")
fill_scale  <- scale_fill_manual(values=cat_colors, guide="none")
color_scale <- scale_color_manual(values=cat_colors, guide="none")

p1_annot <- build_panel1_totals("gray85") +
    labs(subtitle="Red = ALS-associated genes   |   Blue = Synapse-associated genes")
if (nrow(panel1_targets) > 0) {
    p1_annot <- p1_annot +
        geom_segment(data=panel1_targets, aes(x=total_reads, xend=total_reads, y=y_arrow_start, yend=y_arrow_end, color=cat),
                     inherit.aes=FALSE, linewidth=0.4, arrow=arrow(length=unit(0.06, "cm"), angle=15)) +
        geom_point(data=panel1_targets, aes(x=total_reads, y=y_base, fill=cat),
                   inherit.aes=FALSE, shape=21, color="black", size=2, stroke=0.3) +
        geom_text(data=panel1_targets, aes(x=total_reads, y=y_label_pos, label=label_text, angle=rotation_angle, color=cat),
                  inherit.aes=FALSE, hjust=0, size=1.9, fontface="italic") +
        fill_scale + color_scale
}

p2_annot <- build_panel2_totals("gray85")
if (nrow(panel2_targets) > 0) {
    p2_annot <- p2_annot +
        geom_segment(data=panel2_targets, aes(x=total_reads, xend=total_reads, y=y_arrow_start, yend=y_arrow_end, color=cat),
                     inherit.aes=FALSE, linewidth=0.4, arrow=arrow(length=unit(0.06, "cm"), angle=15)) +
        geom_point(data=panel2_targets, aes(x=total_reads, y=y_base, fill=cat),
                   inherit.aes=FALSE, shape=21, color="black", size=2, stroke=0.3) +
        geom_text(data=panel2_targets, aes(x=total_reads, y=y_label_pos, label=label_text, angle=rotation_angle, color=cat),
                  inherit.aes=FALSE, hjust=0, size=1.9, fontface="italic") +
        fill_scale + color_scale
}

plot_totals_annotated <- p1_annot + p2_annot + plot_layout(widths=c(0.76, 0.24))
save_ggplot(plot_totals_annotated, "circ_abundance_distribution", width_in=12.5)

print("Generating cohort average split-axis distribution plot...")

bin_width <- 0.05
breaks_seq <- seq(0, max(df_metrics$avg_reads) + bin_width, by=bin_width)
h <- hist(df_metrics$avg_reads, breaks=breaks_seq, plot=FALSE)

plot_data_avg <- data.frame(
    reads_bin_center = h$mids[h$counts > 0],
    circ_count = h$counts[h$counts > 0]
)

df_avg_sorted <- df_metrics[order(-df_metrics$avg_reads), ]
top_outliers_avg <- head(df_avg_sorted, 5)

max_avg <- max(plot_data_avg$reads_bin_center)
xlim_p2_avg <- c(max_avg * 0.15, max_avg * 1.05)

outliers_p2_avg <- top_outliers_avg[
    top_outliers_avg$avg_reads >= xlim_p2_avg[1] &
    top_outliers_avg$avg_reads <= xlim_p2_avg[2],
]
if (nrow(outliers_p2_avg) > 0) {
    outliers_p2_avg$label_text <- paste0(
        outliers_p2_avg$circ_id, " (", outliers_p2_avg$gene_name, ")"
    )
}

p1_avg <- ggplot(plot_data_avg, aes(x=reads_bin_center, y=circ_count)) +
    geom_segment(aes(xend=reads_bin_center, y=0, yend=circ_count), color="darkorange2", linewidth=1, lineend="square") +
    scale_y_log10(limits=c(1, 1e5), breaks=y_ticks, labels=comma) +
    coord_cartesian(xlim=c(0, 20)) +
    labs(x="Average number of back-spliced reads per sample", y="Number of circular RNAs",
         title="Distribution of circRNA Expression by Average Back-spliced Read Support") +
    theme_classic(base_size=13) +
    theme(plot.title=element_text(size=12, face="bold"))

p2_avg <- ggplot(plot_data_avg, aes(x=reads_bin_center, y=circ_count)) +
    geom_segment(aes(xend=reads_bin_center, y=0, yend=circ_count), color="darkorange2", linewidth=1, lineend="square") +
    scale_y_log10(limits=c(1, 1e5)) +
    scale_x_continuous(breaks=round(seq(xlim_p2_avg[1], xlim_p2_avg[2], length.out=3), 1)) +
    coord_cartesian(xlim=xlim_p2_avg) +
    labs(x=NULL, y=NULL) +
    theme_classic(base_size=13) +
    theme(axis.line.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

if (nrow(outliers_p2_avg) > 0) {
    p2_avg <- p2_avg +
        geom_segment(data=outliers_p2_avg, aes(x=avg_reads, xend=avg_reads, y=3, yend=25000),
                     inherit.aes=FALSE, linewidth=0.5, linetype="dashed", color="black") +
        geom_text(data=outliers_p2_avg, aes(x=avg_reads, y=30000, label=label_text),
                  inherit.aes=FALSE, angle=45, hjust=0, size=1.8, fontface="bold", color="black")
}

plot_averages <- p1_avg + p2_avg + plot_layout(widths=c(0.80, 0.20))
save_ggplot(plot_averages, "circ_abundance_averages")

print("SUCCESS: circ matrix rebuilt; 3 plots (unannotated totals, annotated totals, averages) saved as PNG and SVG via ggplot2/ggsave.")
EOF
