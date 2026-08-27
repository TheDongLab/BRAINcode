#!/bin/bash
#SBATCH --job-name=run_coloc_plot
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/run_coloc_plot.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/run_coloc_plot.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G

# ============================================================
# PARAMETERS
# ============================================================

TYPE="${1:-}"                         # eQTL, sQTL, or cQTL

ROOT="/home/zw529/donglab/data/target_ALS"
MR_DIR="${ROOT}/MR"
METADATA="${ROOT}/targetALS_rnaseq_metadata.csv"

FLANK_FRAC=0.20                       # 20% of trait width on EACH side
FLANK_MIN=2000                        # minimum flank = 2 kb
FLANK_MAX=20000                       # maximum flank = 20 kb

BIGWIG_BIN_BP=25                      # desired resolution
MAX_TRACK_POINTS=2500                 # prevent enormous plots
MAX_HITS=0                            # 0 = all passing colocs
TISSUE_FILTER=""                      # blank = all tissues
OVERWRITE=0                           # 1 = remake existing plots
DPI=180

# ============================================================

set -euo pipefail

case "$TYPE" in
    eQTL|sQTL|cQTL) ;;
    *)
        echo "Usage: sbatch $0 {eQTL|sQTL|cQTL}"
        exit 1
        ;;
esac

COLOC="${MR_DIR}/${TYPE}_SuSiE_coloc_summary.tsv"
OUTDIR="${MR_DIR}/coloc_raw_plots/${TYPE}"

mkdir -p "$OUTDIR"

if [[ ! -f "$COLOC" ]]; then
    echo "ERROR: coloc summary not found:"
    echo "  $COLOC"
    exit 1
fi

if [[ ! -f "$METADATA" ]]; then
    echo "ERROR: metadata not found:"
    echo "  $METADATA"
    exit 1
fi

# deepTools environments normally provide pyBigWig, pysam,
# numpy and matplotlib.
module load deepTools 2>/dev/null || true

export TYPE ROOT MR_DIR METADATA COLOC OUTDIR
export FLANK_FRAC FLANK_MIN FLANK_MAX
export BIGWIG_BIN_BP MAX_TRACK_POINTS MAX_HITS
export TISSUE_FILTER OVERWRITE DPI

python - <<'PY'
import os
import re
import csv
import sys
import glob
import math
import warnings
from pathlib import Path
from collections import defaultdict

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pyBigWig
    import pysam
except ImportError as e:
    sys.exit(
        f"\nERROR: missing Python package: {e}\n"
        "Load an environment containing matplotlib, numpy, pyBigWig and pysam.\n"
    )

warnings.filterwarnings("ignore")

TYPE = os.environ["TYPE"]
ROOT = Path(os.environ["ROOT"])
METADATA = Path(os.environ["METADATA"])
COLOC = Path(os.environ["COLOC"])
OUTDIR = Path(os.environ["OUTDIR"])

FLANK_FRAC = float(os.environ["FLANK_FRAC"])
FLANK_MIN = int(os.environ["FLANK_MIN"])
FLANK_MAX = int(os.environ["FLANK_MAX"])
BIGWIG_BIN_BP = int(os.environ["BIGWIG_BIN_BP"])
MAX_TRACK_POINTS = int(os.environ["MAX_TRACK_POINTS"])
MAX_HITS = int(os.environ["MAX_HITS"])
TISSUE_FILTER = os.environ["TISSUE_FILTER"].strip()
OVERWRITE = int(os.environ["OVERWRITE"])
DPI = int(os.environ["DPI"])


# ============================================================
# GENERAL HELPERS
# ============================================================

def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def clean_id(s):
    return str(s).strip()


def safe_name(s):
    s = str(s)
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", s)
    return s[:180]


def missing(v):
    return v is None or str(v).strip().lower() in {
        "", "na", "nan", "none", ".", "null"
    }


def truthy(v):
    return str(v).strip().lower() in {
        "true", "t", "1", "yes", "y", "pass", "passed"
    }


def split_line(line, delimiter):
    line = line.rstrip("\r\n")
    if delimiter == "\t":
        return line.split("\t")
    return line.split()


def detect_delimiter(first_line):
    return "\t" if "\t" in first_line else None


def find_header(headers, aliases, required=True):
    lookup = {norm(h): h for h in headers}

    for alias in aliases:
        if norm(alias) in lookup:
            return lookup[norm(alias)]

    for h in headers:
        nh = norm(h)
        for alias in aliases:
            na = norm(alias)
            if na and na in nh:
                return h

    if required:
        raise KeyError(
            f"Could not find column matching {aliases}. "
            f"Available columns: {headers}"
        )
    return None


def tissue_equal(a, b):
    return norm(a) == norm(b)


def strip_gene_version(x):
    return re.sub(r"\.\d+$", "", str(x))


# ============================================================
# READ COLOC SUMMARY
# ============================================================

with open(COLOC, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    summary_headers = reader.fieldnames or []

    if "H4_reference_pass" not in summary_headers:
        sys.exit(
            "ERROR: summary does not contain H4_reference_pass.\n"
            f"Columns found:\n{summary_headers}"
        )

    hits = [
        row for row in reader
        if truthy(row.get("H4_reference_pass"))
    ]

print("============================================================")
print(f"QTL TYPE          : {TYPE}")
print(f"COLOC SUMMARY     : {COLOC}")
print(f"H4 PASSING HITS   : {len(hits)}")
print("FILTER            : H4_reference_pass == TRUE")
print("============================================================")

if TISSUE_FILTER:
    hits = [
        h for h in hits
        if tissue_equal(h.get("tissue", ""), TISSUE_FILTER)
    ]
    print(f"After tissue filter ({TISSUE_FILTER}): {len(hits)}")

if MAX_HITS > 0:
    hits = hits[:MAX_HITS]
    print(f"After MAX_HITS={MAX_HITS}: {len(hits)}")

if not hits:
    sys.exit("No passing coloc hits to plot.")


# ============================================================
# DETECT EVENT COLUMN
# ============================================================

if TYPE == "eQTL":
    event_col = find_header(
        summary_headers,
        ["geneid", "gene_id", "gene"]
    )
    symbol_col = find_header(
        summary_headers,
        ["gene_symbol", "symbol"],
        required=False
    )

elif TYPE == "sQTL":
    event_col = find_header(
        summary_headers,
        [
            "junction_id", "junctionid",
            "splicing_id", "splice_id",
            "event_id", "phenotype_id"
        ]
    )
    symbol_col = find_header(
        summary_headers,
        ["gene_symbol", "symbol", "gene"],
        required=False
    )

else:
    event_col = find_header(
        summary_headers,
        [
            "circ_id", "circid",
            "circrna_id", "event_id",
            "phenotype_id"
        ]
    )
    symbol_col = find_header(
        summary_headers,
        ["gene_symbol", "symbol", "gene"],
        required=False
    )

tissue_col = find_header(summary_headers, ["tissue"])


# ============================================================
# METADATA
# ============================================================

with open(METADATA, newline="") as fh:
    meta_reader = csv.DictReader(fh)
    meta_headers = meta_reader.fieldnames or []
    metadata = list(meta_reader)

subject_col = find_header(
    meta_headers,
    ["externalsubjectid", "subject_id", "subjectid"]
)

sample_col = find_header(
    meta_headers,
    ["externalsampleid", "sample_id", "sampleid"]
)

meta_tissue_col = find_header(
    meta_headers,
    ["tissue"],
    required=False
)

rin_col = find_header(
    meta_headers,
    ["rin", "rna_integrity_number"],
    required=False
)


def get_rin(row):
    if rin_col is None:
        return -1.0
    try:
        return float(row.get(rin_col, ""))
    except Exception:
        return -1.0


def metadata_for_qtl_id(qtl_id, tissue):
    """
    QTL matrices may use either subject IDs or sample IDs.

    If qtl_id is already a sample ID, use that exact sample.
    If qtl_id is a subject ID, choose the highest-RIN sample from
    the requested tissue, matching QTL preprocessing behavior.
    """
    qtl_id = str(qtl_id)

    tissue_rows = metadata

    if meta_tissue_col:
        exact_tissue = [
            r for r in metadata
            if tissue_equal(r.get(meta_tissue_col, ""), tissue)
        ]
        if exact_tissue:
            tissue_rows = exact_tissue

    exact_sample = [
        r for r in tissue_rows
        if str(r.get(sample_col, "")).strip() == qtl_id
    ]

    if exact_sample:
        return max(exact_sample, key=get_rin)

    subject_rows = [
        r for r in tissue_rows
        if str(r.get(subject_col, "")).strip() == qtl_id
    ]

    if subject_rows:
        return max(subject_rows, key=get_rin)

    return None


# ============================================================
# MATRIX READER
# Reads only requested rows rather than loading a 200 MB+
# MatrixEQTL matrix into memory.
# ============================================================

def extract_matrix_rows(path, wanted):
    path = Path(path)

    wanted = set(str(x) for x in wanted)
    wanted_nover = {strip_gene_version(x) for x in wanted}

    result = {}

    with open(path) as fh:
        first = fh.readline()
        delim = detect_delimiter(first)
        header = split_line(first, delim)

        samples = header[1:]

        for line in fh:
            if not line.strip():
                continue

            fields = split_line(line, delim)

            if len(fields) < 2:
                continue

            rid = fields[0].strip()

            if (
                rid not in wanted
                and strip_gene_version(rid) not in wanted_nover
            ):
                continue

            vals = fields[1:]

            if len(vals) < len(samples):
                vals += [""] * (len(samples) - len(vals))

            result[rid] = dict(zip(samples, vals[:len(samples)]))

    return result


def find_matrix_row(rows, event):
    if event in rows:
        return rows[event]

    target = strip_gene_version(event)

    for rid, values in rows.items():
        if strip_gene_version(rid) == target:
            return values

    return None


# ============================================================
# LOCATION TABLE
# ============================================================

def load_location_table(path, kind):
    path = Path(path)

    with open(path) as fh:
        first = fh.readline()
        delim = detect_delimiter(first)
        fields = split_line(first, delim)

        # Detect whether first row is a header.
        is_data = (
            len(fields) >= 3
            and re.match(r"^(chr)?[0-9XYM]+$", fields[1], re.I)
            and re.match(r"^\d+$", fields[2])
        )

        if is_data:
            n = len(fields)
            headers = ["id", "chr", "start"]

            if n >= 4:
                headers.append("end")
            if n >= 5:
                headers.append("strand")
            while len(headers) < n:
                headers.append(f"extra{len(headers)}")

            all_lines = [fields]

        else:
            headers = fields
            all_lines = []

        for line in fh:
            if line.strip():
                all_lines.append(split_line(line, delim))

    rows = []

    for fields in all_lines:
        if len(fields) < len(headers):
            fields += [""] * (len(headers) - len(fields))

        rows.append(dict(zip(headers, fields)))

    if kind == "snp":
        id_col = find_header(
            headers, ["snpid", "snp_id", "rsid", "id"]
        )
        chr_col = find_header(
            headers, ["chr", "chrom", "chromosome"]
        )
        pos_col = find_header(
            headers, ["pos", "position", "start"]
        )

        return {
            str(r[id_col]): (
                str(r[chr_col]),
                int(float(r[pos_col]))
            )
            for r in rows
            if not missing(r.get(id_col))
            and not missing(r.get(chr_col))
            and not missing(r.get(pos_col))
        }

    if TYPE == "eQTL":
        id_aliases = ["geneid", "gene_id", "gene", "id"]
    elif TYPE == "sQTL":
        id_aliases = [
            "junction_id", "junctionid",
            "splicing_id", "event_id",
            "phenotype_id", "id"
        ]
    else:
        id_aliases = [
            "circ_id", "circid", "circrna_id",
            "event_id", "phenotype_id", "id"
        ]

    id_col = find_header(headers, id_aliases)
    chr_col = find_header(
        headers, ["chr", "chrom", "chromosome"]
    )
    start_col = find_header(
        headers,
        ["start", "left", "begin", "donor"]
    )
    end_col = find_header(
        headers,
        ["end", "right", "stop", "acceptor"]
    )
    strand_col = find_header(
        headers,
        ["strand"],
        required=False
    )

    out = {}

    for r in rows:
        try:
            rid = str(r[id_col])
            out[rid] = {
                "chr": str(r[chr_col]),
                "start": int(float(r[start_col])),
                "end": int(float(r[end_col])),
                "strand": (
                    str(r[strand_col])
                    if strand_col and not missing(r.get(strand_col))
                    else None
                )
            }
        except Exception:
            continue

    return out


def parse_event_coordinates(event):
    """
    Fallback for IDs such as:
      chr16:72108382-72108443:+
      chr2:24095229:24095962:+
    """
    s = str(event)

    patterns = [
        r"(chr[^:]+):(\d+)-(\d+):([+-])",
        r"(chr[^:]+):(\d+):(\d+):([+-])",
        r"(chr[^:]+):(\d+)-(\d+)",
        r"(chr[^:]+):(\d+):(\d+)"
    ]

    for pat in patterns:
        m = re.search(pat, s)

        if m:
            return {
                "chr": m.group(1),
                "start": int(m.group(2)),
                "end": int(m.group(3)),
                "strand": (
                    m.group(4)
                    if len(m.groups()) >= 4 else None
                )
            }

    return None


def find_location(locations, event):
    if event in locations:
        return locations[event]

    target = strip_gene_version(event)

    for rid, loc in locations.items():
        if strip_gene_version(rid) == target:
            return loc

    return parse_event_coordinates(event)


# ============================================================
# SAMPLE DIRECTORY RESOLUTION
# ============================================================

processed_cache = {}


def resolve_sample_dir(tissue, sample_id):
    processed = ROOT / tissue / "RNAseq" / "Processed"

    candidates = [
        sample_id,
        sample_id.replace("-", "_"),
        sample_id.replace("_", "-")
    ]

    for candidate in candidates:
        p = processed / candidate
        if p.is_dir():
            return p

    key = str(processed)

    if key not in processed_cache:
        lookup = {}

        if processed.is_dir():
            for p in processed.iterdir():
                if p.is_dir():
                    lookup[norm(p.name)] = p

        processed_cache[key] = lookup

    return processed_cache[key].get(norm(sample_id))


# ============================================================
# BIGWIG TRACKS
# ============================================================

def bw_chrom(bw, chrom):
    chroms = bw.chroms()

    if chrom in chroms:
        return chrom

    if chrom.startswith("chr") and chrom[3:] in chroms:
        return chrom[3:]

    test = "chr" + chrom

    if test in chroms:
        return test

    return None


def bigwig_vector(sample_dir, chrom, start, end, strand, n_bins):
    plus = sample_dir / \
        "STAR.Aligned.sortedByCoord.out.plus.normalized.bw"

    minus = sample_dir / \
        "STAR.Aligned.sortedByCoord.out.minus.normalized.bw"

    def read_one(path):
        if not path.exists():
            return None

        bw = pyBigWig.open(str(path))

        try:
            c = bw_chrom(bw, chrom)

            if c is None:
                return None

            vals = bw.stats(
                c,
                int(start),
                int(end),
                nBins=int(n_bins),
                type="mean"
            )

            vals = np.asarray(
                [0.0 if v is None else float(v) for v in vals],
                dtype=float
            )

            # Minus-strand normalized bigWigs are sometimes negative.
            return np.abs(vals)

        finally:
            bw.close()

    if strand == "+":
        return read_one(plus), str(plus)

    if strand == "-":
        return read_one(minus), str(minus)

    a = read_one(plus)
    b = read_one(minus)

    if a is None and b is None:
        return None, None

    if a is None:
        return b, str(minus)

    if b is None:
        return a, str(plus)

    return a + b, f"{plus};{minus}"


# ============================================================
# cQTL REMAPPED BAM TRACKS
# ============================================================

bam_mapped_cache = {}


def ensure_bam_index(path):
    path = str(path)

    bai1 = path + ".bai"
    bai2 = re.sub(r"\.bam$", ".bai", path)

    if os.path.exists(bai1) or os.path.exists(bai2):
        return True

    try:
        print(f"    Indexing BAM: {path}")
        pysam.index(path)
        return True
    except Exception as e:
        print(f"    WARNING: could not index {path}: {e}")
        return False


def resolve_bam_chrom(bam, chrom):
    refs = set(bam.references)

    if chrom in refs:
        return chrom

    if chrom.startswith("chr") and chrom[3:] in refs:
        return chrom[3:]

    test = "chr" + chrom

    if test in refs:
        return test

    return None


def mapped_reads_from_star(sample_dir):
    star = sample_dir / "STAR.Aligned.sortedByCoord.out.bam"

    if not star.exists():
        return None

    key = str(star)

    if key in bam_mapped_cache:
        return bam_mapped_cache[key]

    try:
        if not ensure_bam_index(star):
            return None

        with pysam.AlignmentFile(str(star), "rb") as bam:
            n = float(bam.mapped)

        bam_mapped_cache[key] = n
        return n

    except Exception:
        return None


def bin_array(values, n_bins):
    values = np.asarray(values, dtype=float)

    if len(values) == n_bins:
        return values

    edges = np.linspace(
        0, len(values), n_bins + 1
    ).astype(int)

    out = np.zeros(n_bins, dtype=float)

    for i in range(n_bins):
        a, b = edges[i], edges[i + 1]

        if b > a:
            out[i] = np.mean(values[a:b])

    return out


def remap_bam_vector(sample_dir, chrom, start, end, n_bins):
    bams = sorted(
        sample_dir.glob("*remap.sorted.bam")
    )

    if not bams:
        return None, None

    bam_path = bams[0]

    if not ensure_bam_index(bam_path):
        return None, str(bam_path)

    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            c = resolve_bam_chrom(bam, chrom)

            if c is None:
                return None, str(bam_path)

            cov = bam.count_coverage(
                c,
                int(start),
                int(end),
                quality_threshold=0
            )

            cov = np.sum(
                np.asarray(cov, dtype=float),
                axis=0
            )

        # Normalize remapped regional signal by total mapped
        # reads in the ordinary RNA-seq STAR BAM.
        denominator = mapped_reads_from_star(sample_dir)

        if denominator is None or denominator <= 0:
            print(
                f"    WARNING: STAR mapped-read denominator "
                f"unavailable for {sample_dir.name}"
            )
            return None, str(bam_path)

        rpm = cov * 1_000_000.0 / denominator
        rpm = bin_array(rpm, n_bins)

        return rpm, str(bam_path)

    except Exception as e:
        print(
            f"    WARNING: BAM extraction failed "
            f"for {bam_path}: {e}"
        )
        return None, str(bam_path)


# ============================================================
# SNP SELECTION
# ============================================================

def possible_snps(hit):
    candidates = []

    # Main choice.
    for col in [
        "best_candidate_SNP",
        "best_candidate_snp"
    ]:
        if col in hit and not missing(hit[col]):
            candidates.append(hit[col])

    # QTL-specific top SNP fallback.
    prefixes = {
        "eQTL": ["top_eQTL_SNP", "top_QTL_SNP"],
        "sQTL": ["top_sQTL_SNP", "top_QTL_SNP"],
        "cQTL": ["top_cQTL_SNP", "top_QTL_SNP"]
    }

    for col in prefixes[TYPE]:
        if col in hit and not missing(hit[col]):
            candidates.append(hit[col])

    # Other non-GWAS top SNP fields.
    for k, v in hit.items():
        nk = norm(k)

        if (
            "topsnp" in nk
            and "gwas" not in nk
            and not missing(v)
        ):
            candidates.append(v)

    # Finally try all candidate SNPs.
    if (
        "candidate_SNPs" in hit
        and not missing(hit["candidate_SNPs"])
    ):
        candidates.extend(
            x for x in hit["candidate_SNPs"].split(";")
            if x
        )

    out = []
    seen = set()

    for x in candidates:
        x = str(x).strip()

        if x and x not in seen:
            seen.add(x)
            out.append(x)

    return out


# ============================================================
# GENOTYPE GROUP
# ============================================================

def genotype_group(value):
    try:
        x = float(value)
    except Exception:
        return None

    if not np.isfinite(x) or x < -0.1 or x > 2.1:
        return None

    # MatrixEQTL genotype files from PLINK should usually be 0/1/2.
    # Round minor floating representation differences.
    return int(np.clip(np.rint(x), 0, 2))


# ============================================================
# PROCESS TISSUE BY TISSUE
# ============================================================

hits_by_tissue = defaultdict(list)

for hit in hits:
    tissue = hit.get(tissue_col, "").strip()
    hits_by_tissue[tissue].append(hit)


total_plotted = 0
total_skipped = 0

for tissue, tissue_hits in hits_by_tissue.items():

    print()
    print("============================================================")
    print(f"TISSUE: {tissue}")
    print(f"PASSING HITS: {len(tissue_hits)}")
    print("============================================================")

    qtl_dir = ROOT / tissue / TYPE

    if TYPE == "eQTL":
        phenotype_file = qtl_dir / f"expression_{tissue}.txt"
        trait_location_file = qtl_dir / "gene_location.txt"

    elif TYPE == "sQTL":
        phenotype_file = qtl_dir / f"splicing_{tissue}.txt"
        trait_location_file = qtl_dir / "splicing_location.txt"

    else:
        phenotype_file = qtl_dir / f"circ_{tissue}.txt"
        trait_location_file = qtl_dir / "circ_location.txt"

    snp_file = qtl_dir / f"snp_{tissue}.txt"
    snp_location_file = qtl_dir / "snp_location.txt"

    required = [
        phenotype_file,
        trait_location_file,
        snp_file
    ]

    missing_files = [p for p in required if not p.exists()]

    if missing_files:
        print("SKIPPING TISSUE: missing input files:")
        for p in missing_files:
            print(f"  {p}")
        total_skipped += len(tissue_hits)
        continue

    events = {
        h[event_col].strip()
        for h in tissue_hits
        if not missing(h.get(event_col))
    }

    all_snps = set()

    for h in tissue_hits:
        all_snps.update(possible_snps(h))

    print("Reading requested phenotype rows...")
    phenotype_rows = extract_matrix_rows(
        phenotype_file,
        events
    )

    print("Reading requested genotype rows...")
    genotype_rows = extract_matrix_rows(
        snp_file,
        all_snps
    )

    print("Reading trait locations...")
    trait_locations = load_location_table(
        trait_location_file,
        "trait"
    )

    if snp_location_file.exists():
        snp_locations = load_location_table(
            snp_location_file,
            "snp"
        )
    else:
        snp_locations = {}

    for index, hit in enumerate(tissue_hits, 1):

        event = hit[event_col].strip()

        symbol = (
            hit.get(symbol_col, "").strip()
            if symbol_col else ""
        )

        pp_h4 = hit.get("max_PP_H4", "")

        print()
        print(
            f"[{index}/{len(tissue_hits)}] "
            f"{tissue} | {event} | {symbol}"
        )

        # ----------------------------------------------------
        # Find phenotype
        # ----------------------------------------------------

        phenotype = find_matrix_row(
            phenotype_rows,
            event
        )

        if phenotype is None:
            print(
                f"  SKIP: event not found in "
                f"{phenotype_file.name}"
            )
            total_skipped += 1
            continue

        # ----------------------------------------------------
        # Find usable SNP
        # ----------------------------------------------------

        snp = None
        genotype = None

        for candidate in possible_snps(hit):
            if candidate in genotype_rows:
                snp = candidate
                genotype = genotype_rows[candidate]
                break

        if genotype is None:
            print("  SKIP: no candidate SNP found in genotype matrix")
            print(
                "  Tried: "
                + ", ".join(possible_snps(hit)[:10])
            )
            total_skipped += 1
            continue

        print(f"  SNP: {snp}")

        # ----------------------------------------------------
        # Trait genomic location
        # ----------------------------------------------------

        loc = find_location(
            trait_locations,
            event
        )

        if loc is None:
            print(
                f"  SKIP: unable to determine coordinates "
                f"for {event}"
            )
            total_skipped += 1
            continue

        chrom = loc["chr"]
        trait_start = min(loc["start"], loc["end"])
        trait_end = max(loc["start"], loc["end"])
        strand = loc.get("strand")

        trait_width = max(
            1,
            trait_end - trait_start
        )

        flank = int(
            max(
                FLANK_MIN,
                math.ceil(trait_width * FLANK_FRAC)
            )
        )

        flank = min(
            flank,
            FLANK_MAX
        )

        plot_start = max(
            0,
            trait_start - flank
        )

        plot_end = trait_end + flank
        plot_width = plot_end - plot_start

        n_bins = int(
            min(
                MAX_TRACK_POINTS,
                max(
                    50,
                    math.ceil(
                        plot_width / BIGWIG_BIN_BP
                    )
                )
            )
        )

        x = np.linspace(
            plot_start,
            plot_end,
            n_bins,
            endpoint=False
        )

        x += (
            (plot_end - plot_start)
            / n_bins / 2
        )

        print(
            f"  LOCUS: {chrom}:{trait_start}-{trait_end} "
            f"strand={strand or 'unknown'}"
        )
        print(
            f"  WINDOW: {chrom}:{plot_start}-{plot_end} "
            f"(flank={flank:,} bp)"
        )

        # ----------------------------------------------------
        # SNP coordinate
        # ----------------------------------------------------

        snp_chr = None
        snp_pos = None

        if snp in snp_locations:
            snp_chr, snp_pos = snp_locations[snp]

            if norm(snp_chr) != norm(chrom):
                snp_pos = None

        # ----------------------------------------------------
        # Align genotype + phenotype subjects
        # ----------------------------------------------------

        common_ids = [
            qid for qid in phenotype
            if qid in genotype
        ]

        subject_records = []

        for qid in common_ids:
            g = genotype_group(
                genotype[qid]
            )

            if g is None:
                continue

            try:
                pheno = float(phenotype[qid])
            except Exception:
                continue

            if not np.isfinite(pheno):
                continue

            meta = metadata_for_qtl_id(
                qid,
                tissue
            )

            sample_id = None
            subject_id = None
            sample_dir = None

            if meta:
                sample_id = str(
                    meta.get(sample_col, "")
                ).strip()

                subject_id = str(
                    meta.get(subject_col, "")
                ).strip()

                if sample_id:
                    sample_dir = resolve_sample_dir(
                        tissue,
                        sample_id
                    )

            subject_records.append({
                "qtl_id": qid,
                "subject_id": subject_id,
                "sample_id": sample_id,
                "genotype": g,
                "genotype_raw": genotype[qid],
                "phenotype": pheno,
                "sample_dir": sample_dir,
                "track": None,
                "track_source": None
            })

        counts = {
            g: sum(
                r["genotype"] == g
                for r in subject_records
            )
            for g in [0, 1, 2]
        }

        print(
            "  GENOTYPES: "
            f"0={counts[0]}, "
            f"1={counts[1]}, "
            f"2={counts[2]}"
        )

        if not subject_records:
            print("  SKIP: no aligned subjects")
            total_skipped += 1
            continue

        # ----------------------------------------------------
        # Extract tracks
        # ----------------------------------------------------

        print("  Extracting RNA-seq tracks...")

        for r in subject_records:
            sample_dir = r["sample_dir"]

            if sample_dir is None:
                continue

            if TYPE in {"eQTL", "sQTL"}:
                vec, source = bigwig_vector(
                    sample_dir,
                    chrom,
                    plot_start,
                    plot_end,
                    strand,
                    n_bins
                )
            else:
                vec, source = remap_bam_vector(
                    sample_dir,
                    chrom,
                    plot_start,
                    plot_end,
                    n_bins
                )

            r["track"] = vec
            r["track_source"] = source

        tracks_by_group = {
            g: [
                r["track"]
                for r in subject_records
                if r["genotype"] == g
                and r["track"] is not None
            ]
            for g in [0, 1, 2]
        }

        track_counts = {
            g: len(tracks_by_group[g])
            for g in [0, 1, 2]
        }

        print(
            "  TRACKS: "
            f"0={track_counts[0]}, "
            f"1={track_counts[1]}, "
            f"2={track_counts[2]}"
        )

        # ----------------------------------------------------
        # Output names
        # ----------------------------------------------------

        label = symbol if symbol else event

        basename = safe_name(
            f"{tissue}__{TYPE}__{label}__{snp}"
        )

        pdf = OUTDIR / f"{basename}.pdf"
        png = OUTDIR / f"{basename}.png"
        subject_tsv = OUTDIR / f"{basename}.subjects.tsv"

        if (
            not OVERWRITE
            and pdf.exists()
            and png.exists()
        ):
            print("  EXISTS; skipping plot")
            continue

        # ----------------------------------------------------
        # Subject-level audit table
        # ----------------------------------------------------

        with open(subject_tsv, "w", newline="") as fh:
            writer = csv.writer(
                fh,
                delimiter="\t"
            )

            writer.writerow([
                "tissue",
                "qtl_type",
                "event",
                "gene_symbol",
                "snp",
                "qtl_id",
                "subject_id",
                "sample_id",
                "genotype_group",
                "genotype_raw",
                "phenotype",
                "track_available",
                "track_source"
            ])

            for r in subject_records:
                writer.writerow([
                    tissue,
                    TYPE,
                    event,
                    symbol,
                    snp,
                    r["qtl_id"],
                    r["subject_id"] or "",
                    r["sample_id"] or "",
                    r["genotype"],
                    r["genotype_raw"],
                    r["phenotype"],
                    r["track"] is not None,
                    r["track_source"] or ""
                ])

        # ----------------------------------------------------
        # GLOBAL y-axis ACROSS GENOTYPES
        # ----------------------------------------------------

        all_track_values = []

        for g in [0, 1, 2]:
            for arr in tracks_by_group[g]:
                all_track_values.append(
                    np.asarray(arr, dtype=float)
                )

        if all_track_values:
            ymax = max(
                float(np.nanmax(a))
                for a in all_track_values
                if len(a)
            )

            if not np.isfinite(ymax) or ymax <= 0:
                ymax = 1.0
        else:
            ymax = 1.0

        ymax *= 1.05

        # ----------------------------------------------------
        # PLOT
        # ----------------------------------------------------

        fig = plt.figure(
            figsize=(16, 8)
        )

        gs = fig.add_gridspec(
            2,
            3,
            height_ratios=[1.0, 1.8],
            hspace=0.35,
            wspace=0.08
        )

        # ---------------- RAW PHENOTYPE ----------------

        ax0 = fig.add_subplot(
            gs[0, :]
        )

        rng = np.random.default_rng(12345)

        for g in [0, 1, 2]:
            vals = np.asarray([
                r["phenotype"]
                for r in subject_records
                if r["genotype"] == g
            ])

            if not len(vals):
                continue

            jitter = rng.normal(
                0,
                0.055,
                size=len(vals)
            )

            ax0.scatter(
                np.full(len(vals), g) + jitter,
                vals,
                s=24,
                alpha=0.60
            )

            mean_val = np.mean(vals)

            ax0.plot(
                [g - 0.18, g + 0.18],
                [mean_val, mean_val],
                linewidth=3
            )

        ax0.set_xticks(
            [0, 1, 2],
            [
                f"0 copies (n={counts[0]})",
                f"1 copy (n={counts[1]})",
                f"2 copies (n={counts[2]})"
            ]
        )

        if TYPE == "eQTL":
            ax0.set_ylabel("Raw MatrixEQTL expression phenotype")
        elif TYPE == "sQTL":
            ax0.set_ylabel("Raw MatrixEQTL splicing phenotype")
        else:
            ax0.set_ylabel("Raw MatrixEQTL circRNA phenotype")

        title_event = (
            f"{symbol} ({event})"
            if symbol and symbol != event
            else event
        )

        title = (
            f"{tissue} | {TYPE} | {title_event} | {snp}"
        )

        if pp_h4 and not missing(pp_h4):
            title += f" | PP.H4={pp_h4}"

        ax0.set_title(
            title,
            fontsize=13
        )

        ax0.spines["top"].set_visible(False)
        ax0.spines["right"].set_visible(False)

        # ---------------- COVERAGE ----------------

        genotype_axes = []

        for g in [0, 1, 2]:
            ax = fig.add_subplot(
                gs[1, g]
            )
            genotype_axes.append(ax)

            group_tracks = tracks_by_group[g]

            for arr in group_tracks:
                ax.plot(
                    x,
                    arr,
                    linewidth=0.65,
                    alpha=0.18
                )

            if group_tracks:
                mean_track = np.nanmean(
                    np.vstack(group_tracks),
                    axis=0
                )

                ax.plot(
                    x,
                    mean_track,
                    linewidth=2.4,
                    label="Genotype mean"
                )
            else:
                ax.text(
                    0.5,
                    0.5,
                    "No usable tracks",
                    transform=ax.transAxes,
                    ha="center",
                    va="center"
                )

            # Trait itself
            ax.axvspan(
                trait_start,
                trait_end,
                alpha=0.12
            )

            # Coloc SNP
            if (
                snp_pos is not None
                and plot_start <= snp_pos <= plot_end
            ):
                ax.axvline(
                    snp_pos,
                    linestyle="--",
                    linewidth=1.3,
                    alpha=0.8
                )

            ax.set_xlim(
                plot_start,
                plot_end
            )

            # EXACT SAME Y SCALE FOR ALL GENOTYPES.
            ax.set_ylim(
                0,
                ymax
            )

            ax.set_title(
                f"Genotype {g} | "
                f"subjects n={counts[g]} | "
                f"tracks n={track_counts[g]}"
            )

            ax.set_xlabel(
                f"{chrom} genomic position"
            )

            if g == 0:
                if TYPE == "cQTL":
                    ax.set_ylabel(
                        "Remapped read coverage\n"
                        "(RPM; normalized to STAR mapped reads)"
                    )
                else:
                    ax.set_ylabel(
                        "Normalized RNA-seq read density"
                    )
            else:
                ax.tick_params(
                    labelleft=False
                )

            ax.ticklabel_format(
                axis="x",
                style="plain",
                useOffset=False
            )

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        fig.text(
            0.5,
            0.015,
            (
                f"Trait: {chrom}:{trait_start:,}-{trait_end:,} | "
                f"plot flank: {flank:,} bp each side | "
                f"shaded region = trait locus | "
                f"dashed line = {snp}"
            ),
            ha="center",
            fontsize=9
        )

        fig.savefig(
            pdf,
            bbox_inches="tight"
        )

        fig.savefig(
            png,
            dpi=DPI,
            bbox_inches="tight"
        )

        plt.close(fig)

        print(f"  PDF: {pdf}")
        print(f"  PNG: {png}")
        print(f"  TSV: {subject_tsv}")

        total_plotted += 1


print()
print("============================================================")
print("DONE")
print(f"QTL type        : {TYPE}")
print(f"Plots generated : {total_plotted}")
print(f"Skipped hits    : {total_skipped}")
print(f"Output          : {OUTDIR}")
print("============================================================")
PY
