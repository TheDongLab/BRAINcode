#!/bin/bash
#SBATCH --job-name=run_coloc_plot
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/run_coloc_plot.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/run_coloc_plot.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G

# ============================================================
# PARAMETERS
# ============================================================

TYPE="${1:-}"                         # eQTL, sQTL, cQTL
ROOT="/home/zw529/donglab/data/target_ALS"
REAL_DATA_ROOT="/nfs/roberts/pi/pi_xd96/data/target_ALS"
MR_DIR="${ROOT}/MR"
METADATA="${ROOT}/targetALS_rnaseq_metadata.csv"

FLANK_FRAC=0.20
FLANK_MIN=2000
FLANK_MAX=20000
BIGWIG_BIN_BP=25
MAX_TRACK_POINTS=2500
MAX_HITS=0                            # 0 = all passing colocs
TISSUE_FILTER=""                      # blank = all tissues
OVERWRITE=0
DPI=180

# ============================================================

set -euo pipefail

case "$TYPE" in
    eQTL|sQTL|cQTL) ;;
    *) echo "Usage: sbatch $0 {eQTL|sQTL|cQTL}"; exit 1 ;;
esac

COLOC="${MR_DIR}/${TYPE}_SuSiE_coloc_summary.tsv"
OUTDIR="${MR_DIR}/coloc_raw_plots/${TYPE}"
mkdir -p "$OUTDIR"

[[ -f "$COLOC" ]] || { echo "ERROR: coloc summary not found: $COLOC"; exit 1; }
[[ -f "$METADATA" ]] || { echo "ERROR: metadata not found: $METADATA"; exit 1; }

module load deepTools 2>/dev/null || true

export TYPE ROOT REAL_DATA_ROOT METADATA COLOC OUTDIR
export FLANK_FRAC FLANK_MIN FLANK_MAX BIGWIG_BIN_BP MAX_TRACK_POINTS
export MAX_HITS TISSUE_FILTER OVERWRITE DPI

python - <<'PY'
import os, re, csv, sys, math, warnings
from pathlib import Path
from collections import defaultdict

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    import pyBigWig, pysam
except ImportError as e:
    sys.exit(f"ERROR: missing Python package: {e}\nNeed numpy, matplotlib, pyBigWig and pysam.")

warnings.filterwarnings("ignore")

TYPE = os.environ["TYPE"]
ROOT = Path(os.environ["ROOT"])
REAL_DATA_ROOT = Path(os.environ["REAL_DATA_ROOT"])
METADATA, COLOC, OUTDIR = map(Path, [os.environ["METADATA"], os.environ["COLOC"], os.environ["OUTDIR"]])

FLANK_FRAC, FLANK_MIN, FLANK_MAX = float(os.environ["FLANK_FRAC"]), int(os.environ["FLANK_MIN"]), int(os.environ["FLANK_MAX"])
BIGWIG_BIN_BP, MAX_TRACK_POINTS = int(os.environ["BIGWIG_BIN_BP"]), int(os.environ["MAX_TRACK_POINTS"])
MAX_HITS, TISSUE_FILTER = int(os.environ["MAX_HITS"]), os.environ["TISSUE_FILTER"].strip()
OVERWRITE, DPI = int(os.environ["OVERWRITE"]), int(os.environ["DPI"])

GENOTYPE_LABELS = {0: "Ref/Ref", 1: "Het", 2: "Hom Alt"}
TOP_COLORS = {0: "#90EE90", 1: "#ADD8E6", 2: "#FFB6C1"}
TRACK_COLORS = {0: "#228B22", 1: "#4682B4", 2: "#C75B7A"}

Y_LABELS = {
    "eQTL": "Z-score normalized expression",
    "sQTL": "Percent spliced in",
    "cQTL": "Percent circularized"
}

# ============================================================
# HELPERS
# ============================================================

def norm(s): return re.sub(r"[^a-z0-9]", "", str(s).lower())
def safe_name(s): return re.sub(r"[^A-Za-z0-9._+-]+", "_", str(s))[:180]
def strip_gene_version(x): return re.sub(r"\.\d+$", "", str(x))
def tissue_equal(a, b): return norm(a) == norm(b)
def tissue_display(x): return str(x).replace("_", " ")
def missing(v): return v is None or str(v).strip().lower() in {"", "na", "nan", "none", ".", "null"}
def truthy(v): return str(v).strip().lower() in {"true", "t", "1", "yes", "y", "pass", "passed"}
def detect_delimiter(line): return "\t" if "\t" in line else None
def split_line(line, delim): return line.rstrip("\r\n").split("\t") if delim == "\t" else line.split()

def genomic_formatter(x, pos=None):
    if abs(x) >= 1_000_000: return f"{x / 1_000_000:.3f} Mb"
    if abs(x) >= 1_000: return f"{x / 1_000:.1f} kb"
    return f"{x:.0f}"

def find_header(headers, aliases, required=True):
    lookup = {norm(h): h for h in headers}
    for alias in aliases:
        if norm(alias) in lookup: return lookup[norm(alias)]
    for h in headers:
        if any(norm(a) and norm(a) in norm(h) for a in aliases): return h
    if required: raise KeyError(f"Could not find {aliases}. Available: {headers}")
    return None

# ============================================================
# COLOC SUMMARY
# ============================================================

with open(COLOC, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    summary_headers = reader.fieldnames or []
    if "H4_reference_pass" not in summary_headers: sys.exit("ERROR: summary lacks H4_reference_pass")
    hits = [r for r in reader if truthy(r.get("H4_reference_pass"))]

print("=" * 60)
print(f"QTL TYPE        : {TYPE}")
print(f"H4 PASSING HITS : {len(hits)}")
print("FILTER          : H4_reference_pass == TRUE")
print("=" * 60)

if TISSUE_FILTER:
    hits = [h for h in hits if tissue_equal(h.get("tissue", ""), TISSUE_FILTER)]
if MAX_HITS > 0:
    hits = hits[:MAX_HITS]
if not hits:
    sys.exit("No passing coloc hits to plot.")

tissue_col = find_header(summary_headers, ["tissue"])

if TYPE == "eQTL":
    event_col = find_header(summary_headers, ["geneid", "gene_id", "gene"])
    symbol_col = find_header(summary_headers, ["gene_symbol", "symbol"], required=False)
elif TYPE == "sQTL":
    event_col = find_header(summary_headers, ["junction_id", "junctionid", "splicing_id", "splice_id", "event_id", "phenotype_id"])
    symbol_col = find_header(summary_headers, ["gene_symbol", "symbol", "gene"], required=False)
else:
    event_col = find_header(summary_headers, ["circ_id", "circid", "circrna_id", "event_id", "phenotype_id"])
    symbol_col = find_header(summary_headers, ["gene_symbol", "symbol", "gene"], required=False)

# ============================================================
# METADATA
# ============================================================

with open(METADATA, newline="") as fh:
    reader = csv.DictReader(fh)
    meta_headers, metadata = reader.fieldnames or [], list(reader)

subject_col = find_header(meta_headers, ["externalsubjectid", "subject_id", "subjectid"])
sample_col = find_header(meta_headers, ["externalsampleid", "sample_id", "sampleid"])
meta_tissue_col = find_header(meta_headers, ["tissue"], required=False)
rin_col = find_header(meta_headers, ["rin", "rna_integrity_number"], required=False)

def get_rin(row):
    try: return float(row.get(rin_col, "")) if rin_col else -1.0
    except: return -1.0

def metadata_for_qtl_id(qtl_id, tissue):
    qtl_id = str(qtl_id)
    rows = metadata

    if meta_tissue_col:
        tissue_rows = [r for r in metadata if tissue_equal(r.get(meta_tissue_col, ""), tissue)]
        if tissue_rows: rows = tissue_rows

    exact_sample = [r for r in rows if str(r.get(sample_col, "")).strip() == qtl_id]
    if exact_sample: return max(exact_sample, key=get_rin)

    subject_rows = [r for r in rows if str(r.get(subject_col, "")).strip() == qtl_id]
    return max(subject_rows, key=get_rin) if subject_rows else None

# ============================================================
# MATRIX READER
# ============================================================

def extract_matrix_rows(path, wanted):
    wanted, result = set(map(str, wanted)), {}
    wanted_nover = {strip_gene_version(x) for x in wanted}

    with open(path) as fh:
        first = fh.readline()
        delim = detect_delimiter(first)
        samples = split_line(first, delim)[1:]

        for line in fh:
            if not line.strip(): continue
            fields = split_line(line, delim)
            if len(fields) < 2: continue

            rid = fields[0].strip()
            if rid not in wanted and strip_gene_version(rid) not in wanted_nover: continue

            vals = fields[1:] + [""] * max(0, len(samples) - len(fields[1:]))
            result[rid] = dict(zip(samples, vals[:len(samples)]))

    return result

def find_matrix_row(rows, event):
    if event in rows: return rows[event]
    target = strip_gene_version(event)
    return next((v for rid, v in rows.items() if strip_gene_version(rid) == target), None)

# ============================================================
# LOCATIONS
# ============================================================

def load_location_table(path, kind):
    with open(path) as fh:
        first = fh.readline()
        delim = detect_delimiter(first)
        fields = split_line(first, delim)

        is_data = (
            len(fields) >= 3
            and re.match(r"^(chr)?[0-9XYM]+$", fields[1], re.I)
            and re.match(r"^\d+$", fields[2])
        )

        if is_data:
            headers = ["id", "chr", "start"]
            if len(fields) >= 4: headers.append("end")
            if len(fields) >= 5: headers.append("strand")
            while len(headers) < len(fields): headers.append(f"extra{len(headers)}")
            lines = [fields]
        else:
            headers, lines = fields, []

        lines += [split_line(line, delim) for line in fh if line.strip()]

    rows = []
    for f in lines:
        f += [""] * max(0, len(headers) - len(f))
        rows.append(dict(zip(headers, f)))

    if kind == "snp":
        idc = find_header(headers, ["snpid", "snp_id", "rsid", "id"])
        cc = find_header(headers, ["chr", "chrom", "chromosome"])
        pc = find_header(headers, ["pos", "position", "start"])

        return {
            str(r[idc]): (str(r[cc]), int(float(r[pc])))
            for r in rows
            if not missing(r.get(idc)) and not missing(r.get(cc)) and not missing(r.get(pc))
        }

    aliases = {
        "eQTL": ["geneid", "gene_id", "gene", "id"],
        "sQTL": ["junction_id", "junctionid", "splicing_id", "event_id", "phenotype_id", "id"],
        "cQTL": ["circ_id", "circid", "circrna_id", "event_id", "phenotype_id", "id"]
    }

    idc = find_header(headers, aliases[TYPE])
    cc = find_header(headers, ["chr", "chrom", "chromosome"])
    sc = find_header(headers, ["start", "left", "begin", "donor"])
    ec = find_header(headers, ["end", "right", "stop", "acceptor"])
    strandc = find_header(headers, ["strand"], required=False)

    out = {}
    for r in rows:
        try:
            out[str(r[idc])] = {
                "chr": str(r[cc]), "start": int(float(r[sc])), "end": int(float(r[ec])),
                "strand": str(r[strandc]) if strandc and not missing(r.get(strandc)) else None
            }
        except:
            pass

    return out

def parse_event_coordinates(event):
    for pat in [
        r"(chr[^:]+):(\d+)-(\d+):([+-])",
        r"(chr[^:]+):(\d+):(\d+):([+-])",
        r"(chr[^:]+):(\d+)-(\d+)",
        r"(chr[^:]+):(\d+):(\d+)"
    ]:
        m = re.search(pat, str(event))
        if m:
            return {
                "chr": m.group(1), "start": int(m.group(2)), "end": int(m.group(3)),
                "strand": m.group(4) if len(m.groups()) >= 4 else None
            }
    return None

def find_location(locations, event):
    if event in locations: return locations[event]
    target = strip_gene_version(event)
    return next(
        (loc for rid, loc in locations.items() if strip_gene_version(rid) == target),
        parse_event_coordinates(event)
    )

# ============================================================
# SAMPLE DIRECTORIES
# ============================================================

processed_cache = {}

def resolve_sample_dir(tissue, sample_id):
    processed = ROOT / tissue / "RNAseq" / "Processed"

    for candidate in [sample_id, sample_id.replace("-", "_"), sample_id.replace("_", "-")]:
        p = processed / candidate
        if p.is_dir(): return p

    key = str(processed)
    if key not in processed_cache:
        processed_cache[key] = (
            {norm(p.name): p for p in processed.iterdir() if p.is_dir()}
            if processed.is_dir() else {}
        )

    return processed_cache[key].get(norm(sample_id))

# ============================================================
# BIGWIG — eQTL / sQTL
# ============================================================

def bw_chrom(bw, chrom):
    chroms = bw.chroms()

    candidates = [chrom]
    if chrom.startswith("chr"): candidates.append(chrom[3:])
    else: candidates.append("chr" + chrom)

    for c in candidates:
        if c in chroms: return c

    return None

def bigwig_vector(sample_dir, chrom, start, end, strand, n_bins):
    plus = sample_dir / "STAR.Aligned.sortedByCoord.out.plus.normalized.bw"
    minus = sample_dir / "STAR.Aligned.sortedByCoord.out.minus.normalized.bw"

    def read_one(path):
        if not path.exists(): return None

        bw = pyBigWig.open(str(path))
        try:
            c = bw_chrom(bw, chrom)
            if c is None: return None

            vals = bw.stats(c, int(start), int(end), nBins=int(n_bins), type="mean")
            return np.abs(np.asarray([0.0 if v is None else float(v) for v in vals]))
        finally:
            bw.close()

    if strand == "+": return read_one(plus), str(plus)
    if strand == "-": return read_one(minus), str(minus)

    a, b = read_one(plus), read_one(minus)

    if a is None and b is None: return None, None
    if a is None: return b, str(minus)
    if b is None: return a, str(plus)

    return a + b, f"{plus};{minus}"

# ============================================================
# REAL BAM PATHS — cQTL
# ============================================================

bam_mapped_cache = {}

def real_sample_dir(tissue, sample_id):
    return REAL_DATA_ROOT / tissue / "RNAseq" / "Processed" / sample_id

def real_remap_bam(tissue, sample_id):
    return (
        real_sample_dir(tissue, sample_id)
        / "remap_chimeric"
        / f"{sample_id}.remap.Aligned.sortedByCoord.out.bam"
    )

def real_star_bam(tissue, sample_id):
    return real_sample_dir(tissue, sample_id) / "STAR.Aligned.sortedByCoord.out.bam"

def ensure_bam_index(path):
    path = Path(path)

    if not path.exists():
        print(f"    WARNING: BAM does not exist: {path}")
        return False

    bai1 = Path(str(path) + ".bai")
    bai2 = path.with_suffix(".bai")

    if bai1.exists() or bai2.exists():
        return True

    try:
        print(f"    Indexing BAM: {path}")
        pysam.index(str(path))
        return True
    except Exception as e:
        print(f"    WARNING: could not index {path}: {e}")
        return False

def resolve_bam_chrom(bam, chrom):
    refs = set(bam.references)

    candidates = [chrom]
    if chrom.startswith("chr"): candidates.append(chrom[3:])
    else: candidates.append("chr" + chrom)

    for c in candidates:
        if c in refs: return c

    return None

def mapped_reads_from_star(tissue, sample_id):
    star = real_star_bam(tissue, sample_id)

    if not star.exists():
        print(f"    WARNING: real STAR BAM missing: {star}")
        return None

    key = str(star)
    if key in bam_mapped_cache:
        return bam_mapped_cache[key]

    try:
        if not ensure_bam_index(star): return None

        with pysam.AlignmentFile(str(star), "rb") as bam:
            n = float(bam.mapped)

        bam_mapped_cache[key] = n
        return n

    except Exception as e:
        print(f"    WARNING: could not read STAR mapped count for {star}: {e}")
        return None

def bin_array(values, n_bins):
    values = np.asarray(values, dtype=float)

    if len(values) == n_bins:
        return values

    edges = np.linspace(0, len(values), n_bins + 1).astype(int)
    out = np.zeros(n_bins)

    for i in range(n_bins):
        a, b = edges[i], edges[i + 1]
        if b > a: out[i] = np.mean(values[a:b])

    return out

def remap_bam_vector(tissue, sample_id, chrom, start, end, n_bins):
    # IMPORTANT: use the real BAM directly, not the broken /home symlink.
    bam_path = real_remap_bam(tissue, sample_id)

    if not bam_path.exists():
        print(f"    WARNING: real remap BAM missing: {bam_path}")
        return None, None

    if not ensure_bam_index(bam_path):
        return None, str(bam_path)

    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            c = resolve_bam_chrom(bam, chrom)

            if c is None:
                print(f"    WARNING: chromosome {chrom} not found in {bam_path}")
                return None, str(bam_path)

            cov = np.sum(
                np.asarray(
                    bam.count_coverage(c, int(start), int(end), quality_threshold=0),
                    dtype=float
                ),
                axis=0
            )

        denominator = mapped_reads_from_star(tissue, sample_id)

        if denominator is None or denominator <= 0:
            print(f"    WARNING: no usable STAR read-count denominator for {sample_id}")
            return None, str(bam_path)

        rpm = cov * 1_000_000.0 / denominator
        return bin_array(rpm, n_bins), str(bam_path)

    except Exception as e:
        print(f"    WARNING: BAM extraction failed for {bam_path}: {e}")
        return None, str(bam_path)

# ============================================================
# SNP SELECTION
# ============================================================

def possible_snps(hit):
    candidates = []

    for col in ["best_candidate_SNP", "best_candidate_snp"]:
        if col in hit and not missing(hit[col]): candidates.append(hit[col])

    for col in {
        "eQTL": ["top_eQTL_SNP", "top_QTL_SNP"],
        "sQTL": ["top_sQTL_SNP", "top_QTL_SNP"],
        "cQTL": ["top_cQTL_SNP", "top_QTL_SNP"]
    }[TYPE]:
        if col in hit and not missing(hit[col]): candidates.append(hit[col])

    for k, v in hit.items():
        if "topsnp" in norm(k) and "gwas" not in norm(k) and not missing(v):
            candidates.append(v)

    if "candidate_SNPs" in hit and not missing(hit["candidate_SNPs"]):
        candidates += [x for x in hit["candidate_SNPs"].split(";") if x]

    out = []
    for x in map(str.strip, candidates):
        if x and x not in out: out.append(x)

    return out

def genotype_group(value):
    try: x = float(value)
    except: return None

    if not np.isfinite(x) or x < -0.1 or x > 2.1:
        return None

    return int(np.clip(np.rint(x), 0, 2))

# ============================================================
# PROCESS
# ============================================================

hits_by_tissue = defaultdict(list)

for hit in hits:
    hits_by_tissue[hit.get(tissue_col, "").strip()].append(hit)

total_plotted = total_skipped = 0

for tissue, tissue_hits in hits_by_tissue.items():
    print(f"\n{'=' * 60}\nTISSUE: {tissue}\nPASSING HITS: {len(tissue_hits)}\n{'=' * 60}")

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

    missing_files = [p for p in [phenotype_file, trait_location_file, snp_file] if not p.exists()]

    if missing_files:
        print("SKIPPING TISSUE: missing " + ", ".join(map(str, missing_files)))
        total_skipped += len(tissue_hits)
        continue

    events = {h[event_col].strip() for h in tissue_hits if not missing(h.get(event_col))}
    all_snps = {s for h in tissue_hits for s in possible_snps(h)}

    phenotype_rows = extract_matrix_rows(phenotype_file, events)
    genotype_rows = extract_matrix_rows(snp_file, all_snps)
    trait_locations = load_location_table(trait_location_file, "trait")
    snp_locations = load_location_table(snp_location_file, "snp") if snp_location_file.exists() else {}

    for index, hit in enumerate(tissue_hits, 1):
        event = hit[event_col].strip()
        symbol = hit.get(symbol_col, "").strip() if symbol_col else ""
        pp_h4 = hit.get("max_PP_H4", "")

        print(f"\n[{index}/{len(tissue_hits)}] {tissue} | {event} | {symbol}")

        phenotype = find_matrix_row(phenotype_rows, event)

        if phenotype is None:
            print("  SKIP: phenotype missing")
            total_skipped += 1
            continue

        snp = genotype = None

        for candidate in possible_snps(hit):
            if candidate in genotype_rows:
                snp, genotype = candidate, genotype_rows[candidate]
                break

        if genotype is None:
            print("  SKIP: candidate SNP absent from genotype matrix")
            total_skipped += 1
            continue

        loc = find_location(trait_locations, event)

        if loc is None:
            print("  SKIP: coordinates unavailable")
            total_skipped += 1
            continue

        chrom = loc["chr"]
        trait_start, trait_end = min(loc["start"], loc["end"]), max(loc["start"], loc["end"])
        strand = loc.get("strand")

        flank = min(
            int(max(FLANK_MIN, math.ceil(max(1, trait_end - trait_start) * FLANK_FRAC))),
            FLANK_MAX
        )

        plot_start, plot_end = max(0, trait_start - flank), trait_end + flank
        n_bins = int(
            min(MAX_TRACK_POINTS, max(50, math.ceil((plot_end - plot_start) / BIGWIG_BIN_BP)))
        )

        x = np.linspace(plot_start, plot_end, n_bins, endpoint=False)
        x += (plot_end - plot_start) / n_bins / 2

        snp_pos = None

        if snp in snp_locations:
            snp_chr, p = snp_locations[snp]

            if norm(snp_chr) == norm(chrom):
                snp_pos = p

        # genotype + MatrixEQTL phenotype + metadata + sample directory
        records = []

        for qid in [x for x in phenotype if x in genotype]:
            g = genotype_group(genotype[qid])

            if g is None:
                continue

            try:
                pheno = float(phenotype[qid])
            except:
                continue

            if not np.isfinite(pheno):
                continue

            meta = metadata_for_qtl_id(qid, tissue)

            if meta is None:
                continue

            sample_id = str(meta.get(sample_col, "")).strip()
            subject_id = str(meta.get(subject_col, "")).strip()

            if not sample_id:
                continue

            sample_dir = resolve_sample_dir(tissue, sample_id)

            if sample_dir is None:
                continue

            records.append({
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

        print(f"  Candidate subjects before track QC: {len(records)}")

        for r in records:
            if TYPE in {"eQTL", "sQTL"}:
                r["track"], r["track_source"] = bigwig_vector(
                    r["sample_dir"], chrom, plot_start, plot_end, strand, n_bins
                )
            else:
                r["track"], r["track_source"] = remap_bam_vector(
                    tissue, r["sample_id"], chrom, plot_start, plot_end, n_bins
                )

        # Same complete subjects in phenotype AND sequencing panels.
        complete = [
            r for r in records
            if r["track"] is not None
            and len(r["track"]) == n_bins
            and np.all(np.isfinite(r["track"]))
        ]

        counts = {
            g: sum(r["genotype"] == g for r in complete)
            for g in [0, 1, 2]
        }

        print(f"  Complete subjects: {len(complete)}")
        print(f"  Ref/Ref={counts[0]} | Het={counts[1]} | Hom Alt={counts[2]}")

        if not complete:
            print("  SKIP: no complete subjects")
            total_skipped += 1
            continue

        tracks = {
            g: [r["track"] for r in complete if r["genotype"] == g]
            for g in [0, 1, 2]
        }

        label = symbol if symbol else event
        basename = safe_name(f"{tissue}__{TYPE}__{label}__{snp}")

        pdf = OUTDIR / f"{basename}.pdf"
        png = OUTDIR / f"{basename}.png"
        subject_tsv = OUTDIR / f"{basename}.subjects.tsv"

        if not OVERWRITE and pdf.exists() and png.exists():
            print("  EXISTS; skipping")
            continue

        with open(subject_tsv, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")

            writer.writerow([
                "tissue", "qtl_type", "event", "gene_symbol", "snp", "qtl_id",
                "subject_id", "sample_id", "genotype_group", "genotype_label",
                "genotype_raw", "phenotype", "track_source"
            ])

            for r in complete:
                writer.writerow([
                    tissue, TYPE, event, symbol, snp, r["qtl_id"], r["subject_id"],
                    r["sample_id"], r["genotype"], GENOTYPE_LABELS[r["genotype"]],
                    r["genotype_raw"], r["phenotype"], r["track_source"]
                ])

        ymax = max(float(np.nanmax(r["track"])) for r in complete)
        ymax = 1.0 if not np.isfinite(ymax) or ymax <= 0 else ymax * 1.05

        # ====================================================
        # PLOT
        # ====================================================

        fig = plt.figure(figsize=(16, 8))
        gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.8], hspace=0.35, wspace=0.08)

        # ---------- TOP: exact MatrixEQTL phenotype ----------

        ax0 = fig.add_subplot(gs[0, :])
        rng = np.random.default_rng(12345)

        for g in [0, 1, 2]:
            vals = np.asarray([
                r["phenotype"] for r in complete
                if r["genotype"] == g
            ])

            if not len(vals):
                continue

            mean_val = np.mean(vals)
            jitter = rng.normal(0, 0.055, size=len(vals))

            ax0.scatter(
                np.full(len(vals), g) + jitter, vals,
                s=24, alpha=0.65, color=TOP_COLORS[g], edgecolor="none"
            )

            ax0.plot(
                [g - 0.18, g + 0.18], [mean_val, mean_val],
                linewidth=3, color=TOP_COLORS[g]
            )

            ax0.annotate(
                f"Mean = {mean_val:.4f}",
                xy=(g - 0.18, mean_val), xytext=(0, 7),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=8.5, color="black"
            )

        ax0.set_xticks(
            [0, 1, 2],
            [
                f"Ref/Ref\n(N={counts[0]})",
                f"Het\n(N={counts[1]})",
                f"Hom Alt\n(N={counts[2]})"
            ]
        )

        ax0.set_ylabel(Y_LABELS[TYPE])

        title_event = f"{symbol} ({event})" if symbol and symbol != event else event
        title = f"{tissue_display(tissue)} | {TYPE} | {title_event} | {snp}"

        if pp_h4 and not missing(pp_h4):
            try:
                title += f" | PP.H4={float(pp_h4):.4f}"
            except:
                title += f" | PP.H4={pp_h4}"

        ax0.set_title(title, fontsize=13)
        ax0.spines["top"].set_visible(False)
        ax0.spines["right"].set_visible(False)

        # ---------- BOTTOM: genomic sequencing signal ----------

        for g in [0, 1, 2]:
            ax = fig.add_subplot(gs[1, g])
            group_tracks = tracks[g]

            for arr in group_tracks:
                ax.plot(x, arr, linewidth=0.65, alpha=0.16, color=TRACK_COLORS[g])

            if group_tracks:
                ax.plot(
                    x, np.nanmean(np.vstack(group_tracks), axis=0),
                    linewidth=2.5, color=TRACK_COLORS[g]
                )
            else:
                ax.text(
                    0.5, 0.5, "No subjects",
                    transform=ax.transAxes, ha="center", va="center"
                )

            ax.axvspan(trait_start, trait_end, alpha=0.12, color="#999999")

            if snp_pos is not None and plot_start <= snp_pos <= plot_end:
                ax.axvline(
                    snp_pos, linestyle="--",
                    linewidth=1.3, alpha=0.8, color="black"
                )

            ax.set_xlim(plot_start, plot_end)
            ax.set_ylim(0, ymax)
            ax.set_title(f"{GENOTYPE_LABELS[g]} (N={counts[g]})")
            ax.set_xlabel(f"{chrom} genomic position")

            if g == 0:
                ax.set_ylabel(
                    "Remapped read coverage\n(RPM-normalized)"
                    if TYPE == "cQTL"
                    else "Normalized RNA-seq read density"
                )
            else:
                ax.tick_params(labelleft=False)

            ax.xaxis.set_major_formatter(FuncFormatter(genomic_formatter))
            ax.locator_params(axis="x", nbins=5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        fig.text(
            0.5, 0.015,
            f"Trait: {chrom}:{trait_start:,}-{trait_end:,} | "
            f"flank: {flank:,} bp each side | shaded region = trait locus",
            ha="center", fontsize=9
        )

        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, dpi=DPI, bbox_inches="tight")
        plt.close(fig)

        print(f"  PDF: {pdf}")
        print(f"  PNG: {png}")
        print(f"  TSV: {subject_tsv}")

        total_plotted += 1

print(f"\n{'=' * 60}")
print("DONE")
print(f"QTL type        : {TYPE}")
print(f"Plots generated : {total_plotted}")
print(f"Skipped hits    : {total_skipped}")
print(f"Output          : {OUTDIR}")
print("=" * 60)
PY
