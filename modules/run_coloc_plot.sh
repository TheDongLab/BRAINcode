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
MR_DIR="${ROOT}/MR"
METADATA="${ROOT}/targetALS_rnaseq_metadata.csv"
GTF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf"

FLANK_FRAC=0.20
FLANK_MIN=2000
FLANK_MAX=20000
BIGWIG_BIN_BP=25
MAX_TRACK_POINTS=2500
CIRC_COORD_TOL=3                      # exact first; +/-3 bp fallback
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

[[ -f "$COLOC" ]] || { echo "ERROR: missing $COLOC"; exit 1; }
[[ -f "$METADATA" ]] || { echo "ERROR: missing $METADATA"; exit 1; }
[[ -f "$GTF" ]] || { echo "ERROR: missing $GTF"; exit 1; }

module load deepTools 2>/dev/null || true

export TYPE ROOT METADATA COLOC OUTDIR GTF
export FLANK_FRAC FLANK_MIN FLANK_MAX BIGWIG_BIN_BP MAX_TRACK_POINTS CIRC_COORD_TOL
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
    from matplotlib.patches import Rectangle
    import pyBigWig
except ImportError as e:
    sys.exit(f"ERROR: missing Python package: {e}\nNeed numpy, matplotlib and pyBigWig.")

warnings.filterwarnings("ignore")

TYPE, ROOT = os.environ["TYPE"], Path(os.environ["ROOT"])
METADATA, COLOC, OUTDIR, GTF = map(Path, [os.environ["METADATA"], os.environ["COLOC"], os.environ["OUTDIR"], os.environ["GTF"]])
FLANK_FRAC, FLANK_MIN, FLANK_MAX = float(os.environ["FLANK_FRAC"]), int(os.environ["FLANK_MIN"]), int(os.environ["FLANK_MAX"])
BIGWIG_BIN_BP, MAX_TRACK_POINTS = int(os.environ["BIGWIG_BIN_BP"]), int(os.environ["MAX_TRACK_POINTS"])
CIRC_COORD_TOL = int(os.environ["CIRC_COORD_TOL"])
MAX_HITS, TISSUE_FILTER = int(os.environ["MAX_HITS"]), os.environ["TISSUE_FILTER"].strip()
OVERWRITE, DPI = int(os.environ["OVERWRITE"]), int(os.environ["DPI"])

GENOTYPE_LABELS = {0: "Ref/Ref", 1: "Het", 2: "Hom Alt"}
TOP_COLORS = {0: "#90EE90", 1: "#ADD8E6", 2: "#FFB6C1"}
TRACK_COLORS = {0: "#228B22", 1: "#4682B4", 2: "#C75B7A"}
Y_LABELS = {"eQTL": "Z-score normalized expression", "sQTL": "Percent spliced in", "cQTL": "Percent circularized"}

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

def format_max_sig(value, sig=5):
    if value == 0: return "0"
    decimals = max(0, sig - int(math.floor(math.log10(abs(value)))) - 1)
    return f"{value:.{decimals}f}".rstrip("0").rstrip(".")

def genomic_formatter(x, pos=None):
    if abs(x) >= 1_000_000: return f"{format_max_sig(x / 1_000_000)} Mb"
    if abs(x) >= 1_000: return f"{format_max_sig(x / 1_000)} kb"
    return format_max_sig(x)

def find_header(headers, aliases, required=True):
    lookup = {norm(h): h for h in headers}
    for a in aliases:
        if norm(a) in lookup: return lookup[norm(a)]
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

if TISSUE_FILTER: hits = [h for h in hits if tissue_equal(h.get("tissue", ""), TISSUE_FILTER)]
if MAX_HITS > 0: hits = hits[:MAX_HITS]
if not hits: sys.exit("No passing coloc hits to plot.")

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
    qtl_id, rows = str(qtl_id), metadata

    if meta_tissue_col:
        tissue_rows = [r for r in metadata if tissue_equal(r.get(meta_tissue_col, ""), tissue)]
        if tissue_rows: rows = tissue_rows

    exact = [r for r in rows if str(r.get(sample_col, "")).strip() == qtl_id]
    if exact: return max(exact, key=get_rin)

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
# LOCATION TABLES
# ============================================================

def load_location_table(path, kind):
    with open(path) as fh:
        first = fh.readline()
        delim = detect_delimiter(first)
        fields = split_line(first, delim)
        is_data = len(fields) >= 3 and re.match(r"^(chr)?[0-9XYM]+$", fields[1], re.I) and re.match(r"^\d+$", fields[2])

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
            for r in rows if not missing(r.get(idc)) and not missing(r.get(cc)) and not missing(r.get(pc))
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
                "chr": str(r[cc]),
                "start": int(float(r[sc])),
                "end": int(float(r[ec])),
                "strand": str(r[strandc]) if strandc and not missing(r.get(strandc)) else None
            }
        except: pass

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
# SAMPLE DIRECTORY RESOLUTION
# ============================================================

processed_cache = {}

def resolve_sample_dir(tissue, sample_id):
    processed = ROOT / tissue / "RNAseq" / "Processed"

    for candidate in [sample_id, sample_id.replace("-", "_"), sample_id.replace("_", "-")]:
        p = processed / candidate
        if p.is_dir(): return p

    key = str(processed)
    if key not in processed_cache:
        processed_cache[key] = {norm(p.name): p for p in processed.iterdir() if p.is_dir()} if processed.is_dir() else {}

    return processed_cache[key].get(norm(sample_id))

# ============================================================
# NORMALIZED BIGWIG
# ============================================================

def bw_chrom(bw, chrom):
    chroms = bw.chroms()
    candidates = [chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom]
    return next((c for c in candidates if c in chroms), None)

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
# cQTL QUANTIFICATION — ANNOTATION ONLY
# ============================================================

def circ_quant_file(sample_dir):
    return sample_dir / "circularRNA_known_circ_percentage.txt"

def read_circ_quant(path, chrom, start, end, strand, tolerance=CIRC_COORD_TOL):
    if not path.exists(): return None

    exact, nearby = None, []

    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")

            for row in reader:
                try:
                    rchr = row["chrom"]
                    rstart, rend = int(float(row["start"])), int(float(row["end"]))
                    rstrand = str(row["strand"]).strip()

                    if norm(rchr) != norm(chrom): continue
                    if strand in {"+", "-"} and rstrand != str(strand).strip(): continue

                    result = {
                        "circ_reads": float(row["circ_reads"]),
                        "linear_reads": float(row["linear_reads"]),
                        "circ_percent": float(row["circ_percent"]),
                        "matched_start": rstart, "matched_end": rend, "match_type": None
                    }

                    if rstart == int(start) and rend == int(end):
                        result["match_type"] = "exact"
                        exact = result
                        break

                    if abs(rstart - int(start)) <= tolerance and abs(rend - int(end)) <= tolerance:
                        result["match_type"] = "tolerance"
                        nearby.append(result)

                except (ValueError, KeyError, TypeError):
                    continue

    except Exception as e:
        print(f"    WARNING: could not read circ quant file {path}: {e}")
        return None

    if exact is not None: return exact
    if not nearby: return None

    nearby.sort(key=lambda r: (
        abs(r["matched_start"] - int(start)) + abs(r["matched_end"] - int(end)),
        abs(r["matched_start"] - int(start)),
        abs(r["matched_end"] - int(end))
    ))

    return nearby[0]

# ============================================================
# GENCODE v49 EXON TRACK
# ============================================================

def parse_gtf_attrs(s):
    attrs = {}
    for item in s.strip().split(";"):
        item = item.strip()
        if not item: continue
        p = item.split(" ", 1)
        if len(p) == 2: attrs[p[0]] = p[1].strip().strip('"')
    return attrs

def load_reference_exons(gtf, chrom, start, end):
    exons = []

    with open(gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")

            if len(f) < 9 or f[2] != "exon" or norm(f[0]) != norm(chrom): continue

            try:
                exon_start, exon_end = int(f[3]) - 1, int(f[4])
            except ValueError:
                continue

            if exon_end < start or exon_start > end: continue

            a = parse_gtf_attrs(f[8])

            exons.append({
                "start": exon_start,
                "end": exon_end,
                "strand": f[6],
                "gene": a.get("gene_name", a.get("gene_id", "Unknown")),
                "gene_id": a.get("gene_id", "")
            })

    return exons

def merge_intervals(intervals):
    if not intervals: return []

    intervals = sorted((int(s), int(e)) for s, e in intervals)
    merged = [list(intervals[0])]

    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])

    return merged

def draw_reference_track(ax, exons, plot_start, plot_end, show_ylabel=False):
    genes, strands = defaultdict(list), {}

    for exon in exons:
        gene = exon["gene"] or exon["gene_id"] or "Unknown"
        genes[gene].append((exon["start"], exon["end"]))
        strands[gene] = exon["strand"]

    if not genes:
        ax.text(0.5, 0.45, "No GENCODE v49 exons", transform=ax.transAxes, ha="center", va="center", fontsize=7)

    else:
        gene_items = sorted(genes.items(), key=lambda kv: min(s for s, e in kv[1]))

        for y, (gene, intervals) in enumerate(gene_items):
            merged = merge_intervals(intervals)
            span_start = max(plot_start, min(s for s, e in merged))
            span_end = min(plot_end, max(e for s, e in merged))

            ax.plot([span_start, span_end], [y, y], color="black", linewidth=0.8)

            for s, e in merged:
                s, e = max(s, plot_start), min(e, plot_end)
                if e <= s: continue
                ax.add_patch(Rectangle((s, y - 0.15), e - s, 0.30, facecolor="black", edgecolor="black", linewidth=0.4))

            strand = strands.get(gene, "")
            label = f"{gene} ({strand})" if strand in {"+", "-"} else gene
            ax.text(plot_start, y + 0.22, label, fontsize=6.5, ha="left", va="bottom")

        ax.set_ylim(-0.45, max(0.55, len(gene_items) - 0.25))

    ax.set_xlim(plot_start, plot_end)
    ax.set_yticks([])
    ax.set_ylabel("GENCODE v49" if show_ylabel else "", fontsize=8)
    ax.set_xlabel("")
    ax.xaxis.set_major_formatter(FuncFormatter(genomic_formatter))
    ax.locator_params(axis="x", nbins=4)

    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)

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
        if "topsnp" in norm(k) and "gwas" not in norm(k) and not missing(v): candidates.append(v)

    if "candidate_SNPs" in hit and not missing(hit["candidate_SNPs"]):
        candidates += [x for x in hit["candidate_SNPs"].split(";") if x]

    out = []
    for x in map(str.strip, candidates):
        if x and x not in out: out.append(x)

    return out

def genotype_group(value):
    try: x = float(value)
    except: return None

    if not np.isfinite(x) or x < -0.1 or x > 2.1: return None
    return int(np.clip(np.rint(x), 0, 2))

# ============================================================
# PROCESS
# ============================================================

hits_by_tissue = defaultdict(list)
for hit in hits: hits_by_tissue[hit.get(tissue_col, "").strip()].append(hit)

total_plotted = total_skipped = 0

for tissue, tissue_hits in hits_by_tissue.items():
    print(f"\n{'=' * 60}\nTISSUE: {tissue}\nPASSING HITS: {len(tissue_hits)}\n{'=' * 60}")

    qtl_dir = ROOT / tissue / TYPE

    if TYPE == "eQTL":
        phenotype_file, trait_location_file = qtl_dir / f"expression_{tissue}.txt", qtl_dir / "gene_location.txt"
    elif TYPE == "sQTL":
        phenotype_file, trait_location_file = qtl_dir / f"splicing_{tissue}.txt", qtl_dir / "splicing_location.txt"
    else:
        phenotype_file, trait_location_file = qtl_dir / f"circ_{tissue}.txt", qtl_dir / "circ_location.txt"

    snp_file, snp_location_file = qtl_dir / f"snp_{tissue}.txt", qtl_dir / "snp_location.txt"
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
        event, pp_h4 = hit[event_col].strip(), hit.get("max_PP_H4", "")
        symbol = hit.get(symbol_col, "").strip() if symbol_col else ""

        print(f"\n[{index}/{len(tissue_hits)}] {tissue} | {event} | {symbol}")

        phenotype = find_matrix_row(phenotype_rows, event)
        if phenotype is None:
            print("  SKIP: phenotype missing"); total_skipped += 1; continue

        snp = genotype = None
        for candidate in possible_snps(hit):
            if candidate in genotype_rows:
                snp, genotype = candidate, genotype_rows[candidate]
                break

        if genotype is None:
            print("  SKIP: candidate SNP absent from genotype matrix"); total_skipped += 1; continue

        loc = find_location(trait_locations, event)
        if loc is None:
            print("  SKIP: coordinates unavailable"); total_skipped += 1; continue

        chrom = loc["chr"]
        trait_start, trait_end = min(loc["start"], loc["end"]), max(loc["start"], loc["end"])
        strand = loc.get("strand")

        flank = min(int(max(FLANK_MIN, math.ceil(max(1, trait_end - trait_start) * FLANK_FRAC))), FLANK_MAX)
        plot_start, plot_end = max(0, trait_start - flank), trait_end + flank
        n_bins = int(min(MAX_TRACK_POINTS, max(50, math.ceil((plot_end - plot_start) / BIGWIG_BIN_BP))))

        x = np.linspace(plot_start, plot_end, n_bins, endpoint=False)
        x += (plot_end - plot_start) / n_bins / 2

        snp_pos = None
        if snp in snp_locations:
            snp_chr, p = snp_locations[snp]
            if norm(snp_chr) == norm(chrom): snp_pos = p

        reference_exons = load_reference_exons(GTF, chrom, plot_start, plot_end)
        print(f"  GENCODE v49 overlapping exons: {len(reference_exons)}")

        records = []

        for qid in [q for q in phenotype if q in genotype]:
            g = genotype_group(genotype[qid])
            if g is None: continue

            try: pheno = float(phenotype[qid])
            except: continue
            if not np.isfinite(pheno): continue

            meta = metadata_for_qtl_id(qid, tissue)
            if meta is None: continue

            sample_id = str(meta.get(sample_col, "")).strip()
            subject_id = str(meta.get(subject_col, "")).strip()

            if not sample_id: continue

            sample_dir = resolve_sample_dir(tissue, sample_id)
            if sample_dir is None: continue

            records.append({
                "qtl_id": qid, "subject_id": subject_id, "sample_id": sample_id,
                "sample_fs_id": sample_dir.name, "genotype": g, "genotype_raw": genotype[qid],
                "phenotype": pheno, "sample_dir": sample_dir, "track": None, "track_source": None,
                "circ_reads": None, "linear_reads": None, "circ_percent_file": None,
                "circ_quant_source": None, "circ_match_type": None,
                "circ_matched_start": None, "circ_matched_end": None
            })

        print(f"  Candidate subjects before track QC: {len(records)}")

        for r in records:
            r["track"], r["track_source"] = bigwig_vector(
                r["sample_dir"], chrom, plot_start, plot_end, strand, n_bins
            )

            if TYPE == "cQTL":
                qfile = circ_quant_file(r["sample_dir"])
                r["circ_quant_source"] = str(qfile)
                q = read_circ_quant(qfile, chrom, trait_start, trait_end, strand)

                if q is not None:
                    r["circ_reads"], r["linear_reads"] = q["circ_reads"], q["linear_reads"]
                    r["circ_percent_file"], r["circ_match_type"] = q["circ_percent"], q["match_type"]
                    r["circ_matched_start"], r["circ_matched_end"] = q["matched_start"], q["matched_end"]

        complete = [
            r for r in records
            if r["track"] is not None and len(r["track"]) == n_bins and np.all(np.isfinite(r["track"]))
        ]

        counts = {g: sum(r["genotype"] == g for r in complete) for g in [0, 1, 2]}

        print(f"  Complete subjects used in plots: {len(complete)}")
        print(f"  Ref/Ref={counts[0]} | Het={counts[1]} | Hom Alt={counts[2]}")

        if not complete:
            print("  SKIP: no subjects with usable tracks"); total_skipped += 1; continue

        if TYPE == "cQTL":
            n_quant = sum(r["circ_reads"] is not None for r in complete)
            n_exact = sum(r["circ_match_type"] == "exact" for r in complete)
            n_tol = sum(r["circ_match_type"] == "tolerance" for r in complete)
            print(f"  Circ quant rows found: {n_quant}/{len(complete)}")
            print(f"  Circ coordinate matches: exact={n_exact}, tolerance-rescued={n_tol}")

        tracks = {g: [r["track"] for r in complete if r["genotype"] == g] for g in [0, 1, 2]}

        label = symbol if symbol else event
        basename = safe_name(f"{tissue}__{TYPE}__{label}__{snp}")
        pdf, png, svg = OUTDIR / f"{basename}.pdf", OUTDIR / f"{basename}.png", OUTDIR / f"{basename}.svg"
        subject_tsv = OUTDIR / f"{basename}.subjects.tsv"

        if not OVERWRITE and pdf.exists() and png.exists() and svg.exists():
            print("  EXISTS; skipping"); continue

        with open(subject_tsv, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([
                "tissue", "qtl_type", "event", "gene_symbol", "snp", "qtl_id", "subject_id",
                "metadata_sample_id", "filesystem_sample_id", "genotype_group", "genotype_label",
                "genotype_raw", "phenotype", "track_source", "circ_reads", "linear_reads",
                "circ_percent_file", "circ_match_type", "circ_matched_start", "circ_matched_end",
                "circ_quant_source"
            ])

            for r in complete:
                writer.writerow([
                    tissue, TYPE, event, symbol, snp, r["qtl_id"], r["subject_id"], r["sample_id"],
                    r["sample_fs_id"], r["genotype"], GENOTYPE_LABELS[r["genotype"]],
                    r["genotype_raw"], r["phenotype"], r["track_source"],
                    "" if r["circ_reads"] is None else r["circ_reads"],
                    "" if r["linear_reads"] is None else r["linear_reads"],
                    "" if r["circ_percent_file"] is None else r["circ_percent_file"],
                    "" if r["circ_match_type"] is None else r["circ_match_type"],
                    "" if r["circ_matched_start"] is None else r["circ_matched_start"],
                    "" if r["circ_matched_end"] is None else r["circ_matched_end"],
                    "" if r["circ_quant_source"] is None else r["circ_quant_source"]
                ])

        ymax = max(float(np.nanmax(r["track"])) for r in complete)
        ymax = 1.0 if not np.isfinite(ymax) or ymax <= 0 else ymax * 1.05

        # ====================================================
        # PLOT
        # ============================================================

        fig = plt.figure(figsize=(16, 8.5))
        gs = fig.add_gridspec(
            3, 3,
            height_ratios=[1.0, 1.75, 0.42],
            hspace=0.16,
            wspace=0.08
        )

        # ---------- TOP: MatrixEQTL phenotype ----------

        ax0, rng = fig.add_subplot(gs[0, :]), np.random.default_rng(12345)

        for g in [0, 1, 2]:
            vals = np.asarray([r["phenotype"] for r in complete if r["genotype"] == g])
            if not len(vals): continue

            mean_val, jitter = np.mean(vals), rng.normal(0, 0.055, size=len(vals))

            ax0.scatter(
                np.full(len(vals), g) + jitter, vals,
                s=24, alpha=0.65, color=TOP_COLORS[g], edgecolor="none"
            )
            ax0.plot([g - 0.18, g + 0.18], [mean_val, mean_val], linewidth=3, color=TOP_COLORS[g])
            ax0.annotate(
                f"Mean = {mean_val:.4f}", xy=(g - 0.18, mean_val), xytext=(0, 7),
                textcoords="offset points", ha="left", va="bottom", fontsize=8.5, color="black"
            )

        ax0.set_xticks([])
        ax0.set_ylabel(Y_LABELS[TYPE])

        title_event = f"{symbol} ({event})" if symbol and symbol != event else event
        title = f"{tissue_display(tissue)} | {TYPE} | {title_event} | {snp}"

        if pp_h4 and not missing(pp_h4):
            try: title += f" | PP.H4={float(pp_h4):.4f}"
            except: title += f" | PP.H4={pp_h4}"

        ax0.set_title(title, fontsize=13)
        ax0.spines["top"].set_visible(False)
        ax0.spines["right"].set_visible(False)

        # ---------- MIDDLE + BOTTOM: reads + triplicate GENCODE ----------

        for g in [0, 1, 2]:
            ax = fig.add_subplot(gs[1, g])
            ax_ref = fig.add_subplot(gs[2, g], sharex=ax)
            group_tracks = tracks[g]

            for arr in group_tracks:
                ax.plot(x, arr, linewidth=0.65, alpha=0.16, color=TRACK_COLORS[g])

            if group_tracks:
                ax.plot(x, np.nanmean(np.vstack(group_tracks), axis=0), linewidth=2.5, color=TRACK_COLORS[g])
            else:
                ax.text(0.5, 0.5, "No subjects", transform=ax.transAxes, ha="center", va="center")

            ax.axvspan(trait_start, trait_end, alpha=0.12, color="#999999")

            if snp_pos is not None and plot_start <= snp_pos <= plot_end:
                ax.axvline(snp_pos, linestyle="--", linewidth=1.3, alpha=0.8, color="black")

            if TYPE == "cQTL":
                gr = [r for r in complete if r["genotype"] == g]
                qgr = [r for r in gr if r["circ_reads"] is not None]

                if gr and qgr:
                    cv = np.asarray([r["circ_reads"] for r in qgr])
                    lv = np.asarray([r["linear_reads"] for r in qgr])

                    annotation = (
                        f"Mean circ reads = {np.mean(cv):.2f}\n"
                        f"Mean linear reads = {np.mean(lv):.2f}\n"
                        f"Circ+ samples = {int(np.sum(cv > 0))}/{len(qgr)}\n"
                        f"Circ quant rows = {len(qgr)}/{len(gr)}"
                    )

                    ax.text(
                        0.02, 0.96, annotation,
                        transform=ax.transAxes, ha="left", va="top",
                        fontsize=9, color="black"
                    )

                elif gr:
                    ax.text(
                        0.02, 0.96, f"Circ quant rows = 0/{len(gr)}",
                        transform=ax.transAxes, ha="left", va="top",
                        fontsize=9, color="black"
                    )

            ax.set_xlim(plot_start, plot_end)
            ax.set_ylim(0, ymax)
            ax.set_title(f"{GENOTYPE_LABELS[g]} (N={counts[g]})")
            ax.tick_params(axis="x", labelbottom=False)

            if g == 0: ax.set_ylabel("Normalized RNA-seq read density")
            else: ax.tick_params(labelleft=False)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            # GENCODE track directly under THIS genotype panel.
            draw_reference_track(ax_ref, reference_exons, plot_start, plot_end, show_ylabel=(g == 0))
            ax_ref.set_xlabel(f"{chrom} genomic position")

        fig.text(
            0.5, 0.008,
            f"Trait: {chrom}:{trait_start:,}-{trait_end:,} | "
            f"flank: {flank:,} bp each side | shaded region = trait locus",
            ha="center", fontsize=9
        )

        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, dpi=DPI, bbox_inches="tight")
        fig.savefig(svg, bbox_inches="tight")
        plt.close(fig)

        print(f"  PDF: {pdf}")
        print(f"  PNG: {png}")
        print(f"  SVG: {svg}")
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
