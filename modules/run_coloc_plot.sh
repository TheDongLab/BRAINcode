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
RAW_FILE="${ROOT}/QTL/plink/joint_all_chrs_matrixEQTL.raw"

FLANK_FRAC=0.20
FLANK_MIN=2000
FLANK_MAX=20000
BIGWIG_BIN_BP=25
MAX_TRACK_POINTS=2500
CIRC_COORD_TOL=3
MIN_MEAN_JUNCTION_READS=0.5           # draw linear junctions above this mean
MAX_JUNCTIONS=30                      # strongest junctions per genotype
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

for f in "$COLOC" "$METADATA" "$GTF" "$RAW_FILE"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done

module load deepTools 2>/dev/null || { echo "ERROR: could not load deepTools"; exit 1; }
module load poppler/25.07.0-GCC-13.3.0 2>/dev/null || { echo "ERROR: could not load Poppler"; exit 1; }

PYTHON="$(command -v python)"
command -v pdfseparate >/dev/null || { echo "ERROR: pdfseparate unavailable"; exit 1; }
command -v pdftotext >/dev/null || { echo "ERROR: pdftotext unavailable"; exit 1; }
"$PYTHON" -c 'import numpy,matplotlib,pyBigWig' || { echo "ERROR: Python plotting stack unavailable"; exit 1; }

export TYPE ROOT METADATA COLOC OUTDIR GTF RAW_FILE
export FLANK_FRAC FLANK_MIN FLANK_MAX BIGWIG_BIN_BP MAX_TRACK_POINTS CIRC_COORD_TOL
export MIN_MEAN_JUNCTION_READS MAX_JUNCTIONS MAX_HITS TISSUE_FILTER OVERWRITE DPI

"$PYTHON" - <<'PY'
import os, re, csv, sys, math, warnings, subprocess
from pathlib import Path
from collections import defaultdict

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    from matplotlib.patches import Rectangle, PathPatch
    from matplotlib.path import Path as MplPath
    import pyBigWig
except ImportError as e:
    sys.exit(f"ERROR: missing Python package: {e}\nNeed numpy, matplotlib and pyBigWig.")

warnings.filterwarnings("ignore")

TYPE, ROOT = os.environ["TYPE"], Path(os.environ["ROOT"])
METADATA, COLOC, OUTDIR, GTF, RAW_FILE = map(
    Path, [os.environ["METADATA"], os.environ["COLOC"], os.environ["OUTDIR"], os.environ["GTF"], os.environ["RAW_FILE"]]
)
FLANK_FRAC, FLANK_MIN, FLANK_MAX = float(os.environ["FLANK_FRAC"]), int(os.environ["FLANK_MIN"]), int(os.environ["FLANK_MAX"])
BIGWIG_BIN_BP, MAX_TRACK_POINTS = int(os.environ["BIGWIG_BIN_BP"]), int(os.environ["MAX_TRACK_POINTS"])
CIRC_COORD_TOL, MIN_MEAN_JUNCTION_READS = int(os.environ["CIRC_COORD_TOL"]), float(os.environ["MIN_MEAN_JUNCTION_READS"])
MAX_JUNCTIONS, MAX_HITS = int(os.environ["MAX_JUNCTIONS"]), int(os.environ["MAX_HITS"])
TISSUE_FILTER, OVERWRITE, DPI = os.environ["TISSUE_FILTER"].strip(), int(os.environ["OVERWRITE"]), int(os.environ["DPI"])

GENOTYPE_LABELS = {0: "Ref/Ref", 1: "Het", 2: "Hom Alt"}
TRACK_COLORS = {0: "#228B22", 1: "#4682B4", 2: "#C75B7A"}

# ============================================================
# HELPERS
# ============================================================

def norm(s): return re.sub(r"[^a-z0-9]", "", str(s).lower())
def chr_key(s): return re.sub(r"^chr", "", str(s), flags=re.I).upper()
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
# PLINK COUNTED ALLELES
# ============================================================

raw_alleles = {}
with open(RAW_FILE) as fh:
    for col in fh.readline().split():
        m = re.match(r"^(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)_([ACGT]+)$", col, re.I)
        if m:
            chrom, pos, ref, alt, counted = m.groups()
            raw_alleles[(chr_key(chrom), int(pos))] = {
                "ref": ref.upper(), "alt": alt.upper(), "counted": counted.upper()
            }

print(f"PLINK RAW allele mappings: {len(raw_alleles):,}")

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
                "chr": str(r[cc]), "start": int(float(r[sc])), "end": int(float(r[ec])),
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
# EXISTING COLORED BOXPLOT EXTRACTION
# ============================================================

def find_boxplot_page(pairs_file, event, snp):
    if not pairs_file.exists(): return None
    fallback, page = [], 0

    with open(pairs_file) as fh:
        first = True
        for line in fh:
            if not line.strip(): continue
            fields = line.rstrip("\r\n").split("\t") if "\t" in line else line.split()
            if len(fields) < 2: continue

            if first:
                first = False
                if norm(fields[-1]) in {"snpid", "snp"}:
                    continue

            page += 1
            row_event, row_snp = fields[0].strip(), fields[-1].strip()

            if row_event == event and row_snp == snp:
                return page

            if row_snp == snp and strip_gene_version(row_event) == strip_gene_version(event):
                fallback.append(page)

    return fallback[0] if len(fallback) == 1 else None

def extract_boxplot(source_pdf, pairs_file, event, snp, output_pdf):
    if not source_pdf.exists() or not pairs_file.exists():
        print(f"  BOXPLOT: source missing ({source_pdf} or {pairs_file})")
        return False

    page = find_boxplot_page(pairs_file, event, snp)
    if page is None:
        print(f"  BOXPLOT: no unique event+SNP match for {event} + {snp}")
        return False

    if output_pdf.exists() and not OVERWRITE:
        print(f"  BOXPLOT: {output_pdf} [exists]")
        return True

    pattern = output_pdf.with_name(output_pdf.stem + "_tmp_%d.pdf")

    try:
        subprocess.run(
            ["pdfseparate", "-f", str(page), "-l", str(page), str(source_pdf), str(pattern)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True
        )
        extracted = Path(str(pattern).replace("%d", str(page)))
        if not extracted.exists():
            print(f"  BOXPLOT: extraction failed for page {page}")
            return False

        extracted.replace(output_pdf)

        verify = subprocess.run(
            ["pdftotext", str(output_pdf), "-"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True
        ).stdout

        event_ok = event in verify or strip_gene_version(event) in verify
        snp_ok = snp in verify

        if not event_ok or not snp_ok:
            output_pdf.unlink(missing_ok=True)
            print(f"  BOXPLOT: page {page} failed event+SNP verification")
            return False

        print(f"  BOXPLOT: page {page} -> {output_pdf}")
        return True

    except Exception as e:
        print(f"  BOXPLOT: extraction error: {e}")
        return False

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
# LEAFCUTTER JUNCTION READS
# ============================================================

leafcutter_cache = {}

def leafcutter_file(sample_dir):
    return sample_dir / "leafcutter" / "psi" / f"{sample_dir.name}.leafcutter.PSI.tsv"

def load_leafcutter(path):
    key = str(path)
    if key in leafcutter_cache: return leafcutter_cache[key]
    if not path.exists():
        leafcutter_cache[key] = None
        return None

    junctions = {}
    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    k = (
                        chr_key(row["chrom"]), str(row["strand"]).strip(),
                        int(float(row["start"])), int(float(row["end"]))
                    )
                    junctions[k] = junctions.get(k, 0.0) + float(row["reads"])
                except (KeyError, ValueError, TypeError):
                    continue
    except Exception:
        junctions = None

    leafcutter_cache[key] = junctions
    return junctions

def mean_junction_reads(group, chrom, start, end, strand):
    sums, n_files = defaultdict(float), 0

    for r in group:
        junctions = load_leafcutter(leafcutter_file(r["sample_dir"]))
        if junctions is None: continue
        n_files += 1

        for (jchr, jstrand, js, je), reads in junctions.items():
            if jchr != chr_key(chrom): continue
            if strand in {"+", "-"} and jstrand != strand: continue
            if js < start or je > end: continue
            sums[(js, je)] += reads

    if n_files == 0: return {}, 0

    # Junction absent from an existing PSI file contributes zero:
    # dividing by all available genotype PSI files implements this.
    means = {k: v / n_files for k, v in sums.items()}
    means = {k: v for k, v in means.items() if v >= MIN_MEAN_JUNCTION_READS}

    if MAX_JUNCTIONS > 0 and len(means) > MAX_JUNCTIONS:
        keep = sorted(means, key=means.get, reverse=True)[:MAX_JUNCTIONS]
        means = {k: means[k] for k in keep}

    return means, n_files

# ============================================================
# cQTL QUANTIFICATION
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

                    if chr_key(rchr) != chr_key(chrom): continue
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
# GENCODE v49 TRACK
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
            if len(f) < 9 or f[2] != "exon" or chr_key(f[0]) != chr_key(chrom): continue

            try:
                exon_start, exon_end = int(f[3]) - 1, int(f[4])
            except ValueError:
                continue

            if exon_end < start or exon_start > end: continue
            a = parse_gtf_attrs(f[8])

            exons.append({
                "start": exon_start, "end": exon_end, "strand": f[6],
                "gene": a.get("gene_name", a.get("gene_id", "Unknown")),
                "gene_id": a.get("gene_id", "")
            })

    return exons

def merge_intervals(intervals):
    if not intervals: return []
    intervals = sorted((int(s), int(e)) for s, e in intervals)
    merged = [list(intervals[0])]

    for s, e in intervals[1:]:
        if s <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], e)
        else: merged.append([s, e])

    return merged

def draw_reference_track(ax, exons, plot_start, plot_end):
    genes, strands = defaultdict(list), {}

    for exon in exons:
        gene = exon["gene"] or exon["gene_id"] or "Unknown"
        genes[gene].append((exon["start"], exon["end"]))
        strands[gene] = exon["strand"]

    if not genes:
        ax.text(0.5, 0.45, "No GENCODE v49 exons",
                transform=ax.transAxes, ha="center", va="center", fontsize=7)
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
                ax.add_patch(Rectangle(
                    (s, y - 0.15), e - s, 0.30,
                    facecolor="black", edgecolor="black", linewidth=0.4
                ))

            strand = strands.get(gene, "")
            label = f"{gene} ({strand})" if strand in {"+", "-"} else gene
            ax.text(plot_start, y + 0.22, label, fontsize=6.5, ha="left", va="bottom")

        ax.set_ylim(-0.45, max(0.55, len(gene_items) - 0.25))

    ax.set_xlim(plot_start, plot_end)
    ax.set_yticks([])
    ax.set_ylabel("GENCODE v49", fontsize=8)
    ax.xaxis.set_major_formatter(FuncFormatter(genomic_formatter))
    ax.locator_params(axis="x", nbins=5)

    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)

# ============================================================
# SASHIMI-STYLE ARCS
# ============================================================

def draw_linear_arcs(ax, junctions, x, mean_track, coverage_ymax, color):
    span_total = max(1, x[-1] - x[0])

    for (start, end), reads in sorted(junctions.items(), key=lambda z: z[0][1] - z[0][0]):
        left_y = float(np.interp(start, x, mean_track))
        right_y = float(np.interp(end, x, mean_track))
        span_frac = min(1.0, max(0.0, (end - start) / span_total))
        peak = coverage_ymax * (1.05 + 0.42 * math.sqrt(span_frac))
        mid = (start + end) / 2

        path = MplPath(
            [(start, left_y), (mid, peak), (end, right_y)],
            [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3]
        )
        lw = min(4.0, max(0.7, math.log(reads + 1, 2)))
        ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, linewidth=lw, alpha=0.9))
        ax.text(
            mid, peak, f"{reads:.2f}",
            ha="center", va="bottom", fontsize=7,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.5)
        )

def draw_circ_arc(ax, start, end, mean_reads, coverage_ymax, color):
    if mean_reads <= 0: return
    depth = -0.30 * coverage_ymax
    path = MplPath(
        [(start, 0), (start, depth), (end, depth), (end, 0)],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    )
    lw = min(4.0, max(0.7, math.log(mean_reads + 1, 2)))
    ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, linewidth=lw, alpha=0.95))
    ax.text(
        (start + end) / 2, depth * 0.90, f"{mean_reads:.2f}",
        ha="center", va="top", fontsize=7,
        bbox=dict(facecolor="white", edgecolor="none", pad=0.5)
    )

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

def genotype_group(value, flip=False):
    try: x = float(value)
    except: return None

    if flip: x = 2 - x
    if not np.isfinite(x) or x < -0.1 or x > 2.1: return None

    return int(np.clip(np.rint(x), 0, 2))

# ============================================================
# PROCESS
# ============================================================

hits_by_tissue = defaultdict(list)
for hit in hits:
    hits_by_tissue[hit.get(tissue_col, "").strip()].append(hit)

total_plotted = total_skipped = total_boxplots = 0

for tissue, tissue_hits in hits_by_tissue.items():
    print(f"\n{'=' * 60}\nTISSUE: {tissue}\nPASSING HITS: {len(tissue_hits)}\n{'=' * 60}")

    qtl_dir = ROOT / tissue / TYPE
    results_dir = qtl_dir / "results"

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
    boxplot_pairs = results_dir / f"{tissue}_{TYPE}.top_for_boxplot.txt"
    colored_boxplots = results_dir / f"{tissue}_all_sig_{TYPE}_boxplots_colored.pdf"

    missing_files = [
        p for p in [phenotype_file, trait_location_file, snp_file, snp_location_file]
        if not p.exists()
    ]

    if missing_files:
        print("SKIPPING TISSUE: missing " + ", ".join(map(str, missing_files)))
        total_skipped += len(tissue_hits)
        continue

    events = {h[event_col].strip() for h in tissue_hits if not missing(h.get(event_col))}
    all_snps = {s for h in tissue_hits for s in possible_snps(h)}

    phenotype_rows = extract_matrix_rows(phenotype_file, events)
    genotype_rows = extract_matrix_rows(snp_file, all_snps)
    trait_locations = load_location_table(trait_location_file, "trait")
    snp_locations = load_location_table(snp_location_file, "snp")

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

        if snp not in snp_locations:
            print(f"  SKIP: no genomic position for {snp}"); total_skipped += 1; continue

        snp_chr, snp_pos = snp_locations[snp]
        allele = raw_alleles.get((chr_key(snp_chr), int(snp_pos)))

        if allele is None or allele["counted"] not in {allele["ref"], allele["alt"]}:
            print(f"  SKIP: allele orientation unavailable for {snp}"); total_skipped += 1; continue

        flip_genotype = allele["counted"] == allele["ref"]

        print(
            f"  SNP alleles: REF={allele['ref']} ALT={allele['alt']} "
            f"PLINK_counted={allele['counted']} | "
            f"{'flipping to ALT dosage' if flip_genotype else 'already ALT dosage'}"
        )

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

        reference_exons = load_reference_exons(GTF, chrom, plot_start, plot_end)
        print(f"  GENCODE v49 overlapping exons: {len(reference_exons)}")

        label = symbol if symbol else event
        basename = safe_name(f"{tissue}__{TYPE}__{label}__{snp}")

        boxplot_out = OUTDIR / f"{basename}.boxplot.pdf"
        density_pdf = OUTDIR / f"{basename}.density.pdf"
        density_png = OUTDIR / f"{basename}.density.png"
        density_svg = OUTDIR / f"{basename}.density.svg"
        subject_tsv = OUTDIR / f"{basename}.subjects.tsv"

        if extract_boxplot(colored_boxplots, boxplot_pairs, event, snp, boxplot_out):
            total_boxplots += 1

        records = []

        for qid in [q for q in phenotype if q in genotype]:
            g = genotype_group(genotype[qid], flip=flip_genotype)
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
                "sample_fs_id": sample_dir.name, "genotype": g,
                "genotype_raw": genotype[qid], "phenotype": pheno,
                "sample_dir": sample_dir, "track": None, "track_source": None,
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
            if r["track"] is not None
            and len(r["track"]) == n_bins
            and np.all(np.isfinite(r["track"]))
        ]

        counts = {g: sum(r["genotype"] == g for r in complete) for g in [0, 1, 2]}

        print(f"  Complete subjects used in density plots: {len(complete)}")
        print(f"  Ref/Ref={counts[0]} | Het={counts[1]} | Hom Alt={counts[2]}")

        if not complete:
            print("  SKIP DENSITY: no subjects with usable tracks")
            total_skipped += 1
            continue

        groups = {g: [r for r in complete if r["genotype"] == g] for g in [0, 1, 2]}
        mean_tracks = {
            g: np.nanmean(np.vstack([r["track"] for r in groups[g]]), axis=0)
            if groups[g] else np.zeros(n_bins)
            for g in [0, 1, 2]
        }

        junctions, junction_denoms = {}, {}
        if TYPE in {"sQTL", "cQTL"}:
            for g in [0, 1, 2]:
                junctions[g], junction_denoms[g] = mean_junction_reads(
                    groups[g], chrom, plot_start, plot_end, strand
                )
                if groups[g]:
                    print(
                        f"  {GENOTYPE_LABELS[g]} LeafCutter files: "
                        f"{junction_denoms[g]}/{len(groups[g])}; "
                        f"junctions plotted={len(junctions[g])}"
                    )

        circ_means = {g: 0.0 for g in [0, 1, 2]}
        if TYPE == "cQTL":
            for g in [0, 1, 2]:
                if groups[g]:
                    # Missing circRNA row contributes zero.
                    circ_means[g] = sum(
                        0.0 if r["circ_reads"] is None else float(r["circ_reads"])
                        for r in groups[g]
                    ) / len(groups[g])

            n_quant = sum(r["circ_reads"] is not None for r in complete)
            n_exact = sum(r["circ_match_type"] == "exact" for r in complete)
            n_tol = sum(r["circ_match_type"] == "tolerance" for r in complete)
            print(f"  Circ quant rows found: {n_quant}/{len(complete)}")
            print(f"  Circ coordinate matches: exact={n_exact}, tolerance-rescued={n_tol}")
            print(
                "  Mean circ reads: "
                + " | ".join(f"{GENOTYPE_LABELS[g]}={circ_means[g]:.3f}" for g in [0, 1, 2])
            )

        with open(subject_tsv, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([
                "tissue", "qtl_type", "event", "gene_symbol", "snp", "ref", "alt",
                "plink_counted_allele", "qtl_id", "subject_id", "metadata_sample_id",
                "filesystem_sample_id", "genotype_group", "genotype_label",
                "genotype_raw", "phenotype", "track_source", "circ_reads",
                "linear_reads", "circ_percent_file", "circ_match_type",
                "circ_matched_start", "circ_matched_end", "circ_quant_source"
            ])

            for r in complete:
                writer.writerow([
                    tissue, TYPE, event, symbol, snp, allele["ref"], allele["alt"],
                    allele["counted"], r["qtl_id"], r["subject_id"], r["sample_id"],
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

        if not OVERWRITE and density_pdf.exists() and density_png.exists() and density_svg.exists():
            print("  DENSITY EXISTS; skipping")
            continue

        coverage_ymax = max(float(np.nanmax(mean_tracks[g])) for g in [0, 1, 2] if groups[g])
        coverage_ymax = 1.0 if not np.isfinite(coverage_ymax) or coverage_ymax <= 0 else coverage_ymax * 1.05
        has_linear_arcs = TYPE in {"sQTL", "cQTL"} and any(junctions.get(g) for g in [0, 1, 2])
        has_circ_arcs = TYPE == "cQTL" and any(circ_means[g] > 0 for g in [0, 1, 2])
        upper_ylim = coverage_ymax * (1.58 if has_linear_arcs else 1.08)
        lower_ylim = -0.38 * coverage_ymax if has_circ_arcs else 0
        bar_width = (plot_end - plot_start) / n_bins * 0.92

        # ====================================================
        # DENSITY-ONLY PLOT
        # ============================================================

        fig = plt.figure(figsize=(16, 10))
        gs = fig.add_gridspec(
            4, 1,
            height_ratios=[1.0, 1.0, 1.0, 0.72],
            hspace=0.10
        )

        axes = []

        for g in [0, 1, 2]:
            ax = fig.add_subplot(gs[g, 0])
            axes.append(ax)

            if groups[g]:
                mean_track = mean_tracks[g]
                ax.bar(
                    x, mean_track, width=bar_width,
                    color=TRACK_COLORS[g], alpha=0.82,
                    align="center", linewidth=0
                )

                if TYPE in {"sQTL", "cQTL"}:
                    draw_linear_arcs(
                        ax, junctions[g], x, mean_track,
                        coverage_ymax, TRACK_COLORS[g]
                    )

                if TYPE == "cQTL":
                    draw_circ_arc(
                        ax, trait_start, trait_end,
                        circ_means[g], coverage_ymax,
                        TRACK_COLORS[g]
                    )
            else:
                ax.text(
                    0.5, 0.5, "No subjects",
                    transform=ax.transAxes,
                    ha="center", va="center", fontsize=10
                )

            ax.axvspan(trait_start, trait_end, alpha=0.10, color="#999999")

            if plot_start <= snp_pos <= plot_end:
                ax.axvline(
                    snp_pos, linestyle="--",
                    linewidth=1.3, alpha=0.8, color="black"
                )

            ax.set_xlim(plot_start, plot_end)
            ax.set_ylim(lower_ylim, upper_ylim)
            ax.set_title(
                f"{GENOTYPE_LABELS[g]} (N={counts[g]})",
                loc="left", fontsize=11, color=TRACK_COLORS[g]
            )

            if g < 2:
                ax.tick_params(axis="x", labelbottom=False)
            else:
                ax.tick_params(axis="x", labelbottom=False)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        # ---------- ONE SHARED GENCODE TRACK ----------

        ax_ref = fig.add_subplot(gs[3, 0], sharex=axes[0])
        draw_reference_track(ax_ref, reference_exons, plot_start, plot_end)
        ax_ref.axvspan(trait_start, trait_end, alpha=0.10, color="#999999")

        if plot_start <= snp_pos <= plot_end:
            ax_ref.axvline(
                snp_pos, linestyle="--",
                linewidth=1.0, alpha=0.7, color="black"
            )

        ax_ref.set_xlabel(f"{chrom} genomic position")

        title_event = f"{symbol} ({event})" if symbol and symbol != event else event
        title = (
            f"{tissue_display(tissue)} | {TYPE} | {title_event} | {snp} | "
            f"REF={allele['ref']} ALT={allele['alt']}"
        )

        if pp_h4 and not missing(pp_h4):
            try: title += f" | PP.H4={float(pp_h4):.4f}"
            except: title += f" | PP.H4={pp_h4}"

        fig.suptitle(title, fontsize=13, y=0.995)
        fig.text(
            0.012, 0.57,
            "Mean normalized RNA-seq read density",
            rotation=90, va="center", fontsize=10
        )
        fig.text(
            0.5, 0.008,
            f"Trait: {chrom}:{trait_start:,}-{trait_end:,} | "
            f"flank: {flank:,} bp each side | shaded region = trait locus",
            ha="center", fontsize=9
        )

        fig.savefig(density_pdf, bbox_inches="tight")
        fig.savefig(density_png, dpi=DPI, bbox_inches="tight")
        fig.savefig(density_svg, bbox_inches="tight")
        plt.close(fig)

        print(f"  BOXPLOT PDF : {boxplot_out if boxplot_out.exists() else 'not extracted'}")
        print(f"  DENSITY PDF : {density_pdf}")
        print(f"  DENSITY PNG : {density_png}")
        print(f"  DENSITY SVG : {density_svg}")
        print(f"  SUBJECT TSV : {subject_tsv}")
        total_plotted += 1

print(f"\n{'=' * 60}")
print("DONE")
print(f"QTL type           : {TYPE}")
print(f"Density plots      : {total_plotted}")
print(f"Boxplots extracted : {total_boxplots}")
print(f"Skipped hits       : {total_skipped}")
print(f"Output             : {OUTDIR}")
print("=" * 60)
PY
