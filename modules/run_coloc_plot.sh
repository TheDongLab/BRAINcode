#!/bin/bash
#SBATCH --job-name=run_coloc_plot
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/run_coloc_plot.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/run_coloc_plot.err
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G

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
MIN_MEAN_JUNCTION_RPM=0.01
MAX_JUNCTIONS=30
MAX_HITS=0
TISSUE_FILTER=""
OVERWRITE=0
DPI=180

# ============================================================

set -euo pipefail
case "$TYPE" in eQTL|sQTL|cQTL) ;; *) echo "Usage: sbatch $0 {eQTL|sQTL|cQTL}"; exit 1 ;; esac

COLOC="${MR_DIR}/${TYPE}_SuSiE_coloc_summary.tsv"
OUTDIR="${MR_DIR}/coloc_raw_plots/${TYPE}"
mkdir -p "$OUTDIR"

for f in "$COLOC" "$METADATA" "$GTF" "$RAW_FILE"; do [[ -f "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done

module load deepTools 2>/dev/null || { echo "ERROR: could not load deepTools"; exit 1; }
module load poppler/25.07.0-GCC-13.3.0 2>/dev/null || { echo "ERROR: could not load Poppler"; exit 1; }

PYTHON="$(command -v python)"
command -v pdfseparate >/dev/null || { echo "ERROR: pdfseparate unavailable"; exit 1; }
command -v pdftotext >/dev/null || { echo "ERROR: pdftotext unavailable"; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools unavailable"; exit 1; }
"$PYTHON" -c 'import numpy,matplotlib,pyBigWig' || { echo "ERROR: Python plotting stack unavailable"; exit 1; }

export TYPE ROOT METADATA COLOC OUTDIR GTF RAW_FILE
export FLANK_FRAC FLANK_MIN FLANK_MAX BIGWIG_BIN_BP MAX_TRACK_POINTS CIRC_COORD_TOL
export MIN_MEAN_JUNCTION_RPM MAX_JUNCTIONS MAX_HITS TISSUE_FILTER OVERWRITE DPI

"$PYTHON" - <<'PY'
import os,re,csv,sys,math,warnings,subprocess
from pathlib import Path
from collections import defaultdict

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    from matplotlib.patches import Rectangle,PathPatch
    from matplotlib.path import Path as MplPath
    import pyBigWig
except ImportError as e:
    sys.exit(f"ERROR: missing Python package: {e}")

warnings.filterwarnings("ignore")

TYPE,ROOT=os.environ["TYPE"],Path(os.environ["ROOT"])
METADATA,COLOC,OUTDIR,GTF,RAW_FILE=map(Path,[os.environ["METADATA"],os.environ["COLOC"],os.environ["OUTDIR"],os.environ["GTF"],os.environ["RAW_FILE"]])
FLANK_FRAC,FLANK_MIN,FLANK_MAX=float(os.environ["FLANK_FRAC"]),int(os.environ["FLANK_MIN"]),int(os.environ["FLANK_MAX"])
BIGWIG_BIN_BP,MAX_TRACK_POINTS=int(os.environ["BIGWIG_BIN_BP"]),int(os.environ["MAX_TRACK_POINTS"])
CIRC_COORD_TOL,MIN_MEAN_JUNCTION_RPM=int(os.environ["CIRC_COORD_TOL"]),float(os.environ["MIN_MEAN_JUNCTION_RPM"])
MAX_JUNCTIONS,MAX_HITS=int(os.environ["MAX_JUNCTIONS"]),int(os.environ["MAX_HITS"])
TISSUE_FILTER,OVERWRITE,DPI=os.environ["TISSUE_FILTER"].strip(),int(os.environ["OVERWRITE"]),int(os.environ["DPI"])

GENOTYPE_LABELS={0:"Ref/Ref",1:"Het",2:"Hom Alt"}
TRACK_COLORS={0:"#228B22",1:"#4682B4",2:"#C75B7A"}

# ============================================================
# HELPERS + TARGETALS TISSUE REMAPPING
# ============================================================

def norm(s): return re.sub(r"[^a-z0-9]","",str(s).lower())
def chr_key(s): return re.sub(r"^chr","",str(s),flags=re.I).upper()
def safe_name(s): return re.sub(r"[^A-Za-z0-9._+-]+","_",str(s))[:180]
def strip_gene_version(x): return re.sub(r"\.\d+$","",str(x))
def tissue_display(x): return str(x).replace("_"," ")
def missing(v): return v is None or str(v).strip().lower() in {"","na","nan","none",".","null"}
def truthy(v): return str(v).strip().lower() in {"true","t","1","yes","y","pass","passed"}
def detect_delimiter(line): return "\t" if "\t" in line else None
def split_line(line,delim): return line.rstrip("\r\n").split("\t") if delim=="\t" else line.split()

TISSUE_REMAP={
    "Motor Cortex Lateral":"Motor_Cortex",
    "Motor Cortex Medial":"Motor_Cortex",
    "Lateral Motor Cortex":"Motor_Cortex",
    "Medial Motor Cortex":"Motor_Cortex",
    "Primary Motor Cortex L":"Motor_Cortex",
    "Primary Motor Cortex M":"Motor_Cortex",
    "Cortex_Motor_Unspecified":"Motor_Cortex",
    "Cortex_Motor_BA4":"Motor_Cortex",
    "BA4 Motor Cortex":"Motor_Cortex",
    "Lateral_motor_cortex":"Motor_Cortex",
    "Frontal Cortex":"Frontal_Cortex",
    "Cerebellum":"Cerebellum",
    "Spinal_Cord_Cervical":"Cervical_Spinal_Cord",
    "Cervical Spinal Cord":"Cervical_Spinal_Cord",
    "Cervical_spinal_cord":"Cervical_Spinal_Cord",
    "Spinal_cord_Cervical":"Cervical_Spinal_Cord",
    "Lumbar Spinal Cord":"Lumbar_Spinal_Cord",
    "Thoracic Spinal Cord":"Thoracic_Spinal_Cord",
    "Spinal_Cord_Lumbosacral":"Lumbar_Spinal_Cord",
    "Lumbosacral_Spinal_Cord":"Lumbar_Spinal_Cord",
    "Lumbar_spinal_cord":"Lumbar_Spinal_Cord"
}
TISSUE_REMAP_NORM={norm(k):v for k,v in TISSUE_REMAP.items()}
CANONICAL_TISSUES={"Motor_Cortex","Frontal_Cortex","Cerebellum","Cervical_Spinal_Cord","Lumbar_Spinal_Cord","Thoracic_Spinal_Cord"}
for t in CANONICAL_TISSUES:TISSUE_REMAP_NORM[norm(t)]=t

def canonical_tissue(x):
    x=str(x).strip()
    return TISSUE_REMAP_NORM.get(norm(x),x)

def tissue_equal(a,b):
    return canonical_tissue(a)==canonical_tissue(b)

def format_max_sig(value,sig=5):
    if value==0:return "0"
    decimals=max(0,sig-int(math.floor(math.log10(abs(value))))-1)
    return f"{value:.{decimals}f}".rstrip("0").rstrip(".")

def genomic_formatter(x,pos=None):
    if abs(x)>=1_000_000:return f"{format_max_sig(x/1_000_000)} Mb"
    if abs(x)>=1_000:return f"{format_max_sig(x/1_000)} kb"
    return format_max_sig(x)

def density_formatter(y,pos=None):
    if abs(y)<1e-12:return "0"
    return format_max_sig(abs(y),sig=4)

def find_header(headers,aliases,required=True):
    lookup={norm(h):h for h in headers}
    for a in aliases:
        if norm(a) in lookup:return lookup[norm(a)]
    for h in headers:
        if any(norm(a) and norm(a) in norm(h) for a in aliases):return h
    if required:raise KeyError(f"Could not find {aliases}. Available: {headers}")
    return None

# ============================================================
# PLINK ALLELES
# ============================================================

raw_alleles={}
with open(RAW_FILE) as fh:
    for col in fh.readline().split():
        m=re.match(r"^(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)_([ACGT]+)$",col,re.I)
        if m:
            chrom,pos,ref,alt,counted=m.groups()
            raw_alleles[(chr_key(chrom),int(pos))]={"ref":ref.upper(),"alt":alt.upper(),"counted":counted.upper()}
print(f"PLINK RAW allele mappings: {len(raw_alleles):,}")

# ============================================================
# COLOC SUMMARY
# ============================================================

with open(COLOC,newline="") as fh:
    reader=csv.DictReader(fh,delimiter="\t"); summary_headers=reader.fieldnames or []
    if "H4_reference_pass" not in summary_headers:sys.exit("ERROR: summary lacks H4_reference_pass")
    hits=[r for r in reader if truthy(r.get("H4_reference_pass"))]

print("="*60); print(f"QTL TYPE        : {TYPE}"); print(f"H4 PASSING HITS : {len(hits)}"); print("FILTER          : H4_reference_pass == TRUE"); print("="*60)
if TISSUE_FILTER:hits=[h for h in hits if tissue_equal(h.get("tissue",""),TISSUE_FILTER)]
if MAX_HITS>0:hits=hits[:MAX_HITS]
if not hits:sys.exit("No passing coloc hits to plot.")

tissue_col=find_header(summary_headers,["tissue"])
if TYPE=="eQTL":
    event_col=find_header(summary_headers,["geneid","gene_id","gene"]); symbol_col=find_header(summary_headers,["gene_symbol","symbol"],required=False)
elif TYPE=="sQTL":
    event_col=find_header(summary_headers,["junction_id","junctionid","splicing_id","splice_id","event_id","phenotype_id"]); symbol_col=find_header(summary_headers,["gene_symbol","symbol","gene"],required=False)
else:
    event_col=find_header(summary_headers,["circ_id","circid","circrna_id","event_id","phenotype_id"]); symbol_col=find_header(summary_headers,["gene_symbol","symbol","gene"],required=False)

# ============================================================
# METADATA
# ============================================================

with open(METADATA,newline="") as fh:
    reader=csv.DictReader(fh); meta_headers,metadata=reader.fieldnames or [],list(reader)

subject_col=find_header(meta_headers,["externalsubjectid","subject_id","subjectid"])
sample_col=find_header(meta_headers,["externalsampleid","sample_id","sampleid"])
meta_tissue_col=find_header(meta_headers,["tissue"],required=False)
rin_col=find_header(meta_headers,["rin","rna_integrity_number"],required=False)

def get_rin(row):
    try:return float(row.get(rin_col,"")) if rin_col else -1.
    except:return -1.

def metadata_for_qtl_id(qtl_id,tissue):
    qtl_id,rows=str(qtl_id),metadata
    target_tissue=canonical_tissue(tissue)
    if meta_tissue_col:
        tissue_rows=[r for r in metadata if canonical_tissue(r.get(meta_tissue_col,""))==target_tissue]
        if tissue_rows:rows=tissue_rows
    exact=[r for r in rows if str(r.get(sample_col,"")).strip()==qtl_id]
    if exact:return max(exact,key=get_rin)
    subject_rows=[r for r in rows if str(r.get(subject_col,"")).strip()==qtl_id]
    return max(subject_rows,key=get_rin) if subject_rows else None

# ============================================================
# MATRIX READER
# ============================================================

def extract_matrix_rows(path,wanted):
    wanted,result=set(map(str,wanted)),{}
    wanted_nover={strip_gene_version(x) for x in wanted}
    with open(path) as fh:
        first=fh.readline(); delim=detect_delimiter(first); samples=split_line(first,delim)[1:]
        for line in fh:
            if not line.strip():continue
            fields=split_line(line,delim)
            if len(fields)<2:continue
            rid=fields[0].strip()
            if rid not in wanted and strip_gene_version(rid) not in wanted_nover:continue
            vals=fields[1:]+[""]*max(0,len(samples)-len(fields[1:]))
            result[rid]=dict(zip(samples,vals[:len(samples)]))
    return result

def find_matrix_row(rows,event):
    if event in rows:return rows[event]
    target=strip_gene_version(event)
    return next((v for rid,v in rows.items() if strip_gene_version(rid)==target),None)

# ============================================================
# LOCATION TABLES
# ============================================================

def load_location_table(path,kind):
    with open(path) as fh:
        first=fh.readline(); delim=detect_delimiter(first); fields=split_line(first,delim)
        is_data=len(fields)>=3 and re.match(r"^(chr)?[0-9XYM]+$",fields[1],re.I) and re.match(r"^\d+$",fields[2])
        if is_data:
            headers=["id","chr","start"]
            if len(fields)>=4:headers.append("end")
            if len(fields)>=5:headers.append("strand")
            while len(headers)<len(fields):headers.append(f"extra{len(headers)}")
            lines=[fields]
        else:headers,lines=fields,[]
        lines += [split_line(line,delim) for line in fh if line.strip()]

    rows=[]
    for f in lines:
        f += [""]*max(0,len(headers)-len(f)); rows.append(dict(zip(headers,f)))

    if kind=="snp":
        idc=find_header(headers,["snpid","snp_id","rsid","id"]); cc=find_header(headers,["chr","chrom","chromosome"]); pc=find_header(headers,["pos","position","start"])
        return {str(r[idc]):(str(r[cc]),int(float(r[pc]))) for r in rows if not missing(r.get(idc)) and not missing(r.get(cc)) and not missing(r.get(pc))}

    aliases={"eQTL":["geneid","gene_id","gene","id"],"sQTL":["junction_id","junctionid","splicing_id","event_id","phenotype_id","id"],"cQTL":["circ_id","circid","circrna_id","event_id","phenotype_id","id"]}
    idc=find_header(headers,aliases[TYPE]); cc=find_header(headers,["chr","chrom","chromosome"])
    sc=find_header(headers,["start","left","begin","donor"]); ec=find_header(headers,["end","right","stop","acceptor"])
    strandc=find_header(headers,["strand"],required=False); out={}
    for r in rows:
        try:out[str(r[idc])]={"chr":str(r[cc]),"start":int(float(r[sc])),"end":int(float(r[ec])),"strand":str(r[strandc]) if strandc and not missing(r.get(strandc)) else None}
        except:pass
    return out

def parse_event_coordinates(event):
    for pat in [r"(chr[^:]+):(\d+)-(\d+):([+-])",r"(chr[^:]+):(\d+):(\d+):([+-])",r"(chr[^:]+):(\d+)-(\d+)",r"(chr[^:]+):(\d+):(\d+)"]:
        m=re.search(pat,str(event))
        if m:return {"chr":m.group(1),"start":int(m.group(2)),"end":int(m.group(3)),"strand":m.group(4) if len(m.groups())>=4 else None}
    return None

def find_location(locations,event):
    if event in locations:return locations[event]
    target=strip_gene_version(event)
    return next((loc for rid,loc in locations.items() if strip_gene_version(rid)==target),parse_event_coordinates(event))

# ============================================================
# COLORED BOXPLOT EXTRACTION
# ============================================================

def naive_boxplot_page(pairs_file,event,snp):
    if not pairs_file.exists():return None
    awk_script='NR>1 && $1==E && $NF==S {print NR-1; exit}' if TYPE=="eQTL" else '$1==E && $NF==S {print NR; exit}'
    try:
        r=subprocess.run(["awk","-v",f"E={event}","-v",f"S={snp}",awk_script,str(pairs_file)],
                         capture_output=True,text=True,check=True)
        return int(r.stdout.strip()) if r.stdout.strip() else None
    except:return None

def pdf_match_text(x):
    return re.sub(r"[^a-z0-9]","",str(x).lower())

def search_pdf_window(source_pdf,event,snp,first,last):
    target_event,target_snp=pdf_match_text(event),pdf_match_text(snp)
    try:
        r=subprocess.run(["pdftotext","-f",str(first),"-l",str(last),str(source_pdf),"-"],
                         capture_output=True,text=True,check=True)
    except:return None
    for i,text in enumerate(r.stdout.split("\f")):
        if not text.strip():continue
        cleaned=pdf_match_text(text)
        if target_event in cleaned and target_snp in cleaned:return first+i
    return None

def find_boxplot_page(source_pdf,pairs_file,event,snp):
    naive=naive_boxplot_page(pairs_file,event,snp)
    if naive is None:return None

    page=search_pdf_window(source_pdf,event,snp,naive,naive)
    if page is not None:
        print(f"  BOXPLOT LOOKUP: pair-file page={naive}, actual PDF page={page}, offset=+0")
        return page

    for radius in [10,50,250,1000,5000]:
        first,last=max(1,naive-radius),naive+radius
        page=search_pdf_window(source_pdf,event,snp,first,last)
        if page is not None:
            print(f"  BOXPLOT LOOKUP: pair-file page={naive}, actual PDF page={page}, offset={page-naive:+d}")
            return page

    print(f"  BOXPLOT LOOKUP FAILED: {event} + {snp}; pair-file page={naive}")
    return None

def extract_boxplot(source_pdf,pairs_file,event,snp,output_pdf):
    if not source_pdf.exists():print(f"  BOXPLOT FAILED: colored PDF missing: {source_pdf}"); return False
    if not pairs_file.exists():print(f"  BOXPLOT FAILED: pairs file missing: {pairs_file}"); return False

    page=find_boxplot_page(source_pdf,pairs_file,event,snp)
    if page is None:
        print(f"  BOXPLOT FAILED: could not locate actual PDF page for {event} + {snp}")
        return False

    try:
        pattern=output_pdf.with_name(output_pdf.stem+"_%d.pdf")
        tmp=Path(str(pattern).replace("%d",str(page)))
        if tmp.exists():tmp.unlink()

        subprocess.run(["pdfseparate","-f",str(page),"-l",str(page),str(source_pdf),str(pattern)],
                       check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,text=True)

        if not tmp.exists():
            print(f"  BOXPLOT FAILED: page {page} extraction failed")
            return False

        if output_pdf.exists():output_pdf.unlink()
        tmp.replace(output_pdf)
        print(f"  BOXPLOT: actual page {page} -> {output_pdf}")
        return True
    except Exception as e:
        print(f"  BOXPLOT FAILED: {e}")
        return False

# ============================================================
# SAMPLE DIRECTORY + PRIMARY READ DEPTH
# ============================================================

processed_root_cache={}
processed_sample_cache={}
primary_read_cache={}

def processed_roots_for_tissue(tissue):
    canonical=canonical_tissue(tissue)
    if canonical in processed_root_cache:return processed_root_cache[canonical]

    roots=[]
    preferred=ROOT/canonical/"RNAseq"/"Processed"
    if preferred.is_dir():roots.append(preferred)

    # targetALS contains duplicate/legacy tissue directory spellings.
    # Search every top-level tissue directory whose name remaps to the
    # same canonical tissue, e.g.:
    #   Lumbar_Spinal_Cord/
    #   Lumbar_spinal_cord/
    if ROOT.is_dir():
        for p in ROOT.iterdir():
            if not p.is_dir() or canonical_tissue(p.name)!=canonical:continue
            candidate=p/"RNAseq"/"Processed"
            if candidate.is_dir() and candidate not in roots:roots.append(candidate)

    processed_root_cache[canonical]=roots
    print(f"  RNAseq roots for {canonical}: "+(" | ".join(map(str,roots)) if roots else "NONE"))
    return roots

def resolve_sample_dir(tissue,sample_id):
    canonical=canonical_tissue(tissue)
    cache_key=(canonical,norm(sample_id))
    if cache_key in processed_sample_cache:return processed_sample_cache[cache_key]

    roots=processed_roots_for_tissue(canonical)

    # First try exact/simple hyphen↔underscore forms in every valid tissue root.
    for processed in roots:
        for candidate in [sample_id,sample_id.replace("-","_"),sample_id.replace("_","-")]:
            p=processed/candidate
            if p.is_dir():
                processed_sample_cache[cache_key]=p
                return p

    # Then normalized lookup, still restricted to equivalent tissue roots.
    target=norm(sample_id)
    for processed in roots:
        key=str(processed)
        if key not in processed_sample_cache:
            processed_sample_cache[key]={norm(p.name):p for p in processed.iterdir() if p.is_dir()}
        p=processed_sample_cache[key].get(target)
        if p is not None:
            processed_sample_cache[cache_key]=p
            return p

    processed_sample_cache[cache_key]=None
    return None

def primary_bam(sample_dir): return sample_dir/"STAR.Aligned.sortedByCoord.out.bam"

def primary_read_count(sample_dir):
    bam=primary_bam(sample_dir); key=str(bam)
    if key in primary_read_cache:return primary_read_cache[key]
    if not bam.exists():
        primary_read_cache[key]=None; return None
    try:
        n=int(subprocess.check_output(["samtools","view","-F","0x100","-c",str(bam)],text=True).strip())
        primary_read_cache[key]=n if n>0 else None
    except Exception:
        primary_read_cache[key]=None
    return primary_read_cache[key]

def rpm_value(reads,primary_reads):
    return float(reads)*1_000_000.0/primary_reads if primary_reads and primary_reads>0 else None

# ============================================================
# STRAND-SPECIFIC NORMALIZED BIGWIG
# ============================================================

def bw_chrom(bw,chrom):
    chroms=bw.chroms(); candidates=[chrom,chrom[3:] if chrom.startswith("chr") else "chr"+chrom]
    return next((c for c in candidates if c in chroms),None)

def bigwig_strands(sample_dir,chrom,start,end,n_bins):
    plus=sample_dir/"STAR.Aligned.sortedByCoord.out.plus.normalized.bw"
    minus=sample_dir/"STAR.Aligned.sortedByCoord.out.minus.normalized.bw"

    def read_one(path):
        if not path.exists():return None
        bw=pyBigWig.open(str(path))
        try:
            c=bw_chrom(bw,chrom)
            if c is None:return None
            vals=bw.stats(c,int(start),int(end),nBins=int(n_bins),type="mean")
            return np.abs(np.asarray([0. if v is None else float(v) for v in vals]))
        finally:bw.close()

    return read_one(plus),read_one(minus),str(plus),str(minus)

# ============================================================
# LEAFCUTTER JUNCTIONS — PER-SAMPLE RPM
# ============================================================

leafcutter_cache={}
def leafcutter_file(sample_dir): return sample_dir/"leafcutter"/"psi"/f"{sample_dir.name}.leafcutter.PSI.tsv"

def load_leafcutter(path):
    key=str(path)
    if key in leafcutter_cache:return leafcutter_cache[key]
    if not path.exists():leafcutter_cache[key]=None; return None
    junctions={}
    try:
        with open(path) as fh:
            reader=csv.DictReader(fh,delimiter="\t")
            for row in reader:
                try:
                    k=(chr_key(row["chrom"]),str(row["strand"]).strip(),int(float(row["start"])),int(float(row["end"])))
                    junctions[k]=junctions.get(k,0.)+float(row["reads"])
                except (KeyError,ValueError,TypeError):continue
    except:junctions=None
    leafcutter_cache[key]=junctions
    return junctions

def mean_junction_rpm(group,chrom,start,end):
    sums,n_files=defaultdict(float),0
    for r in group:
        junctions=load_leafcutter(leafcutter_file(r["sample_dir"]))
        primary=r["primary_reads"]
        if junctions is None or not primary:continue
        n_files+=1
        for (jchr,jstrand,js,je),reads in junctions.items():
            if jchr!=chr_key(chrom) or js<start or je>end:continue
            sums[(jstrand,js,je)] += rpm_value(reads,primary)

    if n_files==0:return {},0
    means={k:v/n_files for k,v in sums.items()}
    means={k:v for k,v in means.items() if v>=MIN_MEAN_JUNCTION_RPM}
    if MAX_JUNCTIONS>0 and len(means)>MAX_JUNCTIONS:
        keep=sorted(means,key=means.get,reverse=True)[:MAX_JUNCTIONS]; means={k:means[k] for k in keep}
    return means,n_files

# ============================================================
# cQTL QUANTIFICATION
# ============================================================

def circ_quant_file(sample_dir):return sample_dir/"circularRNA_known_circ_percentage.txt"

def read_circ_quant(path,chrom,start,end,strand,tolerance=CIRC_COORD_TOL):
    if not path.exists():return None
    exact,nearby=None,[]
    try:
        with open(path) as fh:
            reader=csv.DictReader(fh,delimiter="\t")
            for row in reader:
                try:
                    rchr=row["chrom"]; rstart,rend=int(float(row["start"])),int(float(row["end"])); rstrand=str(row["strand"]).strip()
                    if chr_key(rchr)!=chr_key(chrom):continue
                    if strand in {"+","-"} and rstrand!=str(strand).strip():continue
                    result={"circ_reads":float(row["circ_reads"]),"linear_reads":float(row["linear_reads"]),"circ_percent":float(row["circ_percent"]),
                            "matched_start":rstart,"matched_end":rend,"match_type":None}
                    if rstart==int(start) and rend==int(end):
                        result["match_type"]="exact"; exact=result; break
                    if abs(rstart-int(start))<=tolerance and abs(rend-int(end))<=tolerance:
                        result["match_type"]="tolerance"; nearby.append(result)
                except (ValueError,KeyError,TypeError):continue
    except Exception as e:
        print(f"    WARNING: could not read circ quant file {path}: {e}"); return None

    if exact is not None:return exact
    if not nearby:return None
    nearby.sort(key=lambda r:(abs(r["matched_start"]-int(start))+abs(r["matched_end"]-int(end)),abs(r["matched_start"]-int(start)),abs(r["matched_end"]-int(end))))
    return nearby[0]

# ============================================================
# GENCODE
# ============================================================

def parse_gtf_attrs(s):
    attrs={}
    for item in s.strip().split(";"):
        item=item.strip()
        if not item:continue
        p=item.split(" ",1)
        if len(p)==2:attrs[p[0]]=p[1].strip().strip('"')
    return attrs

def load_reference_exons(gtf,chrom,start,end):
    exons=[]
    with open(gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"):continue
            f=line.rstrip("\n").split("\t")
            if len(f)<9 or f[2]!="exon" or chr_key(f[0])!=chr_key(chrom):continue
            try:exon_start,exon_end=int(f[3])-1,int(f[4])
            except ValueError:continue
            if exon_end<start or exon_start>end:continue
            a=parse_gtf_attrs(f[8])
            exons.append({"start":exon_start,"end":exon_end,"strand":f[6],"gene":a.get("gene_name",a.get("gene_id","Unknown")),"gene_id":a.get("gene_id","")})
    return exons

def merge_intervals(intervals):
    if not intervals:return []
    intervals=sorted((int(s),int(e)) for s,e in intervals); merged=[list(intervals[0])]
    for s,e in intervals[1:]:
        if s<=merged[-1][1]:merged[-1][1]=max(merged[-1][1],e)
        else:merged.append([s,e])
    return merged

def draw_reference_track(ax,exons,plot_start,plot_end):
    genes,strands=defaultdict(list),{}
    for exon in exons:
        gene=exon["gene"] or exon["gene_id"] or "Unknown"
        genes[gene].append((exon["start"],exon["end"])); strands[gene]=exon["strand"]

    if not genes:
        ax.text(.5,.45,"No GENCODE v49 exons",transform=ax.transAxes,ha="center",va="center",fontsize=9)
    else:
        gene_items=sorted(genes.items(),key=lambda kv:min(s for s,e in kv[1]))
        for y,(gene,intervals) in enumerate(gene_items):
            merged=merge_intervals(intervals)
            span_start=max(plot_start,min(s for s,e in merged)); span_end=min(plot_end,max(e for s,e in merged))
            ax.plot([span_start,span_end],[y,y],color="black",linewidth=.8)
            for s,e in merged:
                s,e=max(s,plot_start),min(e,plot_end)
                if e>s:ax.add_patch(Rectangle((s,y-.15),e-s,.30,facecolor="black",edgecolor="black",linewidth=.4))
            strand=strands.get(gene,""); label=f"{gene} ({strand})" if strand in {"+","-"} else gene
            ax.text(plot_start,y+.22,label,fontsize=8,ha="left",va="bottom")
        ax.set_ylim(-.45,max(.55,len(gene_items)-.25))

    ax.set_xlim(plot_start,plot_end); ax.set_yticks([])
    ax.set_ylabel("GENCODE v49",fontsize=12,labelpad=13)
    ax.xaxis.set_major_formatter(FuncFormatter(genomic_formatter)); ax.locator_params(axis="x",nbins=5); ax.tick_params(axis="x",labelsize=11)
    for side in ["top","right","left"]:ax.spines[side].set_visible(False)

# ============================================================
# STRAND-AWARE SASHIMI ARCS
# ============================================================

def draw_linear_arcs(ax,junctions,x,mean_plus,mean_minus,coverage_max,color):
    span_total=max(1,x[-1]-x[0])
    for (strand,start,end),rpm in sorted(junctions.items(),key=lambda z:z[0][2]-z[0][1]):
        span_frac=min(1.,max(0.,(end-start)/span_total)); mid=(start+end)/2

        if strand=="-":
            left_y=-float(np.interp(start,x,mean_minus)); right_y=-float(np.interp(end,x,mean_minus))
            base=min(left_y,right_y); apex=-coverage_max*(1.05+.42*math.sqrt(span_frac))
            path=MplPath([(start,left_y),(mid,apex),(end,right_y)],[MplPath.MOVETO,MplPath.CURVE3,MplPath.CURVE3])
            label_y=base+.58*(apex-base)
        else:
            left_y=float(np.interp(start,x,mean_plus)); right_y=float(np.interp(end,x,mean_plus))
            base=max(left_y,right_y); apex=coverage_max*(1.05+.42*math.sqrt(span_frac))
            path=MplPath([(start,left_y),(mid,apex),(end,right_y)],[MplPath.MOVETO,MplPath.CURVE3,MplPath.CURVE3])
            label_y=base+.58*(apex-base)

        lw=min(4.,max(.7,math.log(rpm+1,2)))
        ax.add_patch(PathPatch(path,facecolor="none",edgecolor=color,linewidth=lw,alpha=.9))
        ax.text(mid,label_y,f"{rpm:.3f}",ha="center",va="center",fontsize=7,
                bbox=dict(facecolor="white",edgecolor="none",alpha=.88,pad=.45))

def draw_circ_arc(ax,start,end,mean_rpm,coverage_max,color,strand):
    if mean_rpm<=0:return
    sign=-1 if strand=="-" else 1
    outer=sign*coverage_max*1.72
    path=MplPath([(start,0),(start,outer),(end,outer),(end,0)],[MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4])
    lw=min(4.,max(.7,math.log(mean_rpm+1,2)))
    ax.add_patch(PathPatch(path,facecolor="none",edgecolor=color,linewidth=lw,alpha=.95))
    ax.text((start+end)/2,outer*.63,f"{mean_rpm:.3f}",ha="center",va="center",fontsize=7,
            bbox=dict(facecolor="white",edgecolor="none",alpha=.88,pad=.45))

# ============================================================
# SNP SELECTION
# ============================================================

def possible_snps(hit):
    candidates=[]
    for col in ["best_candidate_SNP","best_candidate_snp"]:
        if col in hit and not missing(hit[col]):candidates.append(hit[col])
    for col in {"eQTL":["top_eQTL_SNP","top_QTL_SNP"],"sQTL":["top_sQTL_SNP","top_QTL_SNP"],"cQTL":["top_cQTL_SNP","top_QTL_SNP"]}[TYPE]:
        if col in hit and not missing(hit[col]):candidates.append(hit[col])
    for k,v in hit.items():
        if "topsnp" in norm(k) and "gwas" not in norm(k) and not missing(v):candidates.append(v)
    if "candidate_SNPs" in hit and not missing(hit["candidate_SNPs"]):candidates += [x for x in hit["candidate_SNPs"].split(";") if x]
    out=[]
    for x in map(str.strip,candidates):
        if x and x not in out:out.append(x)
    return out

def genotype_group(value,flip=False):
    try:x=float(value)
    except:return None
    if flip:x=2-x
    if not np.isfinite(x) or x<-.1 or x>2.1:return None
    return int(np.clip(np.rint(x),0,2))

# ============================================================
# PROCESS
# ============================================================

hits_by_tissue=defaultdict(list)
for hit in hits:hits_by_tissue[canonical_tissue(hit.get(tissue_col,"").strip())].append(hit)

total_plotted=total_skipped=total_boxplots=0

for tissue,tissue_hits in hits_by_tissue.items():
    print(f"\n{'='*60}\nTISSUE: {tissue}\nPASSING HITS: {len(tissue_hits)}\n{'='*60}")
    qtl_dir=ROOT/tissue/TYPE; results_dir=qtl_dir/"results"

    if TYPE=="eQTL":phenotype_file,trait_location_file=qtl_dir/f"expression_{tissue}.txt",qtl_dir/"gene_location.txt"
    elif TYPE=="sQTL":phenotype_file,trait_location_file=qtl_dir/f"splicing_{tissue}.txt",qtl_dir/"splicing_location.txt"
    else:phenotype_file,trait_location_file=qtl_dir/f"circ_{tissue}.txt",qtl_dir/"circ_location.txt"

    snp_file,snp_location_file=qtl_dir/f"snp_{tissue}.txt",qtl_dir/"snp_location.txt"
    boxplot_pairs=results_dir/f"{tissue}_{TYPE}.top_for_boxplot.txt"
    colored_boxplots=results_dir/f"{tissue}_all_sig_{TYPE}_boxplots_colored.pdf"

    missing_files=[p for p in [phenotype_file,trait_location_file,snp_file,snp_location_file] if not p.exists()]
    if missing_files:
        print("SKIPPING TISSUE: missing "+", ".join(map(str,missing_files))); total_skipped+=len(tissue_hits); continue

    # Show all actual RNA-seq roots that will be searched for this canonical tissue.
    processed_roots_for_tissue(tissue)

    events={h[event_col].strip() for h in tissue_hits if not missing(h.get(event_col))}
    all_snps={s for h in tissue_hits for s in possible_snps(h)}
    phenotype_rows=extract_matrix_rows(phenotype_file,events); genotype_rows=extract_matrix_rows(snp_file,all_snps)
    trait_locations=load_location_table(trait_location_file,"trait"); snp_locations=load_location_table(snp_location_file,"snp")

    for index,hit in enumerate(tissue_hits,1):
        event,pp_h4=hit[event_col].strip(),hit.get("max_PP_H4",""); symbol=hit.get(symbol_col,"").strip() if symbol_col else ""
        print(f"\n[{index}/{len(tissue_hits)}] {tissue} | {event} | {symbol}")

        phenotype=find_matrix_row(phenotype_rows,event)
        if phenotype is None:print("  SKIP: phenotype missing"); total_skipped+=1; continue

        snp=genotype=None
        for candidate in possible_snps(hit):
            if candidate in genotype_rows:snp,genotype=candidate,genotype_rows[candidate]; break
        if genotype is None:print("  SKIP: candidate SNP absent from genotype matrix"); total_skipped+=1; continue
        if snp not in snp_locations:print(f"  SKIP: no genomic position for {snp}"); total_skipped+=1; continue

        snp_chr,snp_pos=snp_locations[snp]; allele=raw_alleles.get((chr_key(snp_chr),int(snp_pos)))
        if allele is None or allele["counted"] not in {allele["ref"],allele["alt"]}:
            print(f"  SKIP: allele orientation unavailable for {snp}"); total_skipped+=1; continue

        flip_genotype=allele["counted"]==allele["ref"]
        print(f"  SNP alleles: REF={allele['ref']} ALT={allele['alt']} PLINK_counted={allele['counted']} | {'flipping to ALT dosage' if flip_genotype else 'already ALT dosage'}")

        loc=find_location(trait_locations,event)
        if loc is None:print("  SKIP: coordinates unavailable"); total_skipped+=1; continue

        chrom=loc["chr"]; trait_start,trait_end=min(loc["start"],loc["end"]),max(loc["start"],loc["end"]); strand=loc.get("strand")
        flank=min(int(max(FLANK_MIN,math.ceil(max(1,trait_end-trait_start)*FLANK_FRAC))),FLANK_MAX)
        plot_start,plot_end=max(0,trait_start-flank),trait_end+flank
        n_bins=int(min(MAX_TRACK_POINTS,max(50,math.ceil((plot_end-plot_start)/BIGWIG_BIN_BP))))
        actual_bin_bp=(plot_end-plot_start)/n_bins
        x=np.linspace(plot_start,plot_end,n_bins,endpoint=False)+actual_bin_bp/2
        reference_exons=load_reference_exons(GTF,chrom,plot_start,plot_end)

        label=symbol if symbol else event; basename=safe_name(f"{tissue}__{TYPE}__{label}__{snp}")
        boxplot_out=OUTDIR/f"{basename}.boxplot.pdf"
        density_pdf,density_png,density_svg=OUTDIR/f"{basename}.density.pdf",OUTDIR/f"{basename}.density.png",OUTDIR/f"{basename}.density.svg"
        subject_tsv=OUTDIR/f"{basename}.subjects.tsv"

        if extract_boxplot(colored_boxplots,boxplot_pairs,event,snp,boxplot_out):total_boxplots+=1
        else:print("  WARNING: boxplot extraction failed; continuing density plot.")

        records=[]
        for qid in [q for q in phenotype if q in genotype]:
            g=genotype_group(genotype[qid],flip=flip_genotype)
            if g is None:continue
            try:pheno=float(phenotype[qid])
            except:continue
            if not np.isfinite(pheno):continue

            meta=metadata_for_qtl_id(qid,tissue)
            if meta is None:continue
            sample_id=str(meta.get(sample_col,"")).strip(); subject_id=str(meta.get(subject_col,"")).strip()
            if not sample_id:continue

            sample_dir=resolve_sample_dir(tissue,sample_id)
            if sample_dir is None:
                print(f"    SAMPLE DIR MISSING: {qid} | {sample_id} | tissue={tissue}")
                continue

            records.append({"qtl_id":qid,"subject_id":subject_id,"sample_id":sample_id,"sample_fs_id":sample_dir.name,
                            "genotype":g,"genotype_raw":genotype[qid],"phenotype":pheno,"sample_dir":sample_dir,
                            "primary_reads":None,"plus_track":None,"minus_track":None,"plus_source":None,"minus_source":None,
                            "circ_reads":None,"circ_rpm":None,"linear_reads":None,"circ_percent_file":None,
                            "circ_quant_source":None,"circ_match_type":None,"circ_matched_start":None,"circ_matched_end":None})

        print(f"  Candidate subjects before track QC: {len(records)}")

        for r in records:
            r["primary_reads"]=primary_read_count(r["sample_dir"])
            r["plus_track"],r["minus_track"],r["plus_source"],r["minus_source"]=bigwig_strands(r["sample_dir"],chrom,plot_start,plot_end,n_bins)

            if TYPE=="cQTL":
                qfile=circ_quant_file(r["sample_dir"]); r["circ_quant_source"]=str(qfile)
                q=read_circ_quant(qfile,chrom,trait_start,trait_end,strand)
                if q is not None:
                    r["circ_reads"],r["linear_reads"]=q["circ_reads"],q["linear_reads"]
                    r["circ_percent_file"],r["circ_match_type"]=q["circ_percent"],q["match_type"]
                    r["circ_matched_start"],r["circ_matched_end"]=q["matched_start"],q["matched_end"]
                    r["circ_rpm"]=rpm_value(q["circ_reads"],r["primary_reads"])

        complete=[r for r in records if r["plus_track"] is not None and r["minus_track"] is not None and len(r["plus_track"])==n_bins and len(r["minus_track"])==n_bins
                  and np.all(np.isfinite(r["plus_track"])) and np.all(np.isfinite(r["minus_track"]))]
        counts={g:sum(r["genotype"]==g for r in complete) for g in [0,1,2]}

        print(f"  Complete subjects used in density plots: {len(complete)}")
        print(f"  Ref/Ref={counts[0]} | Het={counts[1]} | Hom Alt={counts[2]}")

        if TYPE in {"sQTL","cQTL"}:
            n_primary=sum(r["primary_reads"] is not None for r in complete)
            print(f"  Samples with primary-read denominator: {n_primary}/{len(complete)}")

        if not complete:
            print("  SKIP DENSITY: no subjects with usable strand-specific tracks")
            total_skipped+=1
            continue

        groups={g:[r for r in complete if r["genotype"]==g] for g in [0,1,2]}
        mean_plus={g:np.nanmean(np.vstack([r["plus_track"] for r in groups[g]]),axis=0) if groups[g] else np.zeros(n_bins) for g in [0,1,2]}
        mean_minus={g:np.nanmean(np.vstack([r["minus_track"] for r in groups[g]]),axis=0) if groups[g] else np.zeros(n_bins) for g in [0,1,2]}

        junctions,junction_denoms={},{}
        if TYPE in {"sQTL","cQTL"}:
            for g in [0,1,2]:
                junctions[g],junction_denoms[g]=mean_junction_rpm(groups[g],chrom,plot_start,plot_end)
                if groups[g]:print(f"  {GENOTYPE_LABELS[g]} normalized LeafCutter samples: {junction_denoms[g]}/{len(groups[g])}; junctions plotted={len(junctions[g])}")

        circ_means={g:0. for g in [0,1,2]}; circ_denoms={g:0 for g in [0,1,2]}
        if TYPE=="cQTL":
            for g in [0,1,2]:
                valid=[r for r in groups[g] if r["primary_reads"]]
                circ_denoms[g]=len(valid)
                if valid:circ_means[g]=sum(0. if r["circ_rpm"] is None else float(r["circ_rpm"]) for r in valid)/len(valid)

            n_quant=sum(r["circ_reads"] is not None for r in complete)
            n_exact=sum(r["circ_match_type"]=="exact" for r in complete); n_tol=sum(r["circ_match_type"]=="tolerance" for r in complete)
            print(f"  Circ quant rows found: {n_quant}/{len(complete)}")
            print(f"  Circ coordinate matches: exact={n_exact}, tolerance-rescued={n_tol}")
            print("  Mean circ RPM: "+" | ".join(f"{GENOTYPE_LABELS[g]}={circ_means[g]:.4f}" for g in [0,1,2]))

        with open(subject_tsv,"w",newline="") as fh:
            writer=csv.writer(fh,delimiter="\t")
            writer.writerow(["tissue","qtl_type","event","gene_symbol","snp","ref","alt","plink_counted_allele","qtl_id","subject_id",
                             "metadata_sample_id","filesystem_sample_id","genotype_group","genotype_label","genotype_raw","phenotype",
                             "primary_reads","plus_track_source","minus_track_source","circ_reads","circ_rpm","linear_reads",
                             "circ_percent_file","circ_match_type","circ_matched_start","circ_matched_end","circ_quant_source"])
            for r in complete:
                writer.writerow([tissue,TYPE,event,symbol,snp,allele["ref"],allele["alt"],allele["counted"],r["qtl_id"],r["subject_id"],
                                 r["sample_id"],r["sample_fs_id"],r["genotype"],GENOTYPE_LABELS[r["genotype"]],r["genotype_raw"],r["phenotype"],
                                 "" if r["primary_reads"] is None else r["primary_reads"],r["plus_source"],r["minus_source"],
                                 "" if r["circ_reads"] is None else r["circ_reads"],"" if r["circ_rpm"] is None else r["circ_rpm"],
                                 "" if r["linear_reads"] is None else r["linear_reads"],"" if r["circ_percent_file"] is None else r["circ_percent_file"],
                                 "" if r["circ_match_type"] is None else r["circ_match_type"],"" if r["circ_matched_start"] is None else r["circ_matched_start"],
                                 "" if r["circ_matched_end"] is None else r["circ_matched_end"],"" if r["circ_quant_source"] is None else r["circ_quant_source"]])

        if not OVERWRITE and density_pdf.exists() and density_png.exists() and density_svg.exists():
            print("  DENSITY EXISTS; skipping")
            continue

        coverage_max=max(max(float(np.nanmax(mean_plus[g])),float(np.nanmax(mean_minus[g]))) for g in [0,1,2] if groups[g])
        coverage_max=1. if not np.isfinite(coverage_max) or coverage_max<=0 else coverage_max*1.05
        has_linear_arcs=TYPE in {"sQTL","cQTL"} and any(junctions.get(g) for g in [0,1,2])
        has_circ_arcs=TYPE=="cQTL" and any(circ_means[g]>0 for g in [0,1,2])
        axis_extent=coverage_max*(1.95 if has_circ_arcs else 1.58 if has_linear_arcs else 1.08)
        bar_width=actual_bin_bp*.92

        fig=plt.figure(figsize=(16,10))
        gs=fig.add_gridspec(4,1,height_ratios=[1.,1.,1.,.72],hspace=.10)
        axes=[]

        for g in [0,1,2]:
            ax=fig.add_subplot(gs[g,0]); axes.append(ax)

            if groups[g]:
                ax.bar(x,mean_plus[g],width=bar_width,color=TRACK_COLORS[g],alpha=.82,align="center",linewidth=0)
                ax.bar(x,-mean_minus[g],width=bar_width,color=TRACK_COLORS[g],alpha=.55,align="center",linewidth=0)
                if TYPE in {"sQTL","cQTL"}:draw_linear_arcs(ax,junctions[g],x,mean_plus[g],mean_minus[g],coverage_max,TRACK_COLORS[g])
                if TYPE=="cQTL":draw_circ_arc(ax,trait_start,trait_end,circ_means[g],coverage_max,TRACK_COLORS[g],strand)
            else:ax.text(.5,.5,"No subjects",transform=ax.transAxes,ha="center",va="center",fontsize=12)

            ax.axhline(0,color="#555555",linewidth=.8)
            ax.axvspan(trait_start,trait_end,alpha=.10,color="#999999")
            if plot_start<=snp_pos<=plot_end:ax.axvline(snp_pos,linestyle="--",linewidth=1.3,alpha=.8,color="black")
            ax.set_xlim(plot_start,plot_end); ax.set_ylim(-axis_extent,axis_extent)
            ax.yaxis.set_major_formatter(FuncFormatter(density_formatter))
            ax.set_title(f"{GENOTYPE_LABELS[g]} (N={counts[g]})",loc="left",fontsize=13,color=TRACK_COLORS[g],y=.84,pad=0)
            ax.tick_params(axis="x",labelbottom=False); ax.tick_params(axis="y",labelsize=11)
            ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
            ax.text(.995,.94,"Forward (+)",transform=ax.transAxes,ha="right",va="top",fontsize=10)
            ax.text(.995,.06,"Reverse (-)",transform=ax.transAxes,ha="right",va="bottom",fontsize=10)

        ax_ref=fig.add_subplot(gs[3,0],sharex=axes[0])
        draw_reference_track(ax_ref,reference_exons,plot_start,plot_end)
        ax_ref.axvspan(trait_start,trait_end,alpha=.10,color="#999999")
        if plot_start<=snp_pos<=plot_end:ax_ref.axvline(snp_pos,linestyle="--",linewidth=1.,alpha=.7,color="black")
        ax_ref.set_xlabel(f"{chrom} genomic position",fontsize=14,labelpad=10)

        title_event=f"{symbol} ({event})" if symbol and symbol!=event else event
        title=f"{tissue_display(tissue)} | {TYPE} | {title_event} | {snp} | REF={allele['ref']} ALT={allele['alt']}"
        if pp_h4 and not missing(pp_h4):
            try:title+=f" | PP.H4={float(pp_h4):.4f}"
            except:title+=f" | PP.H4={pp_h4}"

        fig.suptitle(title,fontsize=17,y=.995)
        fig.text(.008,.57,"Mean RPM-normalized RNA-seq read density",rotation=90,va="center",fontsize=15)

        base_footer=(f"Trait: {chrom}:{trait_start:,}-{trait_end:,} | flank: {flank:,} bp each side | shaded region = trait locus | "
                     f"coverage bin width ≈ {round(actual_bin_bp)} bp\n"
                     "Coverage uses the same RPM normalization as the source BigWigs (1e6 / primary alignments); + strand above 0, - strand below 0")

        if TYPE=="cQTL":
            footer=base_footer+"\nLinear arc labels = mean LeafCutter junction RPM per subject | circular arc label = mean circular-junction RPM per subject"
        elif TYPE=="sQTL":
            footer=base_footer+"\nArc labels = mean LeafCutter junction RPM per subject"
        else:
            footer=base_footer

        fig.text(.5,.004,footer,ha="center",va="bottom",fontsize=9.5)

        fig.savefig(density_pdf,bbox_inches="tight")
        fig.savefig(density_png,dpi=DPI,bbox_inches="tight")
        fig.savefig(density_svg,bbox_inches="tight")
        plt.close(fig)

        print(f"  BOXPLOT PDF : {boxplot_out if boxplot_out.exists() else 'NOT EXTRACTED'}")
        print(f"  DENSITY PDF : {density_pdf}")
        print(f"  DENSITY PNG : {density_png}")
        print(f"  DENSITY SVG : {density_svg}")
        print(f"  SUBJECT TSV : {subject_tsv}")
        print(f"  Coverage bin width: ≈{round(actual_bin_bp)} bp (actual {actual_bin_bp:.2f} bp)")
        total_plotted+=1

print(f"\n{'='*60}")
print("DONE")
print(f"QTL type           : {TYPE}")
print(f"Density plots      : {total_plotted}")
print(f"Boxplots extracted : {total_boxplots}")
print(f"Skipped hits       : {total_skipped}")
print(f"Output             : {OUTDIR}")
print("="*60)
PY
