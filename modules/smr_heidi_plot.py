#!/usr/bin/env python3

import os,re,csv,sys,math,warnings,subprocess
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Rectangle,PathPatch,Patch
from matplotlib.path import Path as MplPath
from scipy.stats import t as student_t
import pyBigWig
warnings.filterwarnings("ignore")

GENOTYPE_LABELS={0:"Ref/Ref",1:"Het",2:"Hom Alt"}
TRACK_COLORS={0:"#228B22",1:"#4682B4",2:"#C75B7A"}
CANONICAL_TISSUES=["Cerebellum","Frontal_Cortex","Cervical_Spinal_Cord","Lumbar_Spinal_Cord","Motor_Cortex"]

def norm(s):return re.sub(r"[^a-z0-9]","",str(s).lower())
def chr_key(s):return re.sub(r"^chr","",str(s),flags=re.I).upper()
def safe_name(s):return re.sub(r"[^A-Za-z0-9._+-]+","_",str(s))[:180]
def strip_gene_version(x):return re.sub(r"\.\d+$","",str(x))
def tissue_display(x):return str(x).replace("_"," ")
def missing(v):return v is None or str(v).strip().lower() in {"","na","nan","none",".","null"}
def truthy(v):return str(v).strip().lower() in {"true","t","1","yes","y","pass","passed"}
def detect_delimiter(line):return "\t" if "\t" in line else None
def split_line(line,d):return line.rstrip("\r\n").split("\t") if d=="\t" else line.split()
def as_float(x):
    try:v=float(x);return v if np.isfinite(v) else None
    except:return None
def format_p(x):return "NA" if as_float(x) is None else f"{float(x):.3e}"
def format_max_sig(v,sig=5):
    if v==0:return "0"
    d=max(0,sig-int(math.floor(math.log10(abs(v))))-1);return f"{v:.{d}f}".rstrip("0").rstrip(".")
def genomic_formatter(x,pos=None):return f"{format_max_sig(x/1e6)} Mb" if abs(x)>=1e6 else f"{format_max_sig(x/1e3)} kb" if abs(x)>=1e3 else format_max_sig(x)
def density_formatter(y,pos=None):return "0" if abs(y)<1e-12 else format_max_sig(abs(y),4)
def find_header(headers,aliases,required=True):
    lookup={norm(h):h for h in headers}
    for a in aliases:
        if norm(a) in lookup:return lookup[norm(a)]
    for h in headers:
        if any(norm(a) in norm(h) for a in aliases if norm(a)):return h
    if required:raise KeyError(f"Could not find {aliases}. Available: {headers}")
    return None

TISSUE_REMAP={
"Motor Cortex Lateral":"Motor_Cortex","Motor Cortex Medial":"Motor_Cortex","Lateral Motor Cortex":"Motor_Cortex","Medial Motor Cortex":"Motor_Cortex",
"Primary Motor Cortex L":"Motor_Cortex","Primary Motor Cortex M":"Motor_Cortex","Cortex_Motor_Unspecified":"Motor_Cortex","Cortex_Motor_BA4":"Motor_Cortex",
"BA4 Motor Cortex":"Motor_Cortex","Lateral_motor_cortex":"Motor_Cortex","Motor_Cortex_Lateral":"Motor_Cortex","Motor_Cortex_Medial":"Motor_Cortex",
"Lateral_Motor_Cortex":"Motor_Cortex","Medial_Motor_Cortex":"Motor_Cortex","Primary_Motor_Cortex_L":"Motor_Cortex","Primary_Motor_Cortex_M":"Motor_Cortex",
"Frontal Cortex":"Frontal_Cortex","Cerebellum":"Cerebellum","Spinal_Cord_Cervical":"Cervical_Spinal_Cord","Cervical Spinal Cord":"Cervical_Spinal_Cord",
"Cervical_spinal_cord":"Cervical_Spinal_Cord","Spinal_cord_Cervical":"Cervical_Spinal_Cord","Lumbar Spinal Cord":"Lumbar_Spinal_Cord",
"Spinal_Cord_Lumbosacral":"Lumbar_Spinal_Cord","Lumbosacral_Spinal_Cord":"Lumbar_Spinal_Cord","Lumbar_spinal_cord":"Lumbar_Spinal_Cord"}
TISSUE_REMAP_NORM={norm(k):v for k,v in TISSUE_REMAP.items()}
for t in CANONICAL_TISSUES:TISSUE_REMAP_NORM[norm(t)]=t
def canonical_tissue(x):x=str(x).strip();return TISSUE_REMAP_NORM.get(norm(x),x)
def tissue_equal(a,b):return canonical_tissue(a)==canonical_tissue(b)

def extract_matrix_rows(path,wanted):
    wanted=set(map(str,wanted));wanted_nv={strip_gene_version(x) for x in wanted};out={}
    with open(path) as f:
        first=f.readline();d=detect_delimiter(first);samples=split_line(first,d)[1:]
        for line in f:
            if not line.strip():continue
            x=split_line(line,d)
            if len(x)<2:continue
            rid=x[0].strip()
            if rid not in wanted and strip_gene_version(rid) not in wanted_nv:continue
            vals=x[1:]+[""]*max(0,len(samples)-len(x[1:]));out[rid]=dict(zip(samples,vals[:len(samples)]))
    return out
def find_matrix_row(rows,event):
    if event in rows:return rows[event]
    e=strip_gene_version(event);return next((v for k,v in rows.items() if strip_gene_version(k)==e),None)

def load_covariate_status(path):
    if not path.exists():return {}
    with open(path) as f:
        first=f.readline();d=detect_delimiter(first);samples=split_line(first,d)[1:]
        for line in f:
            x=split_line(line,d)
            if x and x[0].strip()=="is_als":
                vals=x[1:]+[""]*max(0,len(samples)-len(x[1:]))
                out={}
                for s,v in zip(samples,vals[:len(samples)]):
                    try:out[s]=1 if float(v)>.5 else 0
                    except:out[s]=0
                return out
    return {}

def load_location_table(path,kind,TYPE):
    with open(path) as f:
        first=f.readline();d=detect_delimiter(first);fields=split_line(first,d)
        is_data=len(fields)>=3 and re.match(r"^(chr)?[0-9XYM]+$",fields[1],re.I) and re.match(r"^\d+$",fields[2])
        if is_data:
            headers=["id","chr","start"]+(["end"] if len(fields)>=4 else [])+(["strand"] if len(fields)>=5 else [])
            while len(headers)<len(fields):headers.append(f"extra{len(headers)}")
            lines=[fields]
        else:headers,lines=fields,[]
        lines += [split_line(x,d) for x in f if x.strip()]
    rows=[]
    for x in lines:x+=[""]*max(0,len(headers)-len(x));rows.append(dict(zip(headers,x)))
    if kind=="snp":
        i=find_header(headers,["snpid","snp_id","rsid","id"]);c=find_header(headers,["chr","chrom","chromosome"]);p=find_header(headers,["pos","position","start"])
        return {str(r[i]):(str(r[c]),int(float(r[p]))) for r in rows if not any(missing(r.get(z)) for z in [i,c,p])}
    aliases={"eQTL":["geneid","gene_id","gene","id"],"sQTL":["junction_id","junctionid","splicing_id","event_id","phenotype_id","id"],"cQTL":["circ_id","circid","circrna_id","event_id","phenotype_id","id"]}
    i=find_header(headers,aliases[TYPE]);c=find_header(headers,["chr","chrom","chromosome"]);s=find_header(headers,["start","left","begin","donor"])
    e=find_header(headers,["end","right","stop","acceptor"]);st=find_header(headers,["strand"],False);out={}
    for r in rows:
        try:out[str(r[i])]={"chr":str(r[c]),"start":int(float(r[s])),"end":int(float(r[e])),"strand":str(r[st]) if st and not missing(r.get(st)) else None}
        except:pass
    return out
def parse_event_coordinates(event):
    for p in [r"(chr[^:]+):(\d+)-(\d+):([+-])",r"(chr[^:]+):(\d+):(\d+):([+-])",r"(chr[^:]+):(\d+)-(\d+)",r"(chr[^:]+):(\d+):(\d+)"]:
        m=re.search(p,str(event))
        if m:return {"chr":m.group(1),"start":int(m.group(2)),"end":int(m.group(3)),"strand":m.group(4) if len(m.groups())>=4 else None}
    return None
def find_location(loc,event):
    if event in loc:return loc[event]
    e=strip_gene_version(event);return next((v for k,v in loc.items() if strip_gene_version(k)==e),parse_event_coordinates(event))

def lm_pvalue(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);ok=np.isfinite(x)&np.isfinite(y);x,y=x[ok],y[ok]
    n=len(x)
    if n<3 or np.var(x)<=0:return "NA"
    X=np.column_stack([np.ones(n),x])
    try:
        beta=np.linalg.lstsq(X,y,rcond=None)[0];res=y-X@beta;df=n-2
        s2=np.sum(res**2)/df;cov=s2*np.linalg.inv(X.T@X);se=math.sqrt(cov[1,1])
        if se<=0:return "NA"
        p=2*student_t.sf(abs(beta[1]/se),df);return f"{p:.3e}"
    except:return "NA"

def make_matched_boxplot(complete,out,event,snp,symbol,allele,TYPE):
    if len(complete)<5:return False
    g=np.asarray([r["genotype"] for r in complete],int);y=np.asarray([r["phenotype"] for r in complete],float)
    als=np.asarray([int(r.get("is_als",0)) for r in complete],int)
    ylabel={"eQTL":"Expression (Z-score)","sQTL":"Splicing Level (PSI)","cQTL":"circRNA Percentage (%)"}[TYPE]
    rng=np.random.default_rng(529);fig,axs=plt.subplots(1,3,figsize=(12,5));models=[]
    models.append(("Additive Model",g,[0,1,2],["Ref/Ref","Het","Hom Alt"],"lightgreen"))
    dom=(g>0).astype(int);models.append(("Dominant Model",dom,[0,1],["Ref/Ref","Any Alt"],"lightblue"))
    rec=(g==2).astype(int);models.append(("Recessive Model",rec,[0,1],["Non-Hom Alt","Hom Alt"],"lightpink"))
    for ax,(title,xv,levels,names,color) in zip(axs,models):
        vals=[y[xv==z] for z in levels]
        ax.boxplot(vals,positions=np.arange(1,len(levels)+1),patch_artist=True,showfliers=False,
                   boxprops=dict(facecolor=color),medianprops=dict(color="black"))
        for j,z in enumerate(levels,1):
            idx=np.where(xv==z)[0]
            if not len(idx):continue
            xx=j+rng.uniform(-.15,.15,len(idx))
            for status,marker,size,c in [(1,"o",14,"red"),(0,"D",18,"#9932CC")]:
                ii=idx[als[idx]==status]
                if len(ii):
                    pos=[np.where(idx==k)[0][0] for k in ii]
                    ax.scatter(xx[pos],y[ii],s=size,c=c,marker=marker,linewidths=0)
        counts=[int(np.sum(xv==z)) for z in levels]
        ax.set_xticks(np.arange(1,len(levels)+1));ax.set_xticklabels([f"{n}\n(N={c})" for n,c in zip(names,counts)])
        ax.set_xlabel("Genotype");ax.set_ylabel(ylabel);ax.set_title(f"{title}\n(p = {lm_pvalue(xv,y)})")
        ax.plot([],[],marker="o",linestyle="None",color="red",markersize=4,label="ALS")
        ax.plot([],[],marker="D",linestyle="None",color="#9932CC",markersize=4,label="Control")
        ax.legend(loc="upper left",frameon=False,fontsize=8)
    if TYPE=="eQTL":head=f"Gene: {symbol} ({event})" if symbol and symbol!=event else f"Gene: {event}"
    elif TYPE=="sQTL":head=f"Junction: {event}"
    else:head=f"circRNA: {event}"
    fig.suptitle(f"{head} | SNP: {snp} | REF: {allele['ref']} | ALT: {allele['alt']}",fontsize=11,fontweight="bold")
    fig.tight_layout(rect=[0,0,.995,.91]);fig.savefig(out,bbox_inches="tight");plt.close(fig)
    print(f"  BOXPLOT REGENERATED: N={len(complete)} -> {out}");return True

processed_root_cache={};processed_sample_cache={};primary_read_cache={}
def processed_roots_for_tissue(ROOT,tissue):
    t=canonical_tissue(tissue)
    if t in processed_root_cache:return processed_root_cache[t]
    roots=[];preferred=ROOT/t/"RNAseq"/"Processed"
    if preferred.is_dir():roots.append(preferred)
    for p in ROOT.iterdir():
        if p.is_dir() and canonical_tissue(p.name)==t:
            q=p/"RNAseq"/"Processed"
            if q.is_dir() and q not in roots:roots.append(q)
    processed_root_cache[t]=roots;print(f"  RNAseq roots for {t}: "+(" | ".join(map(str,roots)) if roots else "NONE"));return roots
def resolve_sample_dir(ROOT,tissue,sid):
    key=(canonical_tissue(tissue),norm(sid))
    if key in processed_sample_cache:return processed_sample_cache[key]
    for root in processed_roots_for_tissue(ROOT,tissue):
        for x in [sid,sid.replace("-","_"),sid.replace("_","-")]:
            p=root/x
            if p.is_dir():processed_sample_cache[key]=p;return p
    for root in processed_roots_for_tissue(ROOT,tissue):
        rk=str(root)
        if rk not in processed_sample_cache:processed_sample_cache[rk]={norm(p.name):p for p in root.iterdir() if p.is_dir()}
        p=processed_sample_cache[rk].get(norm(sid))
        if p is not None:processed_sample_cache[key]=p;return p
    processed_sample_cache[key]=None;return None
def primary_read_count(sd):
    bam=sd/"STAR.Aligned.sortedByCoord.out.bam";k=str(bam)
    if k in primary_read_cache:return primary_read_cache[k]
    if not bam.exists():primary_read_cache[k]=None;return None
    try:n=int(subprocess.check_output(["samtools","view","-F","0x100","-c",str(bam)],text=True).strip());primary_read_cache[k]=n if n>0 else None
    except:primary_read_cache[k]=None
    return primary_read_cache[k]
def rpm_value(reads,n):return float(reads)*1e6/n if n and n>0 else None

def bw_chrom(bw,chrom):
    cs=bw.chroms();return next((c for c in [chrom,chrom[3:] if chrom.startswith("chr") else "chr"+chrom] if c in cs),None)
def bigwig_strands(sd,chrom,start,end,n_bins):
    paths=[sd/"STAR.Aligned.sortedByCoord.out.plus.normalized.bw",sd/"STAR.Aligned.sortedByCoord.out.minus.normalized.bw"]
    def read(path):
        if not path.exists():return None
        bw=pyBigWig.open(str(path))
        try:
            c=bw_chrom(bw,chrom)
            if c is None:return None
            return np.abs(np.asarray([0 if v is None else float(v) for v in bw.stats(c,int(start),int(end),nBins=n_bins,type="mean")]))
        finally:bw.close()
    return read(paths[0]),read(paths[1]),str(paths[0]),str(paths[1])

leafcutter_cache={}
def leafcutter_file(sd):return sd/"leafcutter"/"psi"/f"{sd.name}.leafcutter.PSI.tsv"
def load_leafcutter(path):
    k=str(path)
    if k in leafcutter_cache:return leafcutter_cache[k]
    if not path.exists():leafcutter_cache[k]=None;return None
    out={}
    try:
        with open(path) as f:
            for r in csv.DictReader(f,delimiter="\t"):
                try:key=(chr_key(r["chrom"]),r["strand"].strip(),int(float(r["start"])),int(float(r["end"])));out[key]=out.get(key,0)+float(r["reads"])
                except:continue
    except:out=None
    leafcutter_cache[k]=out;return out
def mean_junction_rpm(group,chrom,start,end,MIN_MEAN_JUNCTION_RPM,MAX_JUNCTIONS):
    sums=defaultdict(float);n=0
    for r in group:
        j=load_leafcutter(leafcutter_file(r["sample_dir"]));primary=r["primary_reads"]
        if j is None or not primary:continue
        n+=1
        for (c,st,s,e),reads in j.items():
            if c==chr_key(chrom) and s>=start and e<=end:sums[(st,s,e)]+=rpm_value(reads,primary)
    if not n:return {},0
    means={k:v/n for k,v in sums.items()};means={k:v for k,v in means.items() if v>=MIN_MEAN_JUNCTION_RPM}
    if MAX_JUNCTIONS>0 and len(means)>MAX_JUNCTIONS:
        keep=sorted(means,key=means.get,reverse=True)[:MAX_JUNCTIONS];means={k:means[k] for k in keep}
    return means,n

def circ_quant_file(sd):return sd/"circularRNA_known_circ_percentage.txt"
def read_circ_quant(path,chrom,start,end,strand,tol):
    if not path.exists():return None
    exact=None;nearby=[]
    with open(path) as f:
        for r in csv.DictReader(f,delimiter="\t"):
            try:
                s,e=int(float(r["start"])),int(float(r["end"]));st=r["strand"].strip()
                if chr_key(r["chrom"])!=chr_key(chrom) or (strand in {"+","-"} and st!=strand):continue
                q={"circ_reads":float(r["circ_reads"]),"linear_reads":float(r["linear_reads"]),"circ_percent":float(r["circ_percent"]),"matched_start":s,"matched_end":e,"match_type":None}
                if s==start and e==end:q["match_type"]="exact";exact=q;break
                if abs(s-start)<=tol and abs(e-end)<=tol:q["match_type"]="tolerance";nearby.append(q)
            except:continue
    if exact:return exact
    if not nearby:return None
    nearby.sort(key=lambda r:abs(r["matched_start"]-start)+abs(r["matched_end"]-end));return nearby[0]

def parse_gtf_attrs(s):
    out={}
    for x in s.strip().split(";"):
        p=x.strip().split(" ",1)
        if len(p)==2:out[p[0]]=p[1].strip().strip('"')
    return out
def load_reference_exons(gtf,chrom,start,end):
    out=[]
    with open(gtf) as f:
        for l in f:
            if l.startswith("#"):continue
            x=l.rstrip().split("\t")
            if len(x)<9 or x[2]!="exon" or chr_key(x[0])!=chr_key(chrom):continue
            s,e=int(x[3])-1,int(x[4])
            if e<start or s>end:continue
            a=parse_gtf_attrs(x[8]);out.append({"start":s,"end":e,"strand":x[6],"gene":a.get("gene_name",a.get("gene_id","Unknown"))})
    return out
def merge_intervals(xs):
    if not xs:return []
    xs=sorted(xs);out=[list(xs[0])]
    for s,e in xs[1:]:
        if s<=out[-1][1]:out[-1][1]=max(out[-1][1],e)
        else:out.append([s,e])
    return out
def draw_reference_track(ax,exons,start,end):
    genes=defaultdict(list);strands={}
    for e in exons:genes[e["gene"]].append((e["start"],e["end"]));strands[e["gene"]]=e["strand"]
    if not genes:ax.text(.5,.45,"No GENCODE v49 exons",transform=ax.transAxes,ha="center",fontsize=9)
    else:
        for y,(gene,ints) in enumerate(sorted(genes.items(),key=lambda z:min(x[0] for x in z[1]))):
            m=merge_intervals(ints);a=max(start,min(x[0] for x in m));b=min(end,max(x[1] for x in m));ax.plot([a,b],[y,y],color="black",lw=.8)
            for s,e in m:
                s,e=max(s,start),min(e,end)
                if e>s:ax.add_patch(Rectangle((s,y-.15),e-s,.3,facecolor="black",edgecolor="black",lw=.4))
            st=strands.get(gene,"");ax.text(start,y+.22,f"{gene} ({st})" if st in {"+","-"} else gene,fontsize=8,ha="left")
        ax.set_ylim(-.45,max(.55,len(genes)-.25))
    ax.set_xlim(start,end);ax.set_yticks([]);ax.set_ylabel("GENCODE v49",fontsize=12,labelpad=13)
    ax.xaxis.set_major_formatter(FuncFormatter(genomic_formatter));ax.locator_params(axis="x",nbins=5);ax.tick_params(axis="x",labelsize=11)
    for s in ["top","right","left"]:ax.spines[s].set_visible(False)

def draw_linear_arcs(ax,junctions,x,plus,minus,cmax,color):
    items=sorted(junctions.items(),key=lambda z:(z[0][0],(z[0][1]+z[0][2])/2,z[0][2]-z[0][1]));span=max(1,x[-1]-x[0])
    for i,((st,s,e),rpm) in enumerate(items):
        frac=min(1,max(0,(e-s)/span));mid=(s+e)/2;stagger=(i%3)*.10*cmax
        if st=="-":ly,ry=-float(np.interp(s,x,minus)),-float(np.interp(e,x,minus));base=min(ly,ry);apex=-cmax*(1.05+.42*math.sqrt(frac))-stagger
        else:ly,ry=float(np.interp(s,x,plus)),float(np.interp(e,x,plus));base=max(ly,ry);apex=cmax*(1.05+.42*math.sqrt(frac))+stagger
        path=MplPath([(s,ly),(mid,apex),(e,ry)],[MplPath.MOVETO,MplPath.CURVE3,MplPath.CURVE3]);lw=min(4,max(.7,math.log(rpm+1,2)))
        ax.add_patch(PathPatch(path,facecolor="none",edgecolor=color,lw=lw,alpha=.9))
        ax.text(mid,base+.58*(apex-base),f"{rpm:.3f}",ha="center",va="center",fontsize=7,bbox=dict(facecolor="white",edgecolor="none",alpha=.88,pad=.45))
def draw_circ_arc(ax,start,end,rpm,cmax,color,strand):
    if rpm<=0:return
    outer=(-1 if strand=="-" else 1)*cmax*1.72
    path=MplPath([(start,0),(start,outer),(end,outer),(end,0)],[MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4])
    lw=min(4,max(.7,math.log(rpm+1,2)));ax.add_patch(PathPatch(path,facecolor="none",edgecolor=color,lw=max(1.4,lw),linestyle="--",alpha=.95))
    ax.text((start+end)/2,outer*.63,f"{rpm:.3f}",ha="center",va="center",fontsize=7,bbox=dict(facecolor="white",edgecolor="none",alpha=.88,pad=.45))

def genotype_group(v,flip=False):
    try:x=float(v)
    except:return None
    if flip:x=2-x
    return int(np.clip(np.rint(x),0,2)) if np.isfinite(x) and -.1<=x<=2.1 else None

TYPE,ROOT,MR_DIR=os.environ["TYPE"],Path(os.environ["ROOT"]),Path(os.environ["MR_DIR"])
METADATA,OUTDIR,GTF,RAW_FILE=map(Path,[os.environ["METADATA"],os.environ["OUTDIR"],os.environ["GTF"],os.environ["RAW_FILE"]])
SMR_FDR_THRESHOLD=float(os.environ["SMR_FDR_THRESHOLD"]);FLANK_FRAC=float(os.environ["FLANK_FRAC"]);FLANK_MIN=int(os.environ["FLANK_MIN"]);FLANK_MAX=int(os.environ["FLANK_MAX"])
BIGWIG_BIN_BP=int(os.environ["BIGWIG_BIN_BP"]);MAX_TRACK_POINTS=int(os.environ["MAX_TRACK_POINTS"]);CIRC_COORD_TOL=int(os.environ["CIRC_COORD_TOL"])
MIN_MEAN_JUNCTION_RPM=float(os.environ["MIN_MEAN_JUNCTION_RPM"]);MAX_JUNCTIONS=int(os.environ["MAX_JUNCTIONS"]);MAX_HITS=int(os.environ["MAX_HITS"])
TISSUE_FILTER=os.environ["TISSUE_FILTER"].strip();OVERWRITE=int(os.environ["OVERWRITE"]);DPI=int(os.environ["DPI"])

def bh_fdr(p):
    p=np.asarray(p,float);n=len(p);o=np.argsort(p);q=p[o]*n/np.arange(1,n+1);q=np.minimum.accumulate(q[::-1])[::-1];q=np.minimum(q,1)
    out=np.empty(n);out[o]=q;return out
def load_smr_hits(tissue):
    f=ROOT/tissue/"MR"/f"{TYPE}_SMR_HEIDI"/"results"/f"{tissue}_{TYPE}_SMR_HEIDI.smr"
    if not f.exists():print(f"SMR FILE MISSING: {f}");return []
    with open(f,newline="") as h:r=csv.DictReader(h,delimiter="\t");rows=list(r);headers=r.fieldnames or []
    for c in ["probeID","topSNP","p_SMR"]:
        if c not in headers:sys.exit(f"ERROR: {f} missing {c}")
    valid=[]
    for x in rows:
        p=as_float(x.get("p_SMR"))
        if p is not None and 0<=p<=1:x["_p_SMR"]=p;valid.append(x)
    if not valid:return []
    for x,q in zip(valid,bh_fdr([x["_p_SMR"] for x in valid])):x["_SMR_FDR"]=float(q)
    hits=[x for x in valid if x["_SMR_FDR"]<SMR_FDR_THRESHOLD];hits.sort(key=lambda x:(x["_SMR_FDR"],x["_p_SMR"]))
    print(f"\n{'='*60}\n{tissue} | {TYPE} | valid SMR={len(valid):,} | FDR<{SMR_FDR_THRESHOLD:g}={len(hits):,} | HEIDI filter=NONE\n{'='*60}")
    return hits[:MAX_HITS] if MAX_HITS>0 else hits

hits_by_tissue={}
for t in CANONICAL_TISSUES:
    if TISSUE_FILTER and not tissue_equal(t,TISSUE_FILTER):continue
    h=load_smr_hits(t)
    if h:hits_by_tissue[t]=h
if not hits_by_tissue:sys.exit(f"No {TYPE} SMR hits with FDR < {SMR_FDR_THRESHOLD:g}")

raw_alleles={}
with open(RAW_FILE) as f:
    for col in f.readline().split():
        m=re.match(r"^(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)_([ACGT]+)$",col,re.I)
        if m:
            c,p,r,a,counted=m.groups();raw_alleles[(chr_key(c),int(p))]={"ref":r.upper(),"alt":a.upper(),"counted":counted.upper()}

with open(METADATA,newline="") as f:r=csv.DictReader(f);meta_headers,metadata=r.fieldnames or [],list(r)
subject_col=find_header(meta_headers,["externalsubjectid","subject_id","subjectid"]);sample_col=find_header(meta_headers,["externalsampleid","sample_id","sampleid"])
meta_tissue_col=find_header(meta_headers,["tissue"],False);rin_col=find_header(meta_headers,["rin","rna_integrity_number"],False)
def get_rin(r):
    try:return float(r.get(rin_col,"")) if rin_col else -1
    except:return -1
def metadata_for_qtl_id(qid,tissue):
    rows=metadata;target=canonical_tissue(tissue)
    if meta_tissue_col:
        z=[r for r in rows if canonical_tissue(r.get(meta_tissue_col,""))==target]
        if z:rows=z
    z=[r for r in rows if str(r.get(sample_col,"")).strip()==str(qid)]
    if z:return max(z,key=get_rin)
    z=[r for r in rows if str(r.get(subject_col,"")).strip()==str(qid)]
    return max(z,key=get_rin) if z else None

total_plotted=total_skipped=total_boxplots=0
for tissue,hits in hits_by_tissue.items():
    print(f"\n{'='*60}\nPROCESSING: {tissue} {TYPE} | SMR FDR hits={len(hits)}\n{'='*60}")
    qtl=ROOT/tissue/TYPE
    if TYPE=="eQTL":pheno_file,loc_file=qtl/f"expression_{tissue}.txt",qtl/"gene_location.txt"
    elif TYPE=="sQTL":pheno_file,loc_file=qtl/f"splicing_{tissue}.txt",qtl/"splicing_location.txt"
    else:pheno_file,loc_file=qtl/f"circ_{tissue}.txt",qtl/"circ_location.txt"
    snp_file,snp_loc_file=qtl/f"snp_{tissue}.txt",qtl/"snp_location.txt"
    cov_file=qtl/f"covariates_{tissue}_encoded.txt"
    miss=[p for p in [pheno_file,loc_file,snp_file,snp_loc_file] if not p.exists()]
    if miss:print("SKIP TISSUE: missing "+", ".join(map(str,miss)));total_skipped+=len(hits);continue
    processed_roots_for_tissue(ROOT,tissue)
    status=load_covariate_status(cov_file)
    events={h["probeID"].strip() for h in hits};snps={h["topSNP"].strip() for h in hits}
    pheno_rows=extract_matrix_rows(pheno_file,events);geno_rows=extract_matrix_rows(snp_file,snps)
    trait_locs=load_location_table(loc_file,"trait",TYPE);snp_locs=load_location_table(snp_loc_file,"snp",TYPE)

    for i,hit in enumerate(hits,1):
        event,snp=hit["probeID"].strip(),hit["topSNP"].strip();symbol=str(hit.get("Gene","")).strip()
        p_smr,fdr=hit["_p_SMR"],hit["_SMR_FDR"];p_heidi=as_float(hit.get("p_HEIDI"));nheidi=hit.get("nsnp_HEIDI","")
        print(f"\n[{i}/{len(hits)}] {event} | {symbol} | {snp} | SMR p={p_smr:.3e} FDR={fdr:.4g} HEIDI={format_p(p_heidi)}")
        phenotype=find_matrix_row(pheno_rows,event);genotype=geno_rows.get(snp)
        if phenotype is None or genotype is None or snp not in snp_locs:print("  SKIP: phenotype/genotype/SNP location missing");total_skipped+=1;continue
        snp_chr,snp_pos=snp_locs[snp];allele=raw_alleles.get((chr_key(snp_chr),int(snp_pos)))
        if allele is None or allele["counted"] not in {allele["ref"],allele["alt"]}:print("  SKIP: allele orientation unavailable");total_skipped+=1;continue
        flip=allele["counted"]==allele["ref"];loc=find_location(trait_locs,event)
        if loc is None:print("  SKIP: coordinates unavailable");total_skipped+=1;continue
        chrom=loc["chr"];tstart,tend=sorted([loc["start"],loc["end"]]);strand=loc.get("strand")
        flank=min(int(max(FLANK_MIN,math.ceil(max(1,tend-tstart)*FLANK_FRAC))),FLANK_MAX);pstart,pend=max(0,tstart-flank),tend+flank
        nb=int(min(MAX_TRACK_POINTS,max(50,math.ceil((pend-pstart)/BIGWIG_BIN_BP))));binbp=(pend-pstart)/nb
        x=np.linspace(pstart,pend,nb,endpoint=False)+binbp/2;exons=load_reference_exons(GTF,chrom,pstart,pend)
        label=symbol if symbol and symbol!=event else event;base=safe_name(f"{tissue}__{TYPE}__{label}__{snp}")
        boxout=OUTDIR/f"{base}.boxplot.pdf";pdf=OUTDIR/f"{base}.density.pdf";png=OUTDIR/f"{base}.density.png";svg=OUTDIR/f"{base}.density.svg";tsv=OUTDIR/f"{base}.subjects.tsv"

        records=[]
        for qid in phenotype.keys()&genotype.keys():
            g=genotype_group(genotype[qid],flip)
            try:ph=float(phenotype[qid])
            except:continue
            if g is None or not np.isfinite(ph):continue
            meta=metadata_for_qtl_id(qid,tissue)
            if not meta:continue
            sid=str(meta.get(sample_col,"")).strip();subj=str(meta.get(subject_col,"")).strip();sd=resolve_sample_dir(ROOT,tissue,sid) if sid else None
            if sd is None:
                print(f"    SAMPLE DIR MISSING: {qid} | {sid} | tissue={tissue}")
                continue
            records.append({"qtl_id":qid,"subject_id":subj,"sample_id":sid,"sample_fs_id":sd.name,"genotype":g,"genotype_raw":genotype[qid],"phenotype":ph,"is_als":status.get(qid,0),
                            "sample_dir":sd,"primary_reads":None,"plus_track":None,"minus_track":None,"plus_source":None,"minus_source":None,"circ_reads":None,"circ_rpm":None,"linear_reads":None,
                            "circ_percent_file":None,"circ_quant_source":None,"circ_match_type":None,"circ_matched_start":None,"circ_matched_end":None})
        print(f"  Candidate subjects before track QC: {len(records)}")
        for r in records:
            if TYPE in {"sQTL","cQTL"}:r["primary_reads"]=primary_read_count(r["sample_dir"])
            r["plus_track"],r["minus_track"],r["plus_source"],r["minus_source"]=bigwig_strands(r["sample_dir"],chrom,pstart,pend,nb)
            if TYPE=="cQTL":
                qf=circ_quant_file(r["sample_dir"]);r["circ_quant_source"]=str(qf);q=read_circ_quant(qf,chrom,tstart,tend,strand,CIRC_COORD_TOL)
                if q:
                    r["circ_reads"],r["linear_reads"],r["circ_percent_file"]=q["circ_reads"],q["linear_reads"],q["circ_percent"]
                    r["circ_match_type"],r["circ_matched_start"],r["circ_matched_end"]=q["match_type"],q["matched_start"],q["matched_end"]
                    r["circ_rpm"]=rpm_value(q["circ_reads"],r["primary_reads"])
        complete=[r for r in records if r["plus_track"] is not None and r["minus_track"] is not None and len(r["plus_track"])==nb and len(r["minus_track"])==nb and np.all(np.isfinite(r["plus_track"])) and np.all(np.isfinite(r["minus_track"]))]
        if not complete:print("  SKIP: no usable tracks");total_skipped+=1;continue
        groups={g:[r for r in complete if r["genotype"]==g] for g in [0,1,2]};counts={g:len(groups[g]) for g in groups}
        print(f"  Matched plotted subjects: {len(complete)} | Ref/Ref={counts[0]} | Het={counts[1]} | Hom Alt={counts[2]}")
        if make_matched_boxplot(complete,boxout,event,snp,symbol,allele,TYPE):total_boxplots+=1
        plus={g:np.mean(np.vstack([r["plus_track"] for r in groups[g]]),0) if groups[g] else np.zeros(nb) for g in groups}
        minus={g:np.mean(np.vstack([r["minus_track"] for r in groups[g]]),0) if groups[g] else np.zeros(nb) for g in groups}
        junctions={}
        if TYPE in {"sQTL","cQTL"}:
            for g in groups:junctions[g],_=mean_junction_rpm(groups[g],chrom,pstart,pend,MIN_MEAN_JUNCTION_RPM,MAX_JUNCTIONS)
        circ={g:0. for g in groups}
        if TYPE=="cQTL":
            for g in groups:
                valid=[r for r in groups[g] if r["primary_reads"]]
                if valid:circ[g]=sum(0 if r["circ_rpm"] is None else r["circ_rpm"] for r in valid)/len(valid)

        with open(tsv,"w",newline="") as f:
            w=csv.writer(f,delimiter="\t");w.writerow(["tissue","qtl_type","event","gene_symbol","snp","smr_p","smr_fdr","heidi_p","heidi_nsnp","ref","alt","qtl_id","subject_id","metadata_sample_id","filesystem_sample_id","genotype_group","genotype_label","phenotype","is_als","primary_reads","plus_track_source","minus_track_source","circ_reads","circ_rpm"])
            for r in complete:w.writerow([tissue,TYPE,event,symbol,snp,p_smr,fdr,"" if p_heidi is None else p_heidi,nheidi,allele["ref"],allele["alt"],r["qtl_id"],r["subject_id"],r["sample_id"],r["sample_fs_id"],r["genotype"],GENOTYPE_LABELS[r["genotype"]],r["phenotype"],r["is_als"],r["primary_reads"] or "",r["plus_source"],r["minus_source"],r["circ_reads"] or "",r["circ_rpm"] or ""])

        if not OVERWRITE and pdf.exists() and png.exists() and svg.exists():print("  DENSITY EXISTS; boxplot/TSV refreshed, density skipped");continue
        cmax=max(max(float(np.max(plus[g])),float(np.max(minus[g]))) for g in groups if groups[g]);cmax=1 if cmax<=0 else cmax*1.05
        linear=TYPE in {"sQTL","cQTL"} and any(junctions.get(g) for g in groups);circular=TYPE=="cQTL" and any(circ[g]>0 for g in groups)
        extent=cmax*(2.05 if circular else 1.80 if linear else 1.08);width=binbp*.92
        fig=plt.figure(figsize=(16,10));gs=fig.add_gridspec(4,1,height_ratios=[1,1,1,.72],hspace=.10);axes=[]
        for g in [0,1,2]:
            ax=fig.add_subplot(gs[g,0]);axes.append(ax)
            if groups[g]:
                ax.bar(x,plus[g],width=width,color=TRACK_COLORS[g],alpha=.82,linewidth=0);ax.bar(x,-minus[g],width=width,color=TRACK_COLORS[g],alpha=.55,linewidth=0)
                if TYPE in {"sQTL","cQTL"}:draw_linear_arcs(ax,junctions[g],x,plus[g],minus[g],cmax,TRACK_COLORS[g])
                if TYPE=="cQTL":draw_circ_arc(ax,tstart,tend,circ[g],cmax,TRACK_COLORS[g],strand)
            else:ax.text(.5,.5,"No subjects",transform=ax.transAxes,ha="center",va="center",fontsize=12)
            ax.axhline(0,color="#555555",lw=.8);ax.axvspan(tstart,tend,alpha=.10,color="#999999")
            if pstart<=snp_pos<=pend:ax.axvline(snp_pos,ls="--",lw=1.3,color="black",alpha=.8)
            ax.set_xlim(pstart,pend);ax.set_ylim(-extent,extent);ax.yaxis.set_major_formatter(FuncFormatter(density_formatter))
            ax.set_title(f"{GENOTYPE_LABELS[g]} (N={counts[g]})",loc="left",fontsize=13,color=TRACK_COLORS[g],y=.84,pad=0)
            ax.tick_params(axis="x",labelbottom=False);ax.tick_params(axis="y",labelsize=11);ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
            ax.text(.995,.94,"Forward (+)",transform=ax.transAxes,ha="right",va="top",fontsize=10);ax.text(.995,.06,"Reverse (-)",transform=ax.transAxes,ha="right",va="bottom",fontsize=10)
        ref=fig.add_subplot(gs[3,0],sharex=axes[0]);draw_reference_track(ref,exons,pstart,pend);ref.axvspan(tstart,tend,alpha=.10,color="#999999")
        if pstart<=snp_pos<=pend:ref.axvline(snp_pos,ls="--",lw=1,color="black",alpha=.7)
        ref.set_xlabel(f"{chrom} genomic position",fontsize=14,labelpad=10)
        title_event=f"{symbol} ({event})" if symbol and symbol!=event else event
        fig.suptitle(f"{tissue_display(tissue)} | {TYPE} | {title_event} | {snp} | REF={allele['ref']} ALT={allele['alt']} | SMR FDR={fdr:.4g}",fontsize=17,y=.995)
        fig.text(.008,.57,"Mean RPM-normalized RNA-seq read density",rotation=90,va="center",fontsize=15)
        footer=(f"SMR p={p_smr:.3e} | SMR FDR={fdr:.4g} | HEIDI p={format_p(p_heidi)} | HEIDI SNPs={nheidi if not missing(nheidi) else 'NA'} | HEIDI not used for filtering\n"
                f"Trait: {chrom}:{tstart:,}-{tend:,} | flank: {flank:,} bp each side | shaded region = trait locus | dashed vertical line = {snp} | coverage bin width ≈ {round(binbp)} bp\n"
                "Coverage uses the same RPM normalization as the source BigWigs (1e6 / primary alignments); + strand above 0, - strand below 0")
        if TYPE=="cQTL":footer+="\nLinear splice junctions = solid arcs | circular/back-splice junction = dashed arc | labels = mean junction RPM per subject"
        elif TYPE=="sQTL":footer+="\nArc labels = mean LeafCutter junction RPM per subject"
        fig.text(.5,.004,footer,ha="center",va="bottom",fontsize=9.5)
        fig.savefig(pdf,bbox_inches="tight");fig.savefig(png,dpi=DPI,bbox_inches="tight");fig.savefig(svg,bbox_inches="tight");plt.close(fig)
        print(f"  BOXPLOT: {boxout}\n  DENSITY: {pdf}\n  SUBJECT TSV: {tsv}");total_plotted+=1

print(f"\n{'='*60}\nDONE\nAnalysis           : SMR / HEIDI\nQTL type           : {TYPE}\nSMR FDR threshold  : {SMR_FDR_THRESHOLD}\nHEIDI filtering    : NONE\nDensity plots      : {total_plotted}\nBoxplots regenerated: {total_boxplots}\nSkipped hits       : {total_skipped}\nOutput             : {OUTDIR}\n{'='*60}")
