#!/usr/bin/env bash
#SBATCH --job-name=eQTL_SMR_HEIDI
#SBATCH --output=/home/zw529/donglab/data/target_ALS/MR/eQTL_SMR_HEIDI.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/MR/eQTL_SMR_HEIDI.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

set -euo pipefail
module --force purge
module load PLINK/1.9b_7.11-x86_64

SMR="$HOME/donglab/pipelines/modules/smr/smr-1.4.2-linux-x86_64/smr"
BASE="$HOME/donglab/data/target_ALS"; GLOBAL="$BASE/MR/SMR_HEIDI"
BFILE="$BASE/QTL/plink/joint_all_chrs_filtered_bed"; RAW="$BASE/QTL/plink/joint_all_chrs_matrixEQTL.raw"; BIM="$BFILE.bim"
GWAS_DIR="$HOME/donglab/data/GCST90027163/GWAS"
GWAS="$GWAS_DIR/harmonised/34873335-GCST90027163-MONDO_0004976.h.tsv.gz"
GWAS_ORIG="$GWAS_DIR/GCST90027163_buildGRCh37.tsv.gz"; GWAS_MA="$GWAS_DIR/ALS_GRCh38_SMR.ma"
TISSUES=(Cervical_Spinal_Cord Lumbar_Spinal_Cord Motor_Cortex Frontal_Cortex Cerebellum)

PEQTL_SMR=5e-8; PEQTL_HEIDI=1.57e-3; HEIDI_PASS=0.05
LD_UPPER=0.90; LD_LOWER=0.05; HEIDI_MIN=3; HEIDI_MAX=20; CIS_WIND=2000

mkdir -p "$GLOBAL"
for f in "$SMR" "$BFILE.bed" "$BFILE.bim" "$BFILE.fam" "$RAW" "$GWAS" "$GWAS_ORIG"; do [[ -s "$f" ]] || { echo "ERROR: missing/empty: $f"; exit 1; }; done
command -v python3 >/dev/null || { echo "ERROR: python3 unavailable"; exit 1; }
python3 -c 'import pandas,numpy' 2>/dev/null || { echo "ERROR: pandas/numpy unavailable"; exit 1; }

export SMR BASE GLOBAL BFILE RAW BIM GWAS GWAS_ORIG GWAS_MA PEQTL_SMR PEQTL_HEIDI HEIDI_PASS LD_UPPER LD_LOWER HEIDI_MIN HEIDI_MAX CIS_WIND TISSUE_LIST="${TISSUES[*]}"

# ---------- Create persistent GRCh38 GWAS .ma only if missing/invalid ----------
if [[ -s "$GWAS_MA" ]] && python3 - "$GWAS_MA" <<'PY'
import sys, pandas as pd, numpy as np
f=sys.argv[1]
try:
    d=pd.read_csv(f,sep="\t")
    required=["SNP","A1","A2","freq","b","se","p","n"]
    assert list(d.columns)==required
    assert len(d)>1000 and d["SNP"].nunique()==len(d)
    assert d["A1"].astype(str).str.upper().isin(["A","C","G","T"]).all()
    assert d["A2"].astype(str).str.upper().isin(["A","C","G","T"]).all()
    b=pd.to_numeric(d["b"],errors="coerce"); se=pd.to_numeric(d["se"],errors="coerce"); p=pd.to_numeric(d["p"],errors="coerce")
    assert b.notna().all() and se.notna().all() and (se>0).all() and p.notna().all() and ((p>=0)&(p<=1)).all()
    print(f"VALID: {len(d):,} SNPs")
except Exception as e:
    print("INVALID:",e); sys.exit(1)
PY
then
  echo "Reusing validated GWAS SMR file: $GWAS_MA"
else
  echo "Creating GRCh38 GWAS SMR file: $GWAS_MA"
  rm -f "$GWAS_MA"
  python3 <<'PY'
import os, pandas as pd, numpy as np
g=os.environ["GWAS"]; go=os.environ["GWAS_ORIG"]; out=os.environ["GWAS_MA"]

want=["hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","p_value","hm_effect_allele_frequency"]
d=pd.read_csv(g,sep="\t",compression="gzip",usecols=lambda x:x in want,dtype=str)
need=["hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","p_value"]
miss=set(need)-set(d.columns)
if miss: raise ValueError(f"GWAS missing: {sorted(miss)}")

for c in ["hm_beta","standard_error","p_value"]: d[c]=pd.to_numeric(d[c],errors="coerce")
d["freq"]=pd.to_numeric(d["hm_effect_allele_frequency"],errors="coerce") if "hm_effect_allele_frequency" in d else np.nan

n=pd.read_csv(go,sep="\t",compression="gzip",usecols=["N_effective"])
N=int(round(pd.to_numeric(n["N_effective"],errors="coerce").median()))

d=d.rename(columns={"hm_rsid":"SNP","hm_effect_allele":"A1","hm_other_allele":"A2","hm_beta":"b","standard_error":"se","p_value":"p"})
d=d[d.SNP.notna()&d.A1.notna()&d.A2.notna()&d.b.notna()&d.se.notna()&(d.se>0)&d.p.notna()].copy()
d["A1"]=d.A1.str.upper(); d["A2"]=d.A2.str.upper(); d=d[d.A1.isin(list("ACGT"))&d.A2.isin(list("ACGT"))]
d["n"]=N; d=d.drop_duplicates("SNP")
d[["SNP","A1","A2","freq","b","se","p","n"]].to_csv(out,sep="\t",index=False,na_rep="NA")
print(f"Created {out}: {len(d):,} SNPs; N={N:,}")
PY
fi

# ---------- Per-tissue eQTL -> BESD -> SMR + HEIDI ----------
for TISSUE in "${TISSUES[@]}"; do
  echo; echo "======================================================================"; echo "$TISSUE"; echo "======================================================================"

  EQTL="$BASE/$TISSUE/eQTL"; RESULTS="$EQTL/results"; OUT="$BASE/$TISSUE/MR/SMR_HEIDI"
  INPUT="$OUT/input"; BESDDIR="$OUT/besd"; RESDIR="$OUT/results"
  FULL="$RESULTS/${TISSUE}_eQTL.full_annotated.txt"; TOP="$RESULTS/${TISSUE}_eQTL.top_for_boxplot.txt"
  GENELOC="$EQTL/gene_location.txt"; SNPLOC="$EQTL/snp_location.txt"; COV="$EQTL/covariates_${TISSUE}_encoded.txt"
  MATEQTL="$INPUT/${TISSUE}_MatrixEQTL_ENSG.txt"; KEEP="$INPUT/plink_keep.txt"; FREQ="$INPUT/${TISSUE}_freq"
  BESD="$BESDDIR/${TISSUE}_eQTL"; SMROUT="$RESDIR/${TISSUE}_SMR_HEIDI"

  rm -rf "$OUT"; mkdir -p "$INPUT" "$BESDDIR" "$RESDIR"
  for f in "$FULL" "$TOP" "$GENELOC" "$SNPLOC" "$COV"; do [[ -s "$f" ]] || { echo "ERROR: missing/empty: $f"; exit 1; }; done
  export TISSUE EQTL RESULTS OUT INPUT BESDDIR RESDIR FULL TOP GENELOC SNPLOC COV MATEQTL KEEP FREQ BESD SMROUT

  # Restore original ENSG identity using the same coordinate + symbol logic validated across all five tissues.
  python3 <<'PY'
import os, pandas as pd, numpy as np

T=os.environ["TISSUE"]; fullf=os.environ["FULL"]; topf=os.environ["TOP"]; glf=os.environ["GENELOC"]
locf=os.environ["SNPLOC"]; rawf=os.environ["RAW"]; bimf=os.environ["BIM"]; famf=os.environ["BFILE"]+".fam"
covf=os.environ["COV"]; outf=os.environ["MATEQTL"]; keepf=os.environ["KEEP"]

g=pd.read_csv(glf,sep="\t",header=None,names=["ensembl_id","chr","start","end"],dtype=str)
if len(g) and str(g.iloc[0,0]).lower() in {"geneid","gene_id","ensembl_id"}: g=g.iloc[1:]
g["chrkey"]=g.chr.str.replace("^chr","",regex=True).str.upper()
g["startkey"]=pd.to_numeric(g.start,errors="coerce").astype("Int64"); g["endkey"]=pd.to_numeric(g.end,errors="coerce").astype("Int64")
coord={k:sorted(z.ensembl_id.dropna().astype(str).unique()) for k,z in g.groupby(["chrkey","startkey","endkey"],dropna=False)}

top=pd.read_csv(topf,sep="\t",header=None,usecols=[0,1],names=["ensembl_id","symbol"],dtype=str).drop_duplicates()
if len(top) and str(top.iloc[0,0]).lower() in {"geneid","gene_id","ensembl_id"}: top=top.iloc[1:]
symmap={str(s):set(z.ensembl_id.astype(str)) for s,z in top.groupby("symbol")}

def resolve(symbol,chrom,start,end):
    try: key=(str(chrom).replace("chr","").upper(),int(start),int(end))
    except: return None
    ids=coord.get(key,[])
    if str(symbol) in ids: return str(symbol)
    sids=[x for x in ids if x in symmap.get(str(symbol),set())]
    if len(sids)==1: return sids[0]
    if len(ids)==1: return ids[0]
    return None

with open(rawf) as fh: hdr=fh.readline().split()
meta={"FID","IID","PAT","MAT","SEX","PHENOTYPE"}; rr=[]
for c in hdr:
    if c in meta or "_" not in c: continue
    vid,a=c.rsplit("_",1); rr.append((vid,a.upper()))
raw=pd.DataFrame(rr,columns=["rawid","effect_allele"])

bim=pd.read_csv(bimf,sep=r"\s+",header=None,names=["chr","rawid","cm","pos","a1","a2"],dtype=str)
bim["key"]=bim.chr.str.replace("^chr","",regex=True).str.upper()+":"+bim.pos
loc=pd.read_csv(locf,sep="\t",dtype=str); loc["key"]=loc.chr.str.replace("^chr","",regex=True).str.upper()+":"+loc.pos
vm=loc[["snpid","key"]].drop_duplicates().merge(bim[["rawid","key"]],on="key").merge(raw,on="rawid")
available=set(vm.snpid.astype(str))

fam=pd.read_csv(famf,sep=r"\s+",header=None,usecols=[0,1],names=["FID","IID"],dtype=str)
cov=pd.read_csv(covf,sep="\t",nrows=0); subjects=set(map(str,list(cov.columns)[1:]))
keep=fam[fam.IID.isin(subjects)]
if keep.empty: raise RuntimeError(f"{T}: no covariate subjects matched PLINK FAM")
keep.to_csv(keepf,sep="\t",header=False,index=False)

use=["geneid","snpid","beta","t-stat","p-value","FDR","gene_chr","gene_start","gene_end"]
first=True; total=written=unresolved=0; cache={}

for ch in pd.read_csv(fullf,sep="\t",usecols=use,dtype=str,chunksize=250000):
    total+=len(ch); ch=ch[ch.snpid.isin(available)].copy()
    if ch.empty: continue

    for r in ch[["geneid","gene_chr","gene_start","gene_end"]].drop_duplicates().itertuples(index=False):
        k=(str(r.geneid),str(r.gene_chr),str(r.gene_start),str(r.gene_end))
        if k not in cache: cache[k]=resolve(*k)

    ch["gene"]=[cache[(str(a),str(b),str(c),str(d))] for a,b,c,d in zip(ch.geneid,ch.gene_chr,ch.gene_start,ch.gene_end)]
    unresolved+=int(ch.gene.isna().sum()); ch=ch[ch.gene.notna()].copy()
    if ch.empty: continue

    for c in ["beta","t-stat","p-value","FDR"]: ch[c]=pd.to_numeric(ch[c],errors="coerce")
    ok=ch.beta.notna()&ch["t-stat"].notna()&(ch["t-stat"]!=0)&ch["p-value"].notna()&(ch["p-value"]>=0)&(ch["p-value"]<=1)&ch.FDR.notna()
    ch=ch[ok]

    out=ch[["snpid","gene","beta","t-stat","p-value","FDR"]].rename(columns={"snpid":"SNP"})
    out.to_csv(outf,sep="\t",index=False,mode="w" if first else "a",header=first)
    written+=len(out); first=False

if unresolved: print(f"WARNING: {T}: {unresolved:,} genotype-backed rows lacked recoverable ENSG identity and were excluded.")
print(f"{T}: source rows={total:,}; SMR MatrixEQTL rows={written:,}; eQTL/LD N={len(keep):,}")
PY

  # Build tissue frequency and BESD. MatrixEQTL format is natively supported by SMR.
  plink --bfile "$BFILE" --keep "$KEEP" --freq --threads "${SLURM_CPUS_PER_TASK:-1}" --out "$FREQ" >/dev/null
  N_EQTL=$(wc -l < "$KEEP")
  "$SMR" --eqtl-summary "$MATEQTL" --matrix-eqtl-format --make-besd --add-n "$N_EQTL" --out "$BESD"

  # Complete .esi/.epi metadata while preserving BESD ordering.
  python3 <<'PY'
import os, shutil, pandas as pd, numpy as np
T=os.environ["TISSUE"]; besd=os.environ["BESD"]; rawf=os.environ["RAW"]; bimf=os.environ["BIM"]; locf=os.environ["SNPLOC"]
glf=os.environ["GENELOC"]; topf=os.environ["TOP"]; frqf=os.environ["FREQ"]+".frq"

with open(rawf) as fh: hdr=fh.readline().split()
meta={"FID","IID","PAT","MAT","SEX","PHENOTYPE"}; rr=[]
for c in hdr:
    if c in meta or "_" not in c: continue
    vid,a=c.rsplit("_",1); rr.append((vid,a.upper()))
raw=pd.DataFrame(rr,columns=["rawid","effect"])

bim=pd.read_csv(bimf,sep=r"\s+",header=None,names=["chr","rawid","cm","pos","a1","a2"],dtype=str)
bim["a1"]=bim.a1.str.upper(); bim["a2"]=bim.a2.str.upper()
bim["key"]=bim.chr.str.replace("^chr","",regex=True).str.upper()+":"+bim.pos
x=raw.merge(bim,on="rawid",how="left"); x["other"]=np.where(x.effect==x.a1,x.a2,np.where(x.effect==x.a2,x.a1,None))

loc=pd.read_csv(locf,sep="\t",dtype=str); loc["key"]=loc.chr.str.replace("^chr","",regex=True).str.upper()+":"+loc.pos
m=loc[["snpid","chr","pos","key"]].drop_duplicates().merge(x[["rawid","key","effect","other"]],on="key",how="left")

frq=pd.read_csv(frqf,sep=r"\s+",dtype=str); frq["MAF"]=pd.to_numeric(frq.MAF,errors="coerce")
fm=m.merge(frq[["SNP","A1","A2","MAF"]],left_on="rawid",right_on="SNP",how="left")
fm["freq"]=np.where(fm.effect==fm.A1.str.upper(),fm.MAF,np.where(fm.effect==fm.A2.str.upper(),1-fm.MAF,np.nan))
snpmap=fm.set_index("snpid")

esi=besd+".esi"; shutil.copy2(esi,esi+".preupdate")
e=pd.read_csv(esi,sep=r"\s+",header=None,dtype=str); rows=[]; missing=[]
for s in e.iloc[:,1].astype(str):
    if s not in snpmap.index: missing.append(s); continue
    r=snpmap.loc[s]; r=r.iloc[0] if isinstance(r,pd.DataFrame) else r
    rows.append([str(r["chr"]).replace("chr",""),s,"0",str(r["pos"]),r["effect"],r["other"],r["freq"]])
if missing: raise RuntimeError(f"{T}: {len(missing)} BESD SNPs missing genotype metadata")
pd.DataFrame(rows).to_csv(esi,sep="\t",header=False,index=False,na_rep="NA")

g=pd.read_csv(glf,sep="\t",header=None,names=["gene","chr","start","end"],dtype=str)
if len(g) and str(g.iloc[0,0]).lower() in {"geneid","gene_id","ensembl_id"}: g=g.iloc[1:]
g["bp"]=((pd.to_numeric(g.start,errors="coerce")+pd.to_numeric(g.end,errors="coerce"))/2).round().astype("Int64")
top=pd.read_csv(topf,sep="\t",header=None,usecols=[0,1],names=["gene","symbol"],dtype=str).drop_duplicates()
sym=top.groupby("gene").symbol.agg(lambda z:";".join(sorted(set(map(str,z))))).to_dict()
gm=g.drop_duplicates("gene").set_index("gene")

epi=besd+".epi"; shutil.copy2(epi,epi+".preupdate")
ep=pd.read_csv(epi,sep=r"\s+",header=None,dtype=str); out=[]; missing=[]
for p in ep.iloc[:,1].astype(str):
    if p not in gm.index: missing.append(p); continue
    r=gm.loc[p]; r=r.iloc[0] if isinstance(r,pd.DataFrame) else r
    out.append([str(r["chr"]).replace("chr",""),p,"0",str(r["bp"]),sym.get(p,p),"NA"])
if missing: raise RuntimeError(f"{T}: {len(missing)} BESD probes missing gene_location metadata")
pd.DataFrame(out).to_csv(epi,sep="\t",header=False,index=False)

print(f"{T}: BESD metadata updated: {len(rows):,} SNPs; {len(out):,} probes")
PY

  gzip -f "$MATEQTL"

  # Standard cis-SMR + HEIDI. HEIDI pass (P>=0.05) is applied downstream.
  "$SMR" --bfile "$BFILE" --keep "$KEEP" --gwas-summary "$GWAS_MA" --beqtl-summary "$BESD" \
    --peqtl-smr "$PEQTL_SMR" --peqtl-heidi "$PEQTL_HEIDI" --ld-upper-limit "$LD_UPPER" --ld-lower-limit "$LD_LOWER" \
    --heidi-min-m "$HEIDI_MIN" --heidi-max-m "$HEIDI_MAX" --heidi-mtd 1 --cis-wind "$CIS_WIND" \
    --thread-num "${SLURM_CPUS_PER_TASK:-1}" --out "$SMROUT"

  [[ -s "$SMROUT.smr" ]] || { echo "ERROR: no SMR output for $TISSUE"; exit 1; }
done

# ---------- BH-FDR + HEIDI interpretation ----------
python3 <<'PY'
import os, numpy as np, pandas as pd
BASE=os.environ["BASE"]; GLOBAL=os.environ["GLOBAL"]; T=os.environ["TISSUE_LIST"].split(); HEIDI=float(os.environ["HEIDI_PASS"])

def bh(p):
    p=pd.to_numeric(p,errors="coerce"); out=pd.Series(np.nan,index=p.index,dtype=float); ok=p.notna()
    if not ok.any(): return out
    x=p[ok].clip(0,1); order=x.sort_values().index; v=x.loc[order].to_numpy(); m=len(v); a=np.empty(m); run=1.
    for i in range(m-1,-1,-1): run=min(run,v[i]*m/(i+1)); a[i]=min(run,1.)
    out.loc[order]=a; return out

allx=[]; summary=[]
for t in T:
    f=f"{BASE}/{t}/MR/SMR_HEIDI/results/{t}_SMR_HEIDI.smr"; d=pd.read_csv(f,sep=r"\s+",dtype=str)
    for c in ["p_SMR","p_HEIDI","b_SMR","se_SMR","p_GWAS","p_eQTL"]:
        if c in d: d[c]=pd.to_numeric(d[c],errors="coerce")
    d.insert(0,"tissue",t); d["FDR_SMR_tissue"]=bh(d["p_SMR"])
    d["HEIDI_tested"]=d["p_HEIDI"].notna(); d["HEIDI_pass"]=d["HEIDI_tested"]&(d["p_HEIDI"]>=HEIDI)
    d["SMR_tissue_FDR_pass"]=d["FDR_SMR_tissue"]<.05
    allx.append(d); summary.append({"tissue":t,"SMR_tests":len(d),"SMR_nominal_P0.05":int((d.p_SMR<.05).sum()),
      "SMR_tissue_FDR0.05":int((d.FDR_SMR_tissue<.05).sum()),"HEIDI_tested":int(d.HEIDI_tested.sum()),"HEIDI_P_ge_0.05":int(d.HEIDI_pass.sum())})

x=pd.concat(allx,ignore_index=True); x["FDR_SMR_global"]=bh(x["p_SMR"]); x["SMR_global_FDR_pass"]=x["FDR_SMR_global"]<.05
x["SMR_HEIDI_pass"]=x["SMR_global_FDR_pass"]&x["HEIDI_pass"]

for t in T:
    out=f"{BASE}/{t}/MR/SMR_HEIDI"; d=x[x.tissue==t].sort_values(["FDR_SMR_global","p_SMR"])
    d.to_csv(f"{out}/{t}_SMR_HEIDI_all.tsv.gz",sep="\t",index=False,compression="gzip")
    d[d.SMR_HEIDI_pass].to_csv(f"{out}/{t}_SMR_HEIDI_pass.tsv",sep="\t",index=False)

s=pd.DataFrame(summary)
for i,r in s.iterrows():
    d=x[x.tissue==r.tissue]
    s.loc[i,"SMR_global_FDR0.05"]=int((d.FDR_SMR_global<.05).sum())
    s.loc[i,"SMR_global_FDR0.05_and_HEIDI_pass"]=int(d.SMR_HEIDI_pass.sum())

s.to_csv(f"{GLOBAL}/SMR_HEIDI_summary_all_tissues.tsv",sep="\t",index=False)
x.to_csv(f"{GLOBAL}/all_tissues_SMR_HEIDI.tsv.gz",sep="\t",index=False,compression="gzip")
x[x.SMR_HEIDI_pass].to_csv(f"{GLOBAL}/all_tissues_SMR_HEIDI_pass.tsv",sep="\t",index=False)

print("\n======================================================================\nSMR + HEIDI SUMMARY\n======================================================================")
print(s.to_string(index=False))
print(f"\nTotal SMR tests: {len(x):,}")
print(f"Global SMR FDR<0.05: {(x.FDR_SMR_global<.05).sum():,}")
print(f"Global SMR FDR<0.05 + HEIDI P>={HEIDI}: {x.SMR_HEIDI_pass.sum():,}")
PY

echo; echo "SMR + HEIDI complete."; echo
column -t -s $'\t' "$GLOBAL/SMR_HEIDI_summary_all_tissues.tsv"
