#!/usr/bin/env bash
#SBATCH --job-name=eQTL_MR
#SBATCH --output=/home/zw529/donglab/data/target_ALS/MR/eQTL_MR.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/MR/eQTL_MR.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

set -euo pipefail
module --force purge
command -v python3 >/dev/null || { echo "ERROR: python3 not found"; exit 1; }
python3 -c 'import pandas,numpy' 2>/dev/null || { echo "ERROR: pandas/numpy unavailable"; exit 1; }

BASE="$HOME/donglab/data/target_ALS"; ROOT="$HOME/donglab/data/GCST90027163/GWAS/eQTL_harmonization"; GLOBAL="$BASE/MR"
TISSUES=(Cervical_Spinal_Cord Lumbar_Spinal_Cord Motor_Cortex Frontal_Cortex Cerebellum); MIN_F=10
export BASE ROOT GLOBAL MIN_F TISSUE_LIST="${TISSUES[*]}"

python3 <<'PY'
import os, math
import numpy as np, pandas as pd

BASE,ROOT,GLOBAL=os.environ["BASE"],os.environ["ROOT"],os.environ["GLOBAL"]
TISSUES=os.environ["TISSUE_LIST"].split(); MIN_F=float(os.environ["MIN_F"])
os.makedirs(GLOBAL,exist_ok=True)

def bh(p):
    p=pd.to_numeric(p,errors="coerce"); out=pd.Series(np.nan,index=p.index,dtype=float); ok=p.notna()
    if not ok.any(): return out
    x=p[ok].clip(0,1); order=x.sort_values().index; vals=x.loc[order].values; m=len(vals); adj=np.empty(m); run=1.
    for i in range(m-1,-1,-1):
        run=min(run,vals[i]*m/(i+1)); adj[i]=min(run,1.)
    out.loc[order]=adj
    return out

def p_from_z(z): return math.erfc(abs(z)/math.sqrt(2)) if pd.notna(z) else np.nan
def exp_safe(x):
    if pd.isna(x): return np.nan
    if x>700: return np.inf
    if x<-700: return 0.
    return math.exp(x)

all_results=[]; summaries=[]

for tissue in TISSUES:
    infile=f"{ROOT}/{tissue}/{tissue}_SNP_gene_MR_ready.tsv.gz"
    out=f"{BASE}/{tissue}/MR"; os.makedirs(out,exist_ok=True)
    if not os.path.isfile(infile) or os.path.getsize(infile)==0: raise FileNotFoundError(f"Missing/empty: {infile}")

    print(f"\n{'='*70}\n{tissue}\n{'='*70}")
    d=pd.read_csv(infile,sep="\t",compression="gzip")

    req=["snpid","geneid","beta","se_qtl","p-value","gwas_beta_harmonized","gwas_se","gwas_p"]
    miss=set(req)-set(d.columns)
    if miss: raise ValueError(f"{tissue} input missing: {sorted(miss)}")

    for c in ["beta","se_qtl","p-value","gwas_beta_harmonized","gwas_se","gwas_p"]:
        d[c]=pd.to_numeric(d[c],errors="coerce")

    d=d.rename(columns={"beta":"beta_eqtl","p-value":"p_eqtl"})
    d["tissue"]=tissue; d["analysis_unit"]="SNP_gene_relationship"; d["method"]="SNP_Wald"
    d["F_statistic"]=(d["beta_eqtl"]/d["se_qtl"])**2
    d["strong_instrument"]=d["F_statistic"]>=MIN_F

    # SNP-specific MR: each SNP→gene relationship is estimated independently
    d["beta_MR"]=d["gwas_beta_harmonized"]/d["beta_eqtl"]
    d["se_MR"]=d["gwas_se"]/d["beta_eqtl"].abs()
    d["se_MR_delta"]=np.sqrt(
        d["gwas_se"]**2/d["beta_eqtl"]**2 +
        d["gwas_beta_harmonized"]**2*d["se_qtl"]**2/d["beta_eqtl"]**4
    )
    d["z_MR"]=d["beta_MR"]/d["se_MR"]; d["p_MR"]=d["z_MR"].apply(p_from_z)
    d["OR_MR"]=d["beta_MR"].apply(exp_safe)
    d["CI95_lower_beta"]=d["beta_MR"]-1.96*d["se_MR"]; d["CI95_upper_beta"]=d["beta_MR"]+1.96*d["se_MR"]
    d["OR_CI95_lower"]=d["CI95_lower_beta"].apply(exp_safe); d["OR_CI95_upper"]=d["CI95_upper_beta"].apply(exp_safe)

    d["MR_status"]="PASS"
    d.loc[~d["strong_instrument"],"MR_status"]="WEAK_INSTRUMENT"
    bad=d[["beta_eqtl","se_qtl","gwas_beta_harmonized","gwas_se","beta_MR","se_MR","p_MR"]].isna().any(axis=1)
    d.loc[bad,"MR_status"]="MISSING_VALUE"
    d.loc[(d["se_qtl"]<=0)|(d["gwas_se"]<=0)|(d["beta_eqtl"]==0),"MR_status"]="INVALID_VALUE"

    passed=d[d["MR_status"]=="PASS"].copy()
    passed["FDR_MR_tissue"]=bh(passed["p_MR"])
    passed=passed.sort_values(["FDR_MR_tissue","p_MR","snpid","geneid"])

    # Number of expression targets is annotation only — no gene aggregation
    ntargets=passed.groupby("snpid")["geneid"].nunique()
    passed["n_regulatory_targets"]=passed["snpid"].map(ntargets).astype(int)

    d.to_csv(f"{out}/{tissue}_SNP_MR_all.tsv.gz",sep="\t",index=False,compression="gzip")
    passed.to_csv(f"{out}/{tissue}_SNP_MR_F10.tsv.gz",sep="\t",index=False,compression="gzip")
    passed[passed["p_MR"]<.05].to_csv(f"{out}/{tissue}_SNP_MR_nominal_P0.05.tsv",sep="\t",index=False)
    passed[passed["FDR_MR_tissue"]<.05].to_csv(f"{out}/{tissue}_SNP_MR_FDR0.05.tsv",sep="\t",index=False)

    summaries.append({
        "tissue":tissue,
        "SNP_gene_MR_tests":len(passed),
        "unique_SNPs":passed["snpid"].nunique(),
        "target_genes":passed["geneid"].nunique(),
        "multi_target_SNPs":int((ntargets>=2).sum()),
        "nominal_P0.05":int((passed["p_MR"]<.05).sum()),
        "tissue_FDR0.05":int((passed["FDR_MR_tissue"]<.05).sum())
    })
    all_results.append(passed)

    print(f"SNP→gene MR relationships : {len(passed):,}")
    print(f"Unique SNPs               : {passed['snpid'].nunique():,}")
    print(f"Target genes              : {passed['geneid'].nunique():,}")
    print(f"SNPs regulating >=2 genes : {(ntargets>=2).sum():,}")
    print(f"P<0.05                    : {(passed['p_MR']<.05).sum():,}")
    print(f"Tissue FDR<0.05           : {(passed['FDR_MR_tissue']<.05).sum():,}")

# ==============================================================
# GLOBAL FDR ACROSS ALL SNP × GENE × TISSUE MR RELATIONSHIPS
# ==============================================================
allmr=pd.concat(all_results,ignore_index=True)
allmr["FDR_MR_global"]=bh(allmr["p_MR"])
allmr=allmr.sort_values(["FDR_MR_global","p_MR","tissue","snpid","geneid"])

for tissue in TISSUES:
    out=f"{BASE}/{tissue}/MR"
    x=allmr[allmr["tissue"]==tissue].copy()
    x.to_csv(f"{out}/{tissue}_SNP_MR_with_global_FDR.tsv.gz",sep="\t",index=False,compression="gzip")
    x[x["FDR_MR_global"]<.05].to_csv(f"{out}/{tissue}_SNP_MR_global_FDR0.05.tsv",sep="\t",index=False)

summary=pd.DataFrame(summaries)
summary.to_csv(f"{GLOBAL}/MR_summary_all_tissues.tsv",sep="\t",index=False)
allmr.to_csv(f"{GLOBAL}/all_tissues_SNP_MR.tsv.gz",sep="\t",index=False,compression="gzip")
allmr[allmr["FDR_MR_global"]<.05].to_csv(f"{GLOBAL}/all_tissues_SNP_MR_global_FDR0.05.tsv",sep="\t",index=False)

print(f"\n{'='*70}\nALL-TISSUE SNP-FIRST MR SUMMARY\n{'='*70}")
print(summary.to_string(index=False))
print(f"\nTotal SNP×gene MR tests : {len(allmr):,}")
print(f"Unique SNPs             : {allmr['snpid'].nunique():,}")
print(f"Global FDR<0.05         : {(allmr['FDR_MR_global']<.05).sum():,}")
PY

echo
echo "SNP-first MR analysis complete."
column -t -s $'\t' "$GLOBAL/MR_summary_all_tissues.tsv"
