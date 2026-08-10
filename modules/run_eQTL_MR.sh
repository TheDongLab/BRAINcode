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

BASE,ROOT,GLOBAL=os.environ["BASE"],os.environ["ROOT"],os.environ["GLOBAL"]; TISSUES=os.environ["TISSUE_LIST"].split(); MIN_F=float(os.environ["MIN_F"])
os.makedirs(GLOBAL,exist_ok=True)

def bh(p):
    p=pd.to_numeric(p,errors="coerce"); out=pd.Series(np.nan,index=p.index,dtype=float); ok=p.notna()
    if not ok.any(): return out
    x=p[ok].clip(0,1); order=x.sort_values().index; vals=x.loc[order].values; m=len(vals); adj=np.empty(m); run=1.
    for i in range(m-1,-1,-1): run=min(run,vals[i]*m/(i+1)); adj[i]=min(run,1.)
    out.loc[order]=adj; return out

def p_from_z(z): return math.erfc(abs(z)/math.sqrt(2)) if pd.notna(z) else np.nan
def exp_safe(x): return math.exp(x) if pd.notna(x) and x<700 else np.inf

lead_all=[]; adaptive_all=[]; summaries=[]

for tissue in TISSUES:
    lead_file=f"{ROOT}/{tissue}/{tissue}_eQTL_GWAS_MR_ready.tsv.gz"
    multi_file=f"{ROOT}/{tissue}/{tissue}_eQTL_GWAS_multiSNP_MR_ready.tsv.gz"
    out=f"{BASE}/{tissue}/MR"; os.makedirs(out,exist_ok=True)
    for f in [lead_file,multi_file]:
        if not os.path.isfile(f) or os.path.getsize(f)==0: raise FileNotFoundError(f"Missing/empty: {f}")

    print(f"\n{'='*70}\n{tissue}\n{'='*70}")

    # ==============================================================
    # A. ORIGINAL LEAD-SNP WALD ANALYSIS
    # ==============================================================
    d=pd.read_csv(lead_file,sep="\t",compression="gzip")
    req=["geneid","snpid","beta","se_qtl","p-value","gwas_beta_harmonized","gwas_se","gwas_p"]
    miss=set(req)-set(d.columns)
    if miss: raise ValueError(f"{tissue} lead file missing: {sorted(miss)}")

    for c in ["beta","se_qtl","p-value","gwas_beta_harmonized","gwas_se","gwas_p"]: d[c]=pd.to_numeric(d[c],errors="coerce")
    d=d.rename(columns={"beta":"beta_eqtl","p-value":"p_eqtl"}); d["tissue"]=tissue; d["method"]="Lead_SNP_Wald"
    d["F_statistic"]=(d["beta_eqtl"]/d["se_qtl"])**2; d["strong_instrument"]=d["F_statistic"]>=MIN_F
    d["beta_MR"]=d["gwas_beta_harmonized"]/d["beta_eqtl"]
    d["se_MR"]=d["gwas_se"]/d["beta_eqtl"].abs()
    d["se_MR_delta"]=np.sqrt(d["gwas_se"]**2/d["beta_eqtl"]**2 + d["gwas_beta_harmonized"]**2*d["se_qtl"]**2/d["beta_eqtl"]**4)
    d["z_MR"]=d["beta_MR"]/d["se_MR"]; d["p_MR"]=d["z_MR"].apply(p_from_z)
    d["OR_MR"]=d["beta_MR"].apply(exp_safe); d["CI95_lower_beta"]=d["beta_MR"]-1.96*d["se_MR"]; d["CI95_upper_beta"]=d["beta_MR"]+1.96*d["se_MR"]
    d["OR_CI95_lower"]=d["CI95_lower_beta"].apply(exp_safe); d["OR_CI95_upper"]=d["CI95_upper_beta"].apply(exp_safe)
    d["FDR_MR_tissue"]=bh(d["p_MR"]); d["MR_status"]=np.where(d["strong_instrument"],"PASS","WEAK_INSTRUMENT")
    d.loc[d[["beta_eqtl","se_qtl","gwas_beta_harmonized","gwas_se","beta_MR","se_MR","p_MR"]].isna().any(axis=1),"MR_status"]="MISSING_VALUE"

    lead_primary=d[d["MR_status"]=="PASS"].copy().sort_values(["FDR_MR_tissue","p_MR"])
    d.to_csv(f"{out}/{tissue}_Wald_MR_all.tsv.gz",sep="\t",index=False,compression="gzip")
    lead_primary.to_csv(f"{out}/{tissue}_Wald_MR_F{int(MIN_F)}.tsv.gz",sep="\t",index=False,compression="gzip")
    lead_primary[lead_primary["p_MR"]<.05].to_csv(f"{out}/{tissue}_Wald_MR_nominal_P0.05.tsv",sep="\t",index=False)
    lead_primary[lead_primary["FDR_MR_tissue"]<.05].to_csv(f"{out}/{tissue}_Wald_MR_FDR0.05.tsv",sep="\t",index=False)

    # ==============================================================
    # B. INDEPENDENT-SNP ANALYSIS: 1 SNP=WALD, >=2 SNPs=IVW
    # ==============================================================
    m=pd.read_csv(multi_file,sep="\t",compression="gzip")
    req=["geneid","snpid","beta","se_qtl","p-value","gwas_beta_harmonized","gwas_se","gwas_p"]
    miss=set(req)-set(m.columns)
    if miss: raise ValueError(f"{tissue} multi-SNP file missing: {sorted(miss)}")

    for c in ["beta","se_qtl","p-value","gwas_beta_harmonized","gwas_se","gwas_p"]: m[c]=pd.to_numeric(m[c],errors="coerce")
    m["F_statistic"]=(m["beta"]/m["se_qtl"])**2
    m=m[(m["F_statistic"]>=MIN_F)&m["beta"].notna()&m["gwas_beta_harmonized"].notna()&(m["gwas_se"]>0)].copy()

    rows=[]
    for gene,g in m.groupby("geneid",sort=False):
        g=g.drop_duplicates("snpid").copy(); n=len(g)
        bx=g["beta"].to_numpy(float); by=g["gwas_beta_harmonized"].to_numpy(float); sy=g["gwas_se"].to_numpy(float)

        if n==1:
            beta=by[0]/bx[0]; se=sy[0]/abs(bx[0]); method="Wald"
            Q=np.nan
        else:
            w=1/(sy**2); denom=np.sum(w*bx**2)
            if denom<=0: continue
            beta=np.sum(w*bx*by)/denom; se=math.sqrt(1/denom); method="IVW"
            Q=np.sum(w*(by-beta*bx)**2)

        z=beta/se; p=p_from_z(z)
        rows.append({"tissue":tissue,"geneid":gene,"method":method,"n_instruments":n,"snps":";".join(g["snpid"].astype(str)),
                     "beta_MR":beta,"se_MR":se,"z_MR":z,"p_MR":p,"OR_MR":exp_safe(beta),
                     "CI95_lower_beta":beta-1.96*se,"CI95_upper_beta":beta+1.96*se,
                     "OR_CI95_lower":exp_safe(beta-1.96*se),"OR_CI95_upper":exp_safe(beta+1.96*se),
                     "Q_heterogeneity":Q,"Q_df":n-1 if n>=2 else np.nan,
                     "min_F":g["F_statistic"].min(),"median_F":g["F_statistic"].median(),"mean_F":g["F_statistic"].mean()})

    a=pd.DataFrame(rows)
    a["FDR_MR_tissue"]=bh(a["p_MR"]); a=a.sort_values(["FDR_MR_tissue","p_MR"])
    a.to_csv(f"{out}/{tissue}_IndependentSNP_MR_all.tsv.gz",sep="\t",index=False,compression="gzip")
    a[a["method"]=="IVW"].to_csv(f"{out}/{tissue}_IVW_MR.tsv.gz",sep="\t",index=False,compression="gzip")
    a[a["p_MR"]<.05].to_csv(f"{out}/{tissue}_IndependentSNP_MR_nominal_P0.05.tsv",sep="\t",index=False)
    a[a["FDR_MR_tissue"]<.05].to_csv(f"{out}/{tissue}_IndependentSNP_MR_FDR0.05.tsv",sep="\t",index=False)

    summaries.append({"tissue":tissue,
        "lead_Wald_tests":len(lead_primary),"lead_Wald_nominal_P0.05":int((lead_primary["p_MR"]<.05).sum()),"lead_Wald_FDR0.05":int((lead_primary["FDR_MR_tissue"]<.05).sum()),
        "independent_gene_tests":len(a),"Wald_1SNP_genes":int((a["method"]=="Wald").sum()),"IVW_2plusSNP_genes":int((a["method"]=="IVW").sum()),
        "independent_nominal_P0.05":int((a["p_MR"]<.05).sum()),"independent_FDR0.05":int((a["FDR_MR_tissue"]<.05).sum())})

    lead_all.append(d); adaptive_all.append(a)
    print(f"Lead-SNP Wald: {len(lead_primary):,} tests; P<0.05={(lead_primary['p_MR']<.05).sum():,}; FDR<0.05={(lead_primary['FDR_MR_tissue']<.05).sum():,}")
    print(f"Independent-SNP MR: {len(a):,} genes; Wald={sum(a['method']=='Wald'):,}; IVW={sum(a['method']=='IVW'):,}; P<0.05={(a['p_MR']<.05).sum():,}; FDR<0.05={(a['FDR_MR_tissue']<.05).sum():,}")

# ==============================================================
# GLOBAL FDR: CALCULATE SEPARATELY FOR ORIGINAL WALD AND NEW ANALYSIS
# ==============================================================
lead=pd.concat(lead_all,ignore_index=True); lead["FDR_MR_global"]=bh(lead["p_MR"])
adaptive=pd.concat(adaptive_all,ignore_index=True); adaptive["FDR_MR_global"]=bh(adaptive["p_MR"])

for tissue in TISSUES:
    out=f"{BASE}/{tissue}/MR"
    x=lead[(lead["tissue"]==tissue)&(lead["MR_status"]=="PASS")].sort_values(["FDR_MR_global","p_MR"])
    x.to_csv(f"{out}/{tissue}_Wald_MR_with_global_FDR.tsv.gz",sep="\t",index=False,compression="gzip")
    x[x["FDR_MR_global"]<.05].to_csv(f"{out}/{tissue}_Wald_MR_global_FDR0.05.tsv",sep="\t",index=False)

    y=adaptive[adaptive["tissue"]==tissue].sort_values(["FDR_MR_global","p_MR"])
    y.to_csv(f"{out}/{tissue}_IndependentSNP_MR_with_global_FDR.tsv.gz",sep="\t",index=False,compression="gzip")
    y[y["FDR_MR_global"]<.05].to_csv(f"{out}/{tissue}_IndependentSNP_MR_global_FDR0.05.tsv",sep="\t",index=False)

summary=pd.DataFrame(summaries)
summary.to_csv(f"{GLOBAL}/MR_summary_all_tissues.tsv",sep="\t",index=False)
lead.to_csv(f"{GLOBAL}/all_tissues_lead_Wald_MR.tsv.gz",sep="\t",index=False,compression="gzip")
adaptive.to_csv(f"{GLOBAL}/all_tissues_independentSNP_MR.tsv.gz",sep="\t",index=False,compression="gzip")

print(f"\n{'='*70}\nALL-TISSUE MR SUMMARY\n{'='*70}")
print(summary.to_string(index=False))
print(f"\nLead Wald global FDR<0.05: {(lead['FDR_MR_global']<.05).sum():,}")
print(f"Independent-SNP global FDR<0.05: {(adaptive['FDR_MR_global']<.05).sum():,}")
PY

echo
echo "MR analysis complete."
column -t -s $'\t' "$GLOBAL/MR_summary_all_tissues.tsv"
