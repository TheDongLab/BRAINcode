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
python3 -c 'import pandas' 2>/dev/null || { echo "ERROR: pandas unavailable"; exit 1; }

BASE="$HOME/donglab/data/target_ALS"
HARMONIZED_ROOT="$HOME/donglab/data/GCST90027163/GWAS/eQTL_harmonization"
GLOBAL_SUMMARY="$BASE/MR/MR_summary_all_tissues.tsv"
TISSUES=(Cervical_Spinal_Cord Lumbar_Spinal_Cord Motor_Cortex Frontal_Cortex Cerebellum)
MIN_F=10

export BASE HARMONIZED_ROOT GLOBAL_SUMMARY MIN_F TISSUE_LIST="${TISSUES[*]}"

python3 <<'PY'
import os, math
import pandas as pd

BASE=os.environ["BASE"]
ROOT=os.environ["HARMONIZED_ROOT"]
GLOBAL_SUMMARY=os.environ["GLOBAL_SUMMARY"]
TISSUES=os.environ["TISSUE_LIST"].split()
MIN_F=float(os.environ["MIN_F"])

def bh(p):
    p=pd.to_numeric(p,errors="coerce")
    out=pd.Series(float("nan"),index=p.index,dtype=float)
    good=p.notna()
    if not good.any(): return out
    x=p[good].clip(0,1)
    order=x.sort_values().index
    m=len(order)
    ranked=x.loc[order].values
    adj=[0.0]*m
    running=1.0
    for i in range(m-1,-1,-1):
        running=min(running,ranked[i]*m/(i+1))
        adj[i]=min(running,1.0)
    out.loc[order]=adj
    return out

required=["geneid","snpid","effect_allele","other_allele","beta","se_qtl","p-value","gwas_beta_harmonized","gwas_se","gwas_p","harmonization_keep"]

all_results=[]
summaries=[]

for tissue in TISSUES:
    infile=f"{ROOT}/{tissue}/{tissue}_eQTL_GWAS_MR_ready.tsv.gz"
    outdir=f"{BASE}/{tissue}/MR"
    os.makedirs(outdir,exist_ok=True)

    if not os.path.isfile(infile) or os.path.getsize(infile)==0:
        raise FileNotFoundError(f"Missing or empty input: {infile}")

    print(f"\n{'='*70}")
    print(f"MR: {tissue}")
    print(f"Output: {outdir}")
    print(f"{'='*70}")

    d=pd.read_csv(infile,sep="\t",compression="gzip")

    missing=set(required)-set(d.columns)
    if missing:
        raise ValueError(f"{tissue}: missing columns: {sorted(missing)}")

    d["tissue"]=tissue
    for c in ["beta","se_qtl","p-value","gwas_beta_harmonized","gwas_se","gwas_p"]:
        d[c]=pd.to_numeric(d[c],errors="coerce")

    d=d.rename(columns={"beta":"beta_eqtl","p-value":"p_eqtl"})

    # Instrument strength
    d["F_statistic"]=(d["beta_eqtl"]/d["se_qtl"])**2
    d["strong_instrument"]=d["F_statistic"]>=MIN_F

    # One-SNP Wald ratio
    d["beta_MR"]=d["gwas_beta_harmonized"]/d["beta_eqtl"]

    # First-order SE for reference
    d["se_MR_first_order"]=d["gwas_se"]/d["beta_eqtl"].abs()

    # Delta-method SE using uncertainty in both eQTL and GWAS effects
    d["se_MR_delta"]=(
        (d["gwas_se"]**2/d["beta_eqtl"]**2) +
        (d["gwas_beta_harmonized"]**2*d["se_qtl"]**2/d["beta_eqtl"]**4)
    )**0.5

    d["z_MR"]=d["beta_MR"]/d["se_MR_delta"]
    d["p_MR"]=d["z_MR"].abs().apply(
        lambda z: math.erfc(z/math.sqrt(2)) if pd.notna(z) else float("nan")
    )

    # ALS outcome is on log-odds scale
    d["OR_MR"]=d["beta_MR"].apply(
        lambda x: math.exp(x) if pd.notna(x) and x<700 else float("inf")
    )
    d["CI95_lower_beta"]=d["beta_MR"]-1.96*d["se_MR_delta"]
    d["CI95_upper_beta"]=d["beta_MR"]+1.96*d["se_MR_delta"]
    d["OR_CI95_lower"]=d["CI95_lower_beta"].apply(
        lambda x: math.exp(x) if pd.notna(x) and x<700 else float("inf")
    )
    d["OR_CI95_upper"]=d["CI95_upper_beta"].apply(
        lambda x: math.exp(x) if pd.notna(x) and x<700 else float("inf")
    )

    # BH-FDR within each tissue
    d["FDR_MR_tissue"]=bh(d["p_MR"])

    d["MR_status"]="PASS"
    d.loc[~d["strong_instrument"],"MR_status"]="WEAK_INSTRUMENT"
    d.loc[
        d[["beta_eqtl","se_qtl","gwas_beta_harmonized","gwas_se","beta_MR","se_MR_delta","p_MR"]]
        .isna().any(axis=1),
        "MR_status"
    ]="MISSING_MR_VALUE"
    d.loc[
        (d["se_qtl"]<=0)|(d["gwas_se"]<=0)|(d["se_MR_delta"]<=0),
        "MR_status"
    ]="INVALID_SE"

    all_file=f"{outdir}/{tissue}_Wald_MR_all.tsv.gz"
    primary_file=f"{outdir}/{tissue}_Wald_MR_F{int(MIN_F)}.tsv.gz"
    nominal_file=f"{outdir}/{tissue}_Wald_MR_nominal_P0.05.tsv"
    fdr_file=f"{outdir}/{tissue}_Wald_MR_FDR0.05.tsv"
    summary_file=f"{outdir}/{tissue}_MR_summary.tsv"

    d.to_csv(all_file,sep="\t",index=False,compression="gzip")

    primary=d[d["MR_status"]=="PASS"].copy()
    primary=primary.sort_values(["FDR_MR_tissue","p_MR"])
    primary.to_csv(primary_file,sep="\t",index=False,compression="gzip")

    primary[primary["p_MR"]<0.05].to_csv(
        nominal_file,sep="\t",index=False
    )

    primary[primary["FDR_MR_tissue"]<0.05].to_csv(
        fdr_file,sep="\t",index=False
    )

    n=len(d)
    passed=len(primary)
    strong=int(d["strong_instrument"].sum())
    nominal=int((primary["p_MR"]<0.05).sum())
    fdr05=int((primary["FDR_MR_tissue"]<0.05).sum())
    median_f=d["F_statistic"].median()

    tissue_summary=pd.DataFrame([{
        "tissue":tissue,
        "MR_rows":n,
        "F_ge_10":strong,
        "MR_primary_rows":passed,
        "weak_or_invalid":n-passed,
        "median_F":median_f,
        "nominal_p_lt_0.05":nominal,
        "tissue_FDR_lt_0.05":fdr05
    }])

    tissue_summary.to_csv(
        summary_file,sep="\t",index=False,float_format="%.4g"
    )

    summaries.append(tissue_summary.iloc[0].to_dict())
    all_results.append(d)

    print(f"MR rows: {n:,}")
    print(f"Strong instruments (F >= {MIN_F:g}): {strong:,} ({100*strong/n:.2f}%)")
    print(f"Median F: {median_f:.2f}")
    print(f"Primary MR rows: {passed:,}")
    print(f"Nominal P < 0.05: {nominal:,}")
    print(f"Tissue BH-FDR < 0.05: {fdr05:,}")
    print(f"Primary output: {primary_file}")

# Global FDR across all five tissues
combined=pd.concat(all_results,ignore_index=True)
combined["FDR_MR_global"]=bh(combined["p_MR"])
combined["FDR_MR_global_strong"]=float("nan")

mask=combined["MR_status"]=="PASS"
combined.loc[mask,"FDR_MR_global_strong"]=bh(
    combined.loc[mask,"p_MR"]
)

# Write global-FDR values back into each tissue's MR directory
for tissue in TISSUES:
    outdir=f"{BASE}/{tissue}/MR"
    td=combined[combined["tissue"]==tissue].copy()
    td.to_csv(
        f"{outdir}/{tissue}_Wald_MR_all_with_global_FDR.tsv.gz",
        sep="\t",index=False,compression="gzip"
    )

    strong=td[td["MR_status"]=="PASS"].copy()
    strong=strong.sort_values(["FDR_MR_global_strong","p_MR"])
    strong.to_csv(
        f"{outdir}/{tissue}_Wald_MR_F{int(MIN_F)}_with_global_FDR.tsv.gz",
        sep="\t",index=False,compression="gzip"
    )

    strong[strong["FDR_MR_global_strong"]<0.05].to_csv(
        f"{outdir}/{tissue}_Wald_MR_global_FDR0.05.tsv",
        sep="\t",index=False
    )

summary=pd.DataFrame(summaries)
summary.to_csv(
    GLOBAL_SUMMARY,sep="\t",index=False,float_format="%.4g"
)

print(f"\n{'='*70}")
print("ALL-TISSUE MR SUMMARY")
print(f"{'='*70}")
print(summary.to_string(index=False))

print(f"\nGlobal summary: {GLOBAL_SUMMARY}")
PY

echo
echo "MR analysis complete."
echo
column -t -s $'\t' "$GLOBAL_SUMMARY"
