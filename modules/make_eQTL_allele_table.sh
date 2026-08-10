#!/usr/bin/env bash
#SBATCH --job-name=make_eQTL_alleles
#SBATCH --output=/home/zw529/donglab/data/GCST90027163/GWAS/make_eQTL_alleles.out
#SBATCH --error=/home/zw529/donglab/data/GCST90027163/GWAS/make_eQTL_alleles.err
#SBATCH --time=4:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4

set -euo pipefail
module --force purge
module load PLINK/1.9b_7.11-x86_64
command -v python3 >/dev/null || { echo "ERROR: python3 not found"; exit 1; }
python3 -c 'import pandas,numpy' 2>/dev/null || { echo "ERROR: pandas/numpy unavailable"; exit 1; }

BASE="$HOME/donglab/data/target_ALS"; PLINK_DIR="$BASE/QTL/plink"; BFILE="$PLINK_DIR/joint_all_chrs_filtered_bed"; RAW="$PLINK_DIR/joint_all_chrs_matrixEQTL.raw"; BIM="$BFILE.bim"
GWAS="$HOME/donglab/data/GCST90027163/GWAS/harmonised/34873335-GCST90027163-MONDO_0004976.h.tsv.gz"; OUTROOT="$HOME/donglab/data/GCST90027163/GWAS/eQTL_harmonization"
WALD_SUMMARY="$OUTROOT/all_tissues_retention_summary.tsv"; MULTI_SUMMARY="$OUTROOT/all_tissues_multiSNP_summary.tsv"
TISSUES=(Cervical_Spinal_Cord Lumbar_Spinal_Cord Motor_Cortex Frontal_Cortex Cerebellum); MIN_F=10; LD_R2=0.01

mkdir -p "$OUTROOT"
for f in "$RAW" "$BFILE.bed" "$BFILE.bim" "$BFILE.fam" "$GWAS"; do [[ -s "$f" ]] || { echo "ERROR: missing/empty: $f"; exit 1; }; done
export BASE BFILE RAW BIM GWAS OUTROOT WALD_SUMMARY MULTI_SUMMARY MIN_F LD_R2 TISSUE_LIST="${TISSUES[*]}" PLINK_BIN="$(command -v plink)"

python3 <<'PY'
import os, subprocess, tempfile
import numpy as np, pandas as pd

BASE,BFILE,RAW,BIM,GWAS,OUTROOT=[os.environ[x] for x in ["BASE","BFILE","RAW","BIM","GWAS","OUTROOT"]]
WALD_SUMMARY,MULTI_SUMMARY,PLINK=os.environ["WALD_SUMMARY"],os.environ["MULTI_SUMMARY"],os.environ["PLINK_BIN"]
TISSUES=os.environ["TISSUE_LIST"].split(); MIN_F=float(os.environ["MIN_F"]); LD_R2=float(os.environ["LD_R2"])
META={"FID","IID","PAT","MAT","SEX","PHENOTYPE"}; COMP={"A":"T","T":"A","C":"G","G":"C"}

def norm_chr(s): return s.astype(str).str.replace(r"^chr","",regex=True,case=False).str.upper()
def pct(n,d): return 100*n/d if d else 0
def req(df,cols,label):
    m=set(cols)-set(df.columns)
    if m: raise ValueError(f"{label} missing columns: {sorted(m)}")
def qtl_stats(d):
    d=d.copy()
    for c in ["beta","t-stat","p-value"]: d[c]=pd.to_numeric(d[c],errors="coerce")
    d["se_qtl"]=(d["beta"]/d["t-stat"]).abs(); d["F_statistic"]=(d["beta"]/d["se_qtl"])**2
    return d
def harmonize(qtl,alleles,gwas):
    x=qtl.merge(alleles,on="snpid",how="left",validate="many_to_one").merge(gwas,on="snpid",how="left",validate="many_to_one")
    x["found_in_gwas"]=x["gwas_chr"].notna()
    x["allele_coordinate_match"]=(norm_chr(x["qtl_chr"])==norm_chr(x["allele_chr"]))&(x["qtl_pos"].astype(str)==x["allele_pos"].astype(str))
    x["gwas_coordinate_match"]=x["found_in_gwas"]&(norm_chr(x["qtl_chr"])==norm_chr(x["gwas_chr"]))&(x["qtl_pos"].astype(str)==x["gwas_pos"].astype(str))
    x["effect_allele"]=x["effect_allele"].str.upper(); x["other_allele"]=x["other_allele"].str.upper()
    x["palindromic"]=(x["effect_allele"].fillna("")+x["other_allele"].fillna("")).isin(["AT","TA","CG","GC"])
    def rel(r):
        if not r["found_in_gwas"]: return "not_found_in_GWAS"
        qe,qo,ge,go=r["effect_allele"],r["other_allele"],r["gwas_effect_allele"],r["gwas_other_allele"]
        if pd.isna(qe) or pd.isna(qo): return "missing_QTL_alleles"
        if pd.isna(ge) or pd.isna(go): return "missing_GWAS_alleles"
        if qe==ge and qo==go: return "exact"
        if qe==go and qo==ge: return "swap_beta"
        if qe==COMP.get(ge) and qo==COMP.get(go): return "strand"
        if qe==COMP.get(go) and qo==COMP.get(ge): return "strand_and_swap_beta"
        return "incompatible"
    x["harmonization_relation"]=x.apply(rel,axis=1)
    x["beta_flip_required"]=x["harmonization_relation"].isin(["swap_beta","strand_and_swap_beta"]); x["strand_flip_required"]=x["harmonization_relation"].isin(["strand","strand_and_swap_beta"])
    x["gwas_beta_harmonized"]=x["gwas_beta_original"]; x.loc[x["beta_flip_required"],"gwas_beta_harmonized"]=-x.loc[x["beta_flip_required"],"gwas_beta_original"]
    x["harmonized_effect_allele"]=x["effect_allele"]; x["harmonized_other_allele"]=x["other_allele"]; x["rejection_reason"]=""
    x.loc[x["harmonization_relation"]=="not_found_in_GWAS","rejection_reason"]="not_found_in_GWAS"
    x.loc[x["harmonization_relation"]=="missing_QTL_alleles","rejection_reason"]="missing_QTL_alleles"
    x.loc[x["harmonization_relation"]=="missing_GWAS_alleles","rejection_reason"]="missing_GWAS_alleles"
    x.loc[x["harmonization_relation"]=="incompatible","rejection_reason"]="incompatible_alleles"
    x.loc[(x["rejection_reason"]=="")&x["palindromic"],"rejection_reason"]="palindromic_without_frequency_check"
    x.loc[(x["rejection_reason"]=="")&~x["allele_coordinate_match"],"rejection_reason"]="QTL_allele_coordinate_mismatch"
    x.loc[(x["rejection_reason"]=="")&~x["gwas_coordinate_match"],"rejection_reason"]="GRCh38_coordinate_mismatch"
    x.loc[(x["rejection_reason"]=="")&(x["allele_status"]!="ok"),"rejection_reason"]="QTL_allele_status_not_ok"
    x.loc[(x["rejection_reason"]=="")&(x["beta"].isna()|x["se_qtl"].isna()|x["gwas_beta_harmonized"].isna()|x["gwas_se"].isna()),"rejection_reason"]="missing_beta_or_standard_error"
    x.loc[(x["rejection_reason"]=="")&((x["se_qtl"]<=0)|(x["gwas_se"]<=0)),"rejection_reason"]="nonpositive_standard_error"
    x["harmonization_keep"]=x["rejection_reason"]==""
    return x

print("Reading RAW header...")
with open(RAW) as f: cols=f.readline().split()
records=[]
for c in [x for x in cols if x not in META]:
    if "_" not in c: records.append((c,None,"missing_counted_allele_suffix"))
    else:
        vid,a=c.rsplit("_",1); records.append((vid,a.upper(),"ok"))
raw_meta=pd.DataFrame(records,columns=["raw_variant_id","effect_allele","raw_parse_status"]); print(f"RAW genotype columns: {len(raw_meta):,}")

print("Reading BIM...")
bim=pd.read_csv(BIM,sep=r"\s+",header=None,names=["bim_chr","bim_id","cm","bim_pos","bim_a1","bim_a2"],dtype=str)
bim["bim_a1"]=bim["bim_a1"].str.upper(); bim["bim_a2"]=bim["bim_a2"].str.upper()
vm=raw_meta.merge(bim,left_on="raw_variant_id",right_on="bim_id",how="left",validate="one_to_one")
vm["other_allele"]=vm.apply(lambda r:r["bim_a2"] if r["effect_allele"]==r["bim_a1"] else r["bim_a1"] if r["effect_allele"]==r["bim_a2"] else None,axis=1)
vm["allele_status"]=vm.apply(lambda r:"ok" if r["effect_allele"] in {r["bim_a1"],r["bim_a2"]} else "missing_bim_match" if pd.isna(r["bim_id"]) else "counted_allele_not_in_bim",axis=1)
vm["key"]=norm_chr(vm["bim_chr"])+":"+vm["bim_pos"].astype(str)

print("Reading GWAS...")
gc=["hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","standard_error","p_value"]
gwas=pd.read_csv(GWAS,sep="\t",compression="gzip",usecols=gc,dtype={"hm_rsid":str,"hm_chrom":str,"hm_pos":str,"hm_other_allele":str,"hm_effect_allele":str})
gwas=gwas.rename(columns={"hm_rsid":"snpid","hm_chrom":"gwas_chr","hm_pos":"gwas_pos","hm_other_allele":"gwas_other_allele","hm_effect_allele":"gwas_effect_allele","hm_beta":"gwas_beta_original","standard_error":"gwas_se","p_value":"gwas_p"})
gwas["gwas_effect_allele"]=gwas["gwas_effect_allele"].str.upper(); gwas["gwas_other_allele"]=gwas["gwas_other_allele"].str.upper()
for c in ["gwas_beta_original","gwas_se","gwas_p"]: gwas[c]=pd.to_numeric(gwas[c],errors="coerce")
gwas=gwas.drop_duplicates("snpid")

wald_summary=[]; multi_summary=[]

for tissue in TISSUES:
    E=f"{BASE}/{tissue}/eQTL"; R=f"{E}/results"; O=f"{OUTROOT}/{tissue}"; os.makedirs(O,exist_ok=True)
    leads_path=f"{R}/{tissue}_eQTL.lead_snps.txt"; fdr_path=f"{R}/{tissue}_eQTL.FDR0.05.txt"; loc_path=f"{E}/snp_location.txt"
    for f in [leads_path,fdr_path,loc_path]:
        if not os.path.isfile(f) or os.path.getsize(f)==0: raise FileNotFoundError(f"Missing/empty: {f}")
    print(f"\n{'='*70}\n{tissue}\n{'='*70}")

    loc=pd.read_csv(loc_path,sep="\t",dtype=str); req(loc,["snpid","chr","pos"],f"{tissue} locations"); loc=loc.drop_duplicates("snpid"); loc["key"]=norm_chr(loc["chr"])+":"+loc["pos"].astype(str)
    leads=pd.read_csv(leads_path,sep="\t",dtype=str); req(leads,["geneid","snpid","beta","t-stat","p-value","chr","pos"],f"{tissue} leads"); leads=leads.rename(columns={"chr":"qtl_chr","pos":"qtl_pos"}); leads=qtl_stats(leads)
    fdr=pd.read_csv(fdr_path,sep="\t",dtype=str); req(fdr,["geneid","snpid","beta","t-stat","p-value","FDR"],f"{tissue} FDR0.05")
    fdr=fdr.merge(loc[["snpid","chr","pos"]],on="snpid",how="left",validate="many_to_one").rename(columns={"chr":"qtl_chr","pos":"qtl_pos"}); fdr=qtl_stats(fdr)

    needed=set(leads["snpid"].dropna())|set(fdr["snpid"].dropna())
    alleles=loc[loc["snpid"].isin(needed)].merge(vm[["key","effect_allele","other_allele","bim_a1","bim_a2","raw_variant_id","allele_status"]],on="key",how="left")
    alleles=alleles.rename(columns={"chr":"allele_chr","pos":"allele_pos","bim_a1":"bim_A1","bim_a2":"bim_A2"})[["snpid","allele_chr","allele_pos","effect_allele","other_allele","bim_A1","bim_A2","raw_variant_id","allele_status"]].drop_duplicates("snpid")

    # ---------------- WALD BRANCH: unchanged lead SNP path ----------------
    lead_a=alleles[alleles["snpid"].isin(set(leads["snpid"]))].copy()
    lead_a.to_csv(f"{E}/qtl_alleles_GRCh38.tsv",sep="\t",index=False)
    wald=harmonize(leads,lead_a,gwas); wald_ready=wald[wald["harmonization_keep"]].copy()
    wald.to_csv(f"{O}/{tissue}_eQTL_GWAS_all.tsv.gz",sep="\t",index=False,compression="gzip")
    wald_ready.to_csv(f"{O}/{tissue}_eQTL_GWAS_MR_ready.tsv.gz",sep="\t",index=False,compression="gzip")
    wald[~wald["harmonization_keep"]].to_csv(f"{O}/{tissue}_eQTL_GWAS_rejected.tsv.gz",sep="\t",index=False,compression="gzip")
    wald_summary.append({"tissue":tissue,"lead_gene_snp_rows":len(leads),"unique_lead_snps":leads["snpid"].nunique(),"retained_MR_ready":len(wald_ready),"retention_pct":pct(len(wald_ready),len(leads))})
    print(f"Wald: {len(leads):,} -> {len(wald_ready):,} MR-ready ({pct(len(wald_ready),len(leads)):.2f}%)")

    # ---------------- MULTI-SNP / IVW BRANCH ----------------
    fdr_a=alleles[alleles["snpid"].isin(set(fdr["snpid"]))].copy()
    multi=harmonize(fdr,fdr_a,gwas); multi["strong_instrument"]=multi["F_statistic"]>=MIN_F
    pre=multi[multi["harmonization_keep"]&multi["strong_instrument"]].sort_values(["geneid","p-value","snpid"]).drop_duplicates(["geneid","raw_variant_id"]).copy()
    multi.to_csv(f"{O}/{tissue}_eQTL_GWAS_multiSNP_all.tsv.gz",sep="\t",index=False,compression="gzip")

    if len(pre)==0:
        print("No multi-SNP candidates survived harmonization/F threshold."); continue

    with tempfile.TemporaryDirectory(prefix=f"{tissue}_LD_",dir=OUTROOT) as tmp:
        ids=os.path.join(tmp,"ids.txt"); prefix=os.path.join(tmp,"geno")
        with open(ids,"w") as fh:
            for v in sorted(set(pre["raw_variant_id"].dropna())): fh.write(v+"\n")
        subprocess.run([PLINK,"--bfile",BFILE,"--extract",ids,"--recode","A-transpose","--threads",os.environ.get("SLURM_CPUS_PER_TASK","1"),"--out",prefix],check=True)
        traw=pd.read_csv(prefix+".traw",sep="\t",low_memory=False); req(traw,["SNP"],f"{tissue} .traw")
        geno_ids=traw["SNP"].astype(str).tolist(); G=traw.iloc[:,6:].to_numpy(dtype=np.float32); del traw
        means=np.nanmean(G,axis=1); means[np.isnan(means)]=0; rr,cc=np.where(np.isnan(G))
        if len(rr): G[rr,cc]=means[rr]
        G-=means[:,None]; norms=np.sqrt(np.sum(G*G,axis=1)); variable=norms>0; G[variable]/=norms[variable,None]; G[~variable]=0
        idx={v:i for i,v in enumerate(geno_ids)}

        pre["ld_selected"]=False; pre["max_r2_to_selected"]=np.nan; pre["ld_prune_reason"]=""
        for gene,g in pre.groupby("geneid",sort=False):
            selected=[]
            for ri,r in g.sort_values(["p-value","snpid"]).iterrows():
                v=str(r["raw_variant_id"])
                if v not in idx: pre.at[ri,"ld_prune_reason"]="missing_PLINK"; continue
                gi=idx[v]
                if not variable[gi]: pre.at[ri,"ld_prune_reason"]="zero_variance"; continue
                if not selected:
                    pre.at[ri,"ld_selected"]=True; pre.at[ri,"ld_prune_reason"]="selected"; selected.append(gi); continue
                r2=(G[selected]@G[gi])**2; maxr2=float(np.max(r2)); pre.at[ri,"max_r2_to_selected"]=maxr2
                if maxr2<LD_R2: pre.at[ri,"ld_selected"]=True; pre.at[ri,"ld_prune_reason"]="selected"; selected.append(gi)
                else: pre.at[ri,"ld_prune_reason"]=f"LD_r2_ge_{LD_R2}"

    selected=pre[pre["ld_selected"]].copy()
    counts=selected.groupby("geneid").size(); selected["n_independent_instruments"]=selected["geneid"].map(counts).astype(int)
    selected["recommended_MR_method"]=np.where(selected["n_independent_instruments"]>=2,"IVW","Wald")
    ivw=selected[selected["n_independent_instruments"]>=2].copy()

    pre.to_csv(f"{O}/{tissue}_eQTL_GWAS_multiSNP_LD_pruned.tsv.gz",sep="\t",index=False,compression="gzip")
    selected.to_csv(f"{O}/{tissue}_eQTL_GWAS_multiSNP_MR_ready.tsv.gz",sep="\t",index=False,compression="gzip")
    ivw.to_csv(f"{O}/{tissue}_eQTL_GWAS_IVW_ready.tsv.gz",sep="\t",index=False,compression="gzip")

    gs=pd.DataFrame({"n_before_LD":pre.groupby("geneid").size(),"n_after_LD":selected.groupby("geneid").size()}).fillna(0).astype(int).reset_index()
    gs["recommended_MR_method"]=np.where(gs["n_after_LD"]>=2,"IVW",np.where(gs["n_after_LD"]==1,"Wald","No_instrument"))
    gs.to_csv(f"{O}/{tissue}_multiSNP_gene_summary.tsv",sep="\t",index=False)

    one=int((gs["n_after_LD"]==1).sum()); nivw=int((gs["n_after_LD"]>=2).sum())
    multi_summary.append({"tissue":tissue,"FDR0.05_rows":len(fdr),"FDR0.05_genes":fdr["geneid"].nunique(),"strong_harmonized_before_LD":len(pre),"independent_SNP_rows":len(selected),"genes_with_1_SNP":one,"genes_with_2plus_SNPs_IVW":nivw,"LD_r2_threshold":LD_R2,"minimum_F":MIN_F})
    print(f"Multi-SNP: {len(pre):,} pre-LD -> {len(selected):,} independent SNP rows; Wald genes={one:,}; IVW genes={nivw:,}")

pd.DataFrame(wald_summary).to_csv(WALD_SUMMARY,sep="\t",index=False,float_format="%.2f")
pd.DataFrame(multi_summary).to_csv(MULTI_SUMMARY,sep="\t",index=False,float_format="%.4g")

print("\nWALD SUMMARY")
print(pd.DataFrame(wald_summary).to_string(index=False))
print("\nMULTI-SNP SUMMARY")
print(pd.DataFrame(multi_summary).to_string(index=False))
PY

echo
echo "=== WALD SUMMARY ==="
column -t -s $'\t' "$WALD_SUMMARY"
echo
echo "=== MULTI-SNP / IVW SUMMARY ==="
column -t -s $'\t' "$MULTI_SUMMARY"
