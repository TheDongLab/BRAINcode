#!/usr/bin/env bash
#SBATCH --job-name=make_eQTL_alleles
#SBATCH --output=/home/zw529/donglab/data/GCST90027163/GWAS/make_eQTL_alleles.out
#SBATCH --error=/home/zw529/donglab/data/GCST90027163/GWAS/make_eQTL_alleles.err
#SBATCH --time=9:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

set -euo pipefail
module purge
module load BCFtools
command -v python3 >/dev/null || { echo "ERROR: python3 not found"; exit 1; }
python3 -c 'import pandas' 2>/dev/null || { echo "ERROR: Python pandas is unavailable"; exit 1; }

PLINK_DIR="$HOME/donglab/data/target_ALS/QTL/plink"; RAW="$PLINK_DIR/joint_all_chrs_matrixEQTL.raw"; BIM="$PLINK_DIR/joint_all_chrs_filtered_bed.bim"
GWAS="$HOME/donglab/data/GCST90027163/GWAS/harmonised/34873335-GCST90027163-MONDO_0004976.h.tsv.gz"
GWAS_OUT="$HOME/donglab/data/GCST90027163/GWAS/eQTL_harmonization"; SUMMARY="$GWAS_OUT/all_tissues_retention_summary.tsv"
TISSUES=(Cervical_Spinal_Cord Lumbar_Spinal_Cord Motor_Cortex Frontal_Cortex Cerebellum)

mkdir -p "$GWAS_OUT"
for f in "$RAW" "$BIM" "$GWAS"; do [[ -s "$f" ]] || { echo "ERROR: missing or empty file: $f"; exit 1; }; done

export RAW BIM GWAS GWAS_OUT SUMMARY
export TISSUE_LIST="${TISSUES[*]}"

python3 <<'PY'
import os
import pandas as pd

raw=os.environ["RAW"]; bim_path=os.environ["BIM"]; gwas_path=os.environ["GWAS"]; out_root=os.environ["GWAS_OUT"]; summary_path=os.environ["SUMMARY"]
tissues=os.environ["TISSUE_LIST"].split(); base=os.path.expanduser("~/donglab/data/target_ALS")

print("Reading RAW genotype header...")
with open(raw) as f: raw_cols=f.readline().split()
variant_cols=[x for x in raw_cols if x not in {"FID","IID","PAT","MAT","SEX","PHENOTYPE"}]
records=[]
for col in variant_cols:
    if "_" not in col: records.append((col,None,"missing_counted_allele_suffix"))
    else:
        variant_id,counted=col.rsplit("_",1); records.append((variant_id,counted.upper(),"ok"))
raw_meta=pd.DataFrame(records,columns=["raw_variant_id","effect_allele","raw_parse_status"])

print("Reading BIM...")
bim=pd.read_csv(bim_path,sep=r"\s+",header=None,names=["bim_chr","bim_id","cm","bim_pos","bim_a1","bim_a2"],dtype=str)
bim["bim_a1"]=bim["bim_a1"].str.upper(); bim["bim_a2"]=bim["bim_a2"].str.upper()
dat=raw_meta.merge(bim,left_on="raw_variant_id",right_on="bim_id",how="left",validate="one_to_one")
dat["other_allele"]=dat.apply(lambda r:r.bim_a2 if r.effect_allele==r.bim_a1 else r.bim_a1 if r.effect_allele==r.bim_a2 else None,axis=1)
dat["allele_status"]=dat.apply(lambda r:"ok" if r.effect_allele in {r.bim_a1,r.bim_a2} else "missing_bim_match" if pd.isna(r.bim_id) else "counted_allele_not_in_bim",axis=1)
dat["key"]=dat["bim_chr"].str.replace("^chr","",regex=True).str.upper()+":"+dat["bim_pos"]

print("Reading harmonized GRCh38 GWAS...")
gwas_cols=["hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","standard_error","p_value"]
gwas=pd.read_csv(gwas_path,sep="\t",compression="gzip",usecols=gwas_cols,dtype={"hm_rsid":str,"hm_chrom":str,"hm_pos":str,"hm_other_allele":str,"hm_effect_allele":str})
gwas=gwas.rename(columns={"hm_rsid":"snpid","hm_chrom":"gwas_chr","hm_pos":"gwas_pos","hm_other_allele":"gwas_other_allele","hm_effect_allele":"gwas_effect_allele","hm_beta":"gwas_beta_original","standard_error":"gwas_se","p_value":"gwas_p"})
gwas["gwas_effect_allele"]=gwas["gwas_effect_allele"].str.upper(); gwas["gwas_other_allele"]=gwas["gwas_other_allele"].str.upper()
gwas=gwas.drop_duplicates("snpid")

comp={"A":"T","T":"A","C":"G","G":"C"}
summaries=[]

for tissue in tissues:
    eqtl_dir=f"{base}/{tissue}/eQTL"; leads=f"{eqtl_dir}/results/{tissue}_eQTL.lead_snps.txt"; snp_loc=f"{eqtl_dir}/snp_location.txt"
    tissue_out=f"{out_root}/{tissue}"; allele_out=f"{eqtl_dir}/qtl_alleles_GRCh38.tsv"
    all_out=f"{tissue_out}/{tissue}_eQTL_GWAS_all.tsv.gz"; kept_out=f"{tissue_out}/{tissue}_eQTL_GWAS_MR_ready.tsv.gz"; rejected_out=f"{tissue_out}/{tissue}_eQTL_GWAS_rejected.tsv.gz"
    os.makedirs(tissue_out,exist_ok=True)
    for f in [leads,snp_loc]:
        if not os.path.isfile(f) or os.path.getsize(f)==0: raise FileNotFoundError(f"Missing or empty file: {f}")

    print(f"\n{'='*70}\nProcessing {tissue}\n{'='*70}")
    loc=pd.read_csv(snp_loc,sep="\t",dtype=str).rename(columns={"snpid":"final_snpid","chr":"final_chr","pos":"final_pos"})
    loc["key"]=loc["final_chr"].str.replace("^chr","",regex=True).str.upper()+":"+loc["final_pos"]
    tissue_dat=dat.merge(loc[["key","final_snpid","final_chr","final_pos"]],on="key",how="left")

    lead=pd.read_csv(leads,sep="\t",dtype=str)
    lead["beta"]=pd.to_numeric(lead["beta"],errors="coerce"); lead["t-stat"]=pd.to_numeric(lead["t-stat"],errors="coerce")
    lead["se_qtl"]=(lead["beta"]/lead["t-stat"]).abs()
    lead_ids=set(lead["snpid"])

    alleles=tissue_dat[tissue_dat["final_snpid"].isin(lead_ids)].copy()
    alleles=alleles[["final_snpid","final_chr","final_pos","effect_allele","other_allele","bim_a1","bim_a2","raw_variant_id","allele_status"]]
    alleles.columns=["snpid","chr","pos","effect_allele","other_allele","bim_A1","bim_A2","raw_variant_id","allele_status"]
    alleles=alleles.drop_duplicates("snpid")
    alleles.to_csv(allele_out,sep="\t",index=False)

    qtl=lead.merge(alleles,on="snpid",how="left")
    x=qtl.merge(gwas,on="snpid",how="left")
    x["found_in_gwas"]=x["gwas_chr"].notna()
    x["coordinate_match"]=x["found_in_gwas"]&(x["chr"].str.replace("^chr","",regex=True).str.upper()==x["gwas_chr"].str.replace("^chr","",regex=True).str.upper())&(x["pos"].astype(str)==x["gwas_pos"].astype(str))
    x["effect_allele"]=x["effect_allele"].str.upper(); x["other_allele"]=x["other_allele"].str.upper()
    x["palindromic"]=(x["effect_allele"]+x["other_allele"]).isin(["AT","TA","CG","GC"])

    def relation(r):
        if not r["found_in_gwas"]: return "not_found_in_GWAS"
        qe,qo,ge,go=r["effect_allele"],r["other_allele"],r["gwas_effect_allele"],r["gwas_other_allele"]
        if pd.isna(qe) or pd.isna(qo): return "missing_QTL_alleles"
        if qe==ge and qo==go: return "exact"
        if qe==go and qo==ge: return "swap_beta"
        if qe==comp.get(ge) and qo==comp.get(go): return "strand"
        if qe==comp.get(go) and qo==comp.get(ge): return "strand_and_swap_beta"
        return "incompatible"

    x["harmonization_relation"]=x.apply(relation,axis=1)
    x["beta_flip_required"]=x["harmonization_relation"].isin(["swap_beta","strand_and_swap_beta"])
    x["strand_flip_required"]=x["harmonization_relation"].isin(["strand","strand_and_swap_beta"])
    x["gwas_beta_harmonized"]=pd.to_numeric(x["gwas_beta_original"],errors="coerce")
    x.loc[x["beta_flip_required"],"gwas_beta_harmonized"]*=-1

    x["rejection_reason"]=""
    x.loc[x["harmonization_relation"]=="not_found_in_GWAS","rejection_reason"]="not_found_in_GWAS"
    x.loc[x["harmonization_relation"]=="missing_QTL_alleles","rejection_reason"]="missing_QTL_alleles"
    x.loc[x["harmonization_relation"]=="incompatible","rejection_reason"]="incompatible_alleles"
    x.loc[(x["rejection_reason"]=="")&x["palindromic"],"rejection_reason"]="palindromic_without_frequency_check"
    x.loc[(x["rejection_reason"]=="")&~x["coordinate_match"],"rejection_reason"]="GRCh38_coordinate_mismatch"
    x.loc[(x["rejection_reason"]=="")&(x["allele_status"]!="ok"),"rejection_reason"]="QTL_allele_status_not_ok"
    x.loc[(x["rejection_reason"]=="")&(x["gwas_beta_harmonized"].isna()|x["gwas_se"].isna()|x["se_qtl"].isna()),"rejection_reason"]="missing_beta_or_standard_error"
    x["harmonization_keep"]=x["rejection_reason"]==""

    x.to_csv(all_out,sep="\t",index=False,compression="gzip")
    x[x["harmonization_keep"]].to_csv(kept_out,sep="\t",index=False,compression="gzip")
    x[~x["harmonization_keep"]].to_csv(rejected_out,sep="\t",index=False,compression="gzip")

    requested=len(lead_ids); allele_records=len(alleles); matched=int(x["found_in_gwas"].sum()); palindromic=int((x["found_in_gwas"]&x["palindromic"]).sum()); retained=int(x["harmonization_keep"].sum())
    summaries.append({"tissue":tissue,"lead_eQTLs_requested":requested,"allele_records":allele_records,"matched_in_GWAS":matched,"matched_palindromic":palindromic,"retained_MR_ready":retained,"allele_recovery_pct":100*allele_records/requested if requested else 0,"GWAS_match_pct":100*matched/requested if requested else 0,"retention_pct":100*retained/requested if requested else 0})

    print(f"Lead eQTLs requested: {requested:,}")
    print(f"Allele records recovered: {allele_records:,} ({100*allele_records/requested:.2f}%)")
    print(f"Matched in GWAS: {matched:,} ({100*matched/requested:.2f}%)")
    print(f"Matched palindromic excluded: {palindromic:,}")
    print(f"MR-ready retained: {retained:,} ({100*retained/requested:.2f}%)")
    print("Harmonization relationships:")
    print(x["harmonization_relation"].value_counts(dropna=False).to_string())
    print(f"MR-ready output: {kept_out}")

summary=pd.DataFrame(summaries)
summary.to_csv(summary_path,sep="\t",index=False,float_format="%.2f")
print(f"\n{'='*70}\nALL-TISSUE RETENTION SUMMARY\n{'='*70}")
print(summary.to_string(index=False))
print(f"\nSummary written to: {summary_path}")
PY

echo
echo "Finished all tissues."
column -t -s $'\t' "$SUMMARY"
