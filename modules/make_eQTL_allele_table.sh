#!/usr/bin/env bash
#SBATCH --job-name=make_eQTL_alleles
#SBATCH --output=/home/zw529/donglab/data/GCST90027163/GWAS/make_eQTL_alleles.out
#SBATCH --error=/home/zw529/donglab/data/GCST90027163/GWAS/make_eQTL_alleles.err
#SBATCH --time=9:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

set -euo pipefail
module --force purge
module load BCFtools
command -v python3 >/dev/null || { echo "ERROR: python3 not found"; exit 1; }
python3 -c 'import pandas' 2>/dev/null || { echo "ERROR: pandas unavailable"; exit 1; }

PLINK_DIR="$HOME/donglab/data/target_ALS/QTL/plink"; RAW="$PLINK_DIR/joint_all_chrs_matrixEQTL.raw"; BIM="$PLINK_DIR/joint_all_chrs_filtered_bed.bim"
GWAS="$HOME/donglab/data/GCST90027163/GWAS/harmonised/34873335-GCST90027163-MONDO_0004976.h.tsv.gz"
OUTROOT="$HOME/donglab/data/GCST90027163/GWAS/eQTL_harmonization"; SUMMARY="$OUTROOT/all_tissues_retention_summary.tsv"
TISSUES=(Cervical_Spinal_Cord Lumbar_Spinal_Cord Motor_Cortex Frontal_Cortex Cerebellum)

mkdir -p "$OUTROOT"
for f in "$RAW" "$BIM" "$GWAS"; do [[ -s "$f" ]] || { echo "ERROR: missing or empty file: $f"; exit 1; }; done
export RAW BIM GWAS OUTROOT SUMMARY TISSUE_LIST="${TISSUES[*]}"

python3 <<'PY'
import os
import pandas as pd

raw,bim_path,gwas_path,outroot,summary_path=[os.environ[x] for x in ["RAW","BIM","GWAS","OUTROOT","SUMMARY"]]
tissues=os.environ["TISSUE_LIST"].split(); base=os.path.expanduser("~/donglab/data/target_ALS")
meta_cols={"FID","IID","PAT","MAT","SEX","PHENOTYPE"}; comp={"A":"T","T":"A","C":"G","G":"C"}

def norm_chr(s): return s.astype(str).str.replace(r"^chr","",regex=True,case=False).str.upper()
def pct(n,d): return 100*n/d if d else 0
def require_columns(df,required,label):
    missing=set(required)-set(df.columns)
    if missing: raise ValueError(f"{label} missing columns: {sorted(missing)}")

print("Reading RAW genotype header...")
with open(raw) as f: raw_cols=f.readline().split()
variant_cols=[x for x in raw_cols if x not in meta_cols]
records=[]
for col in variant_cols:
    if "_" not in col: records.append((col,None,"missing_counted_allele_suffix"))
    else:
        variant_id,counted=col.rsplit("_",1); records.append((variant_id,counted.upper(),"ok"))
raw_meta=pd.DataFrame(records,columns=["raw_variant_id","effect_allele","raw_parse_status"])

print(f"RAW genotype columns: {len(variant_cols):,}")
print("Reading BIM...")
bim=pd.read_csv(bim_path,sep=r"\s+",header=None,names=["bim_chr","bim_id","cm","bim_pos","bim_a1","bim_a2"],dtype=str)
bim["bim_a1"]=bim["bim_a1"].str.upper(); bim["bim_a2"]=bim["bim_a2"].str.upper()
variant_meta=raw_meta.merge(bim,left_on="raw_variant_id",right_on="bim_id",how="left",validate="one_to_one")
variant_meta["other_allele"]=variant_meta.apply(lambda r:r["bim_a2"] if r["effect_allele"]==r["bim_a1"] else r["bim_a1"] if r["effect_allele"]==r["bim_a2"] else None,axis=1)
variant_meta["allele_status"]=variant_meta.apply(lambda r:"ok" if r["effect_allele"] in {r["bim_a1"],r["bim_a2"]} else "missing_bim_match" if pd.isna(r["bim_id"]) else "counted_allele_not_in_bim",axis=1)
variant_meta["key"]=norm_chr(variant_meta["bim_chr"])+":"+variant_meta["bim_pos"].astype(str)

print("Reading harmonized GRCh38 GWAS...")
gwas_cols=["hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","standard_error","p_value"]
gwas=pd.read_csv(gwas_path,sep="\t",compression="gzip",usecols=gwas_cols,dtype={"hm_rsid":str,"hm_chrom":str,"hm_pos":str,"hm_other_allele":str,"hm_effect_allele":str})
gwas=gwas.rename(columns={"hm_rsid":"snpid","hm_chrom":"gwas_chr","hm_pos":"gwas_pos","hm_other_allele":"gwas_other_allele","hm_effect_allele":"gwas_effect_allele","hm_beta":"gwas_beta_original","standard_error":"gwas_se","p_value":"gwas_p"})
gwas["gwas_effect_allele"]=gwas["gwas_effect_allele"].str.upper(); gwas["gwas_other_allele"]=gwas["gwas_other_allele"].str.upper()
gwas["gwas_beta_original"]=pd.to_numeric(gwas["gwas_beta_original"],errors="coerce"); gwas["gwas_se"]=pd.to_numeric(gwas["gwas_se"],errors="coerce"); gwas["gwas_p"]=pd.to_numeric(gwas["gwas_p"],errors="coerce")
gwas=gwas.drop_duplicates("snpid")

summaries=[]
for tissue in tissues:
    eqtl_dir=f"{base}/{tissue}/eQTL"; leads_path=f"{eqtl_dir}/results/{tissue}_eQTL.lead_snps.txt"; loc_path=f"{eqtl_dir}/snp_location.txt"
    tissue_out=f"{outroot}/{tissue}"; allele_out=f"{eqtl_dir}/qtl_alleles_GRCh38.tsv"
    all_out=f"{tissue_out}/{tissue}_eQTL_GWAS_all.tsv.gz"; kept_out=f"{tissue_out}/{tissue}_eQTL_GWAS_MR_ready.tsv.gz"; rejected_out=f"{tissue_out}/{tissue}_eQTL_GWAS_rejected.tsv.gz"
    os.makedirs(tissue_out,exist_ok=True)
    for f in [leads_path,loc_path]:
        if not os.path.isfile(f) or os.path.getsize(f)==0: raise FileNotFoundError(f"Missing or empty file: {f}")

    print(f"\n{'='*70}\nProcessing {tissue}\n{'='*70}")
    loc=pd.read_csv(loc_path,sep="\t",dtype=str)
    require_columns(loc,["snpid","chr","pos"],f"{tissue} SNP-location file")
    loc=loc.rename(columns={"snpid":"final_snpid","chr":"allele_chr","pos":"allele_pos"})
    loc["key"]=norm_chr(loc["allele_chr"])+":"+loc["allele_pos"].astype(str)

    leads=pd.read_csv(leads_path,sep="\t",dtype=str)
    require_columns(leads,["geneid","snpid","beta","t-stat","p-value","chr","pos"],f"{tissue} lead-eQTL file")
    leads=leads.rename(columns={"chr":"qtl_chr","pos":"qtl_pos"})
    leads["beta"]=pd.to_numeric(leads["beta"],errors="coerce"); leads["t-stat"]=pd.to_numeric(leads["t-stat"],errors="coerce"); leads["p-value"]=pd.to_numeric(leads["p-value"],errors="coerce")
    leads["se_qtl"]=(leads["beta"]/leads["t-stat"]).abs()
    lead_ids=set(leads["snpid"].dropna())

    tissue_variants=variant_meta.merge(loc[["key","final_snpid","allele_chr","allele_pos"]],on="key",how="left")
    alleles=tissue_variants[tissue_variants["final_snpid"].isin(lead_ids)][["final_snpid","allele_chr","allele_pos","effect_allele","other_allele","bim_a1","bim_a2","raw_variant_id","allele_status"]].copy()
    alleles.columns=["snpid","allele_chr","allele_pos","effect_allele","other_allele","bim_A1","bim_A2","raw_variant_id","allele_status"]
    alleles=alleles.drop_duplicates("snpid")
    alleles.to_csv(allele_out,sep="\t",index=False)

    qtl=leads.merge(alleles,on="snpid",how="left",validate="one_to_one")
    x=qtl.merge(gwas,on="snpid",how="left",validate="one_to_one")
    x["found_in_gwas"]=x["gwas_chr"].notna()
    x["allele_coordinate_match"]=(norm_chr(x["qtl_chr"])==norm_chr(x["allele_chr"]))&(x["qtl_pos"].astype(str)==x["allele_pos"].astype(str))
    x["gwas_coordinate_match"]=x["found_in_gwas"]&(norm_chr(x["qtl_chr"])==norm_chr(x["gwas_chr"]))&(x["qtl_pos"].astype(str)==x["gwas_pos"].astype(str))
    x["effect_allele"]=x["effect_allele"].str.upper(); x["other_allele"]=x["other_allele"].str.upper()
    x["palindromic"]=(x["effect_allele"].fillna("")+x["other_allele"].fillna("")).isin(["AT","TA","CG","GC"])

    def relation(r):
        if not r["found_in_gwas"]: return "not_found_in_GWAS"
        qe,qo,ge,go=r["effect_allele"],r["other_allele"],r["gwas_effect_allele"],r["gwas_other_allele"]
        if pd.isna(qe) or pd.isna(qo): return "missing_QTL_alleles"
        if pd.isna(ge) or pd.isna(go): return "missing_GWAS_alleles"
        if qe==ge and qo==go: return "exact"
        if qe==go and qo==ge: return "swap_beta"
        if qe==comp.get(ge) and qo==comp.get(go): return "strand"
        if qe==comp.get(go) and qo==comp.get(ge): return "strand_and_swap_beta"
        return "incompatible"

    x["harmonization_relation"]=x.apply(relation,axis=1)
    x["beta_flip_required"]=x["harmonization_relation"].isin(["swap_beta","strand_and_swap_beta"])
    x["strand_flip_required"]=x["harmonization_relation"].isin(["strand","strand_and_swap_beta"])
    x["gwas_beta_harmonized"]=x["gwas_beta_original"]
    x.loc[x["beta_flip_required"],"gwas_beta_harmonized"]=-x.loc[x["beta_flip_required"],"gwas_beta_original"]
    x["harmonized_effect_allele"]=x["effect_allele"]; x["harmonized_other_allele"]=x["other_allele"]

    x["rejection_reason"]=""
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

    preferred=["geneid","snpid","qtl_chr","qtl_pos","effect_allele","other_allele","beta","se_qtl","p-value","gwas_chr","gwas_pos","gwas_effect_allele","gwas_other_allele","gwas_beta_original","gwas_beta_harmonized","gwas_se","gwas_p","harmonization_relation","beta_flip_required","strand_flip_required","palindromic","allele_coordinate_match","gwas_coordinate_match","harmonization_keep","rejection_reason"]
    x=x[[c for c in preferred if c in x.columns]+[c for c in x.columns if c not in preferred]]
    x.to_csv(all_out,sep="\t",index=False,compression="gzip")
    x[x["harmonization_keep"]].to_csv(kept_out,sep="\t",index=False,compression="gzip")
    x[~x["harmonization_keep"]].to_csv(rejected_out,sep="\t",index=False,compression="gzip")

    requested=len(lead_ids); allele_records=int(alleles["snpid"].nunique()); matched=int(x["found_in_gwas"].sum()); exact=int((x["harmonization_relation"]=="exact").sum())
    swapped=int((x["harmonization_relation"]=="swap_beta").sum()); strand=int((x["harmonization_relation"]=="strand").sum()); strand_swap=int((x["harmonization_relation"]=="strand_and_swap_beta").sum())
    palindromic=int((x["found_in_gwas"]&x["palindromic"]).sum()); retained=int(x["harmonization_keep"].sum()); missing_gwas=int((~x["found_in_gwas"]).sum())

    summaries.append({"tissue":tissue,"lead_eQTLs_requested":requested,"allele_records":allele_records,"matched_in_GWAS":matched,"exact":exact,"swap_beta":swapped,"strand":strand,"strand_and_swap_beta":strand_swap,"matched_palindromic":palindromic,"not_found_in_GWAS":missing_gwas,"retained_MR_ready":retained,"allele_recovery_pct":pct(allele_records,requested),"GWAS_match_pct":pct(matched,requested),"retention_pct":pct(retained,requested)})

    print(f"Lead eQTLs requested: {requested:,}")
    print(f"Allele records recovered: {allele_records:,} ({pct(allele_records,requested):.2f}%)")
    print(f"Matched in GWAS: {matched:,} ({pct(matched,requested):.2f}%)")
    print(f"Exact: {exact:,}; beta reversed: {swapped+strand_swap:,}; strand only: {strand:,}")
    print(f"Matched palindromic excluded: {palindromic:,}")
    print(f"Not found in GWAS: {missing_gwas:,}")
    print(f"MR-ready retained: {retained:,} ({pct(retained,requested):.2f}%)")
    print(f"Allele table: {allele_out}")
    print(f"MR-ready output: {kept_out}")

summary=pd.DataFrame(summaries)
summary.to_csv(summary_path,sep="\t",index=False,float_format="%.2f")
print(f"\n{'='*70}\nALL-TISSUE RETENTION SUMMARY\n{'='*70}")
print(summary.to_string(index=False,float_format=lambda x:f"{x:.2f}"))
print(f"\nSummary written to: {summary_path}")
PY

echo
echo "Finished all tissues."
column -t -s $'\t' "$SUMMARY"
