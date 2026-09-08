#!/bin/bash
#SBATCH --job-name=ALS_variant_evidence
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/ALS_variant_evidence.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/ALS_variant_evidence.err
#SBATCH --time=06:00:00
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

set -euo pipefail

ROOT="$HOME/donglab/data/target_ALS"
GWAS="$HOME/donglab/references/GWAS/ALS/harmonised/34873335-GCST90027163-MONDO_0004976.h.tsv.gz"
OUTDIR="$ROOT/QTL/ALS_variant_validation_evidence"

mkdir -p "$OUTDIR"

python -u - "$GWAS" "$ROOT" "$OUTDIR" <<'PY'
import sys,re,gzip
from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np

GWAS=Path(sys.argv[1])
ROOT=Path(sys.argv[2])
OUTDIR=Path(sys.argv[3])

TISSUES=[
    "Cerebellum",
    "Frontal_Cortex",
    "Cervical_Spinal_Cord",
    "Lumbar_Spinal_Cord",
    "Motor_Cortex"
]

QTLS=["eQTL","sQTL","cQTL"]

def norm(x):
    return re.sub(r"[^a-z0-9]","",str(x).lower())

def find_col(cols, aliases, required=False):
    lookup={norm(c):c for c in cols}
    for a in aliases:
        if norm(a) in lookup:
            return lookup[norm(a)]
    for c in cols:
        nc=norm(c)
        for a in aliases:
            if norm(a) and norm(a) in nc:
                return c
    if required:
        raise RuntimeError(f"Could not find {aliases}; columns={list(cols)}")
    return None

def read_table(path):
    try:
        return pd.read_csv(path,sep="\t",dtype=str,low_memory=False)
    except:
        return pd.read_csv(path,sep=r"\s+",dtype=str,engine="python")

def clean_snp(x):
    if x is None:return None
    x=str(x).strip()
    if not x or x.lower() in {"nan","na","none","."}:return None
    return x

# ------------------------------------------------------------------
# 1. LOAD ALS GWAS
# ------------------------------------------------------------------
print(f"Loading GWAS: {GWAS}",flush=True)
gwas=read_table(GWAS)

snp_col=find_col(
    gwas.columns,
    ["rsid","rs_id","variant_id","variant","snp","snpid"],
    required=True
)
p_col=find_col(gwas.columns,["p_value","pvalue","pval","p"])
beta_col=find_col(gwas.columns,["beta","effect","effect_size"])
se_col=find_col(gwas.columns,["standard_error","stderr","se"])
ea_col=find_col(gwas.columns,["effect_allele","effectallele","ea"])
oa_col=find_col(gwas.columns,["other_allele","otherallele","nea","non_effect_allele"])
chr_col=find_col(gwas.columns,["chromosome","chrom","chr"])
pos_col=find_col(gwas.columns,["base_pair_location","position","pos","bp"])

gwas["SNP"]=gwas[snp_col].map(clean_snp)
gwas=gwas[gwas["SNP"].notna()].drop_duplicates("SNP").copy()

out=pd.DataFrame({"SNP":gwas["SNP"]})

def copy_gwas(col,new):
    if col:
        out[new]=gwas[col].values
    else:
        out[new]=np.nan

copy_gwas(chr_col,"GWAS_chr")
copy_gwas(pos_col,"GWAS_pos")
copy_gwas(ea_col,"GWAS_effect_allele")
copy_gwas(oa_col,"GWAS_other_allele")
copy_gwas(beta_col,"GWAS_beta")
copy_gwas(se_col,"GWAS_se")
copy_gwas(p_col,"GWAS_p")

out=out.set_index("SNP",drop=False)

# evidence accumulators
genes=defaultdict(set)

for q in QTLS:
    out[f"{q}_significant"]=False
    out[f"{q}_tissues"]=""
    out[f"{q}_traits"]=""
    out[f"{q}_best_FDR"]=np.nan

    out[f"{q}_SMR_significant"]=False
    out[f"{q}_SMR_tissues"]=""
    out[f"{q}_best_p_SMR"]=np.nan
    out[f"{q}_HEIDI_pass"]=False

    out[f"{q}_coloc_significant"]=False
    out[f"{q}_coloc_tissues"]=""
    out[f"{q}_best_PPH4"]=np.nan

    out[f"{q}_SuSiE_credible_set"]=False
    out[f"{q}_best_PIP"]=np.nan

# temporary sets
qtl_tissues={q:defaultdict(set) for q in QTLS}
qtl_traits={q:defaultdict(set) for q in QTLS}
smr_tissues={q:defaultdict(set) for q in QTLS}
coloc_tissues={q:defaultdict(set) for q in QTLS}

# ------------------------------------------------------------------
# 2. SIGNIFICANT QTL RESULTS
# ------------------------------------------------------------------
for tissue in TISSUES:
    for q in QTLS:
        f=ROOT/tissue/q/"results"/f"{tissue}_{q}.FDR0.05.txt"
        if not f.exists():
            print(f"MISSING QTL: {f}",flush=True)
            continue

        print(f"Reading {q}: {tissue}",flush=True)
        d=read_table(f)

        sc=find_col(d.columns,["snpid","snp_id","snp","rsid","variant"],True)
        fc=find_col(d.columns,["fdr","qvalue","q_value","padj"])
        tc=find_col(
            d.columns,
            ["geneid","gene_id","junction_id","circ_id","phenotype_id","probeid","event_id"]
        )
        gc=find_col(d.columns,["gene_symbol","gene_name","symbol"])

        for _,r in d.iterrows():
            snp=clean_snp(r.get(sc))
            if snp not in out.index:
                continue

            out.at[snp,f"{q}_significant"]=True
            qtl_tissues[q][snp].add(tissue)

            if tc and pd.notna(r.get(tc)):
                qtl_traits[q][snp].add(str(r[tc]))

            if gc and pd.notna(r.get(gc)):
                genes[snp].add(str(r[gc]))

            if q=="eQTL" and tc and pd.notna(r.get(tc)):
                genes[snp].add(str(r[tc]))

            if fc:
                try:
                    x=float(r[fc])
                    old=out.at[snp,f"{q}_best_FDR"]
                    if pd.isna(old) or x<old:
                        out.at[snp,f"{q}_best_FDR"]=x
                except:
                    pass

# ------------------------------------------------------------------
# 3. SMR + HEIDI
# ------------------------------------------------------------------
for tissue in TISSUES:
    for q in QTLS:
        f=ROOT/tissue/"MR"/f"{q}_SMR_HEIDI"/"results"/f"{tissue}_{q}_SMR_HEIDI.smr"
        if not f.exists():
            print(f"MISSING SMR: {f}",flush=True)
            continue

        print(f"Reading SMR: {q} {tissue}",flush=True)
        d=read_table(f)

        sc=find_col(d.columns,["topSNP","top_snp","snp","rsid"])
        pc=find_col(d.columns,["p_SMR","psmr"])
        hc=find_col(d.columns,["p_HEIDI","pheidi"])
        gc=find_col(d.columns,["Gene","gene","gene_symbol"])

        if sc is None:
            continue

        # Recalculate tissue BH FDR directly from current .smr results
        if pc:
            p=pd.to_numeric(d[pc],errors="coerce")
            valid=p.notna()
            ranked=p[valid].rank(method="min")
            m=valid.sum()

            # standard BH monotonic correction
            vals=p[valid].sort_values()
            adj=(vals*m/np.arange(1,m+1)).clip(upper=1)
            adj=np.minimum.accumulate(adj[::-1])[::-1]
            fdr=pd.Series(np.nan,index=d.index)
            fdr.loc[vals.index]=adj
        else:
            fdr=pd.Series(np.nan,index=d.index)

        for idx,r in d.iterrows():
            snp=clean_snp(r.get(sc))
            if snp not in out.index:
                continue

            if gc and pd.notna(r.get(gc)):
                genes[snp].add(str(r[gc]))

            try:
                psmr=float(r[pc]) if pc else np.nan
                old=out.at[snp,f"{q}_best_p_SMR"]
                if np.isfinite(psmr) and (pd.isna(old) or psmr<old):
                    out.at[snp,f"{q}_best_p_SMR"]=psmr
            except:
                pass

            # tissue-level FDR < 0.05
            try:
                if float(fdr.loc[idx])<0.05:
                    out.at[snp,f"{q}_SMR_significant"]=True
                    smr_tissues[q][snp].add(tissue)
            except:
                pass

            # HEIDI pass = p >= 0.05
            if hc:
                try:
                    if float(r[hc])>=0.05:
                        out.at[snp,f"{q}_HEIDI_pass"]=True
                except:
                    pass

# ------------------------------------------------------------------
# 4. SuSiE / COLOC
#
# We search recursively because your three scripts write several
# per-locus/global files under each *_SuSiE_coloc directory.
# ------------------------------------------------------------------
for tissue in TISSUES:
    for q in QTLS:
        droot=ROOT/tissue/"MR"/f"{q}_SuSiE_coloc"
        if not droot.exists():
            continue

        for f in droot.rglob("*.tsv*"):
            try:
                d=read_table(f)
            except:
                continue

            sc=find_col(d.columns,["candidate_snp","snp","snpid","rsid","variant"])
            if sc is None:
                continue

            pph4c=find_col(d.columns,["PPH4","PP.H4","H4"])
            h4pass=find_col(d.columns,["H4_reference_pass","H4_pass","coloc_pass"])
            pipc=find_col(d.columns,["PIP","susie_pip","posterior_inclusion_probability"])
            csc=find_col(d.columns,["credible_set","in_credible_set","cs"])

            # Ignore unrelated tables that have SNP but none of the coloc/SuSiE evidence
            if not any([pph4c,h4pass,pipc,csc]):
                continue

            for _,r in d.iterrows():
                snp=clean_snp(r.get(sc))
                if snp not in out.index:
                    continue

                if pph4c:
                    try:
                        x=float(r[pph4c])
                        old=out.at[snp,f"{q}_best_PPH4"]
                        if pd.isna(old) or x>old:
                            out.at[snp,f"{q}_best_PPH4"]=x
                    except:
                        pass

                passed=False
                if h4pass:
                    v=str(r[h4pass]).strip().lower()
                    passed=v in {"true","t","1","yes","pass"}

                # If no explicit pass column exists, use PPH4 >= 0.80
                # only as a fallback.
                elif pph4c:
                    try:
                        passed=float(r[pph4c])>=0.80
                    except:
                        pass

                if passed:
                    out.at[snp,f"{q}_coloc_significant"]=True
                    coloc_tissues[q][snp].add(tissue)

                if pipc:
                    try:
                        x=float(r[pipc])
                        old=out.at[snp,f"{q}_best_PIP"]
                        if pd.isna(old) or x>old:
                            out.at[snp,f"{q}_best_PIP"]=x
                    except:
                        pass

                if csc:
                    v=str(r[csc]).strip().lower()
                    if v not in {"","nan","na","none","false","0","no"}:
                        out.at[snp,f"{q}_SuSiE_credible_set"]=True

# ------------------------------------------------------------------
# 5. COLLAPSE SETS
# ------------------------------------------------------------------
for q in QTLS:
    for snp,v in qtl_tissues[q].items():
        out.at[snp,f"{q}_tissues"]=";".join(sorted(v))
    for snp,v in qtl_traits[q].items():
        out.at[snp,f"{q}_traits"]=";".join(sorted(v))
    for snp,v in smr_tissues[q].items():
        out.at[snp,f"{q}_SMR_tissues"]=";".join(sorted(v))
    for snp,v in coloc_tissues[q].items():
        out.at[snp,f"{q}_coloc_tissues"]=";".join(sorted(v))

out["candidate_genes"]=[
    ";".join(sorted(genes.get(s,set())))
    for s in out.index
]

# ------------------------------------------------------------------
# 6. CURRENTLY UNAVAILABLE ANNOTATIONS
# ------------------------------------------------------------------
out["TF_binding_disruption"]="NA"
out["regulatory_annotation"]="NA"
out["coding_consequence"]="NA"

# ------------------------------------------------------------------
# 7. SIMPLE EVIDENCE COUNTS
#
# Keep the individual columns as the scientifically useful output.
# Score is only a prioritization convenience.
# ------------------------------------------------------------------
score_cols=[
    "eQTL_significant",
    "sQTL_significant",
    "cQTL_significant",
    "eQTL_SMR_significant",
    "sQTL_SMR_significant",
    "cQTL_SMR_significant",
    "eQTL_HEIDI_pass",
    "sQTL_HEIDI_pass",
    "cQTL_HEIDI_pass",
    "eQTL_coloc_significant",
    "sQTL_coloc_significant",
    "cQTL_coloc_significant",
    "eQTL_SuSiE_credible_set",
    "sQTL_SuSiE_credible_set",
    "cQTL_SuSiE_credible_set"
]

out["evidence_count"]=out[score_cols].astype(int).sum(axis=1)

out["any_QTL"]=(
    out["eQTL_significant"] |
    out["sQTL_significant"] |
    out["cQTL_significant"]
)

out["any_SMR"]=(
    out["eQTL_SMR_significant"] |
    out["sQTL_SMR_significant"] |
    out["cQTL_SMR_significant"]
)

out["any_coloc"]=(
    out["eQTL_coloc_significant"] |
    out["sQTL_coloc_significant"] |
    out["cQTL_coloc_significant"]
)

out["any_evidence"]=out["evidence_count"]>0

# priority is intentionally transparent, not a biological claim
out["validation_priority"]=np.select(
    [
        (out["any_coloc"] & out["any_SMR"] & out["any_QTL"]),
        (out["any_coloc"] & out["any_QTL"]),
        (out["any_SMR"] & out["any_QTL"]),
        (out["any_QTL"])
    ],
    ["HIGH","HIGH","MEDIUM","LOW"],
    default="NONE"
)

# ------------------------------------------------------------------
# 8. OUTPUT
# ------------------------------------------------------------------
out=out.reset_index(drop=True)

allfile=OUTDIR/"ALS_GWAS_variant_validation_evidence.ALL.tsv.gz"
evfile=OUTDIR/"ALS_GWAS_variant_validation_evidence.WITH_EVIDENCE.tsv"
rankfile=OUTDIR/"ALS_GWAS_variant_validation_evidence.RANKED.tsv"

out.to_csv(allfile,sep="\t",index=False,compression="gzip")

evidence=out[out["any_evidence"]].copy()
evidence.to_csv(evfile,sep="\t",index=False)

ranked=evidence.sort_values(
    ["evidence_count","GWAS_p"],
    ascending=[False,True],
    na_position="last"
)
ranked.to_csv(rankfile,sep="\t",index=False)

print()
print("="*72)
print("DONE")
print(f"GWAS SNPs             : {len(out):,}")
print(f"SNPs with QTL evidence: {out['any_QTL'].sum():,}")
print(f"SNPs with SMR evidence: {out['any_SMR'].sum():,}")
print(f"SNPs with coloc       : {out['any_coloc'].sum():,}")
print(f"SNPs with any evidence: {out['any_evidence'].sum():,}")
print()
print(f"ALL    : {allfile}")
print(f"EVIDENCE: {evfile}")
print(f"RANKED : {rankfile}")
print("="*72)
PY
