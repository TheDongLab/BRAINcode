#!/usr/bin/env python3
from pathlib import Path
import gzip, subprocess
from collections import defaultdict

B=Path.home()/"donglab/data/target_ALS/WGS_LR/repeat_comparison_gangSTR_vs_TRGT"
M=B/"inputs/master.tsv"
G=B/"gangstr/NEUAD700YFB.gangstr.vcf"
T=B/"trgt/NEUAD700YFB.trgt.vcf.gz"
O=B/"comparison"
O.mkdir(exist_ok=True)
TMP=O/"tmp"
TMP.mkdir(exist_ok=True)

def fmt(fmt,sample):
    return dict(zip(fmt.split(":"),sample.split(":")))

def called(gt):
    return gt not in ("",".","./.",".|.")

def pair(x):
    try:
        z=sorted(float(v) for v in x.split(","))
        return z if len(z)==2 else None
    except:
        return None

print("1/6 Extracting GangSTR...")
gfile=TMP/"gangstr.tsv"
with G.open() as f,gfile.open("w") as o:
    for l in f:
        if l.startswith("#"): continue
        a=l.rstrip().split("\t")
        info=dict(x.split("=",1) for x in a[7].split(";") if "=" in x)
        s=fmt(a[8],a[9])
        key=f"{a[0]}:{a[1]}:{info['END']}"
        o.write("\t".join([
            key,s.get("GT","."),s.get("REPCN","."),
            s.get("REPCI","."),s.get("DP","."),
            s.get("Q",".")
        ])+"\n")

print("2/6 Extracting TRGT...")
tfile=TMP/"trgt.tsv"
with gzip.open(T,"rt") as f,tfile.open("w") as o:
    for l in f:
        if l.startswith("#"): continue
        a=l.rstrip().split("\t")
        info=dict(x.split("=",1) for x in a[7].split(";") if "=" in x)
        s=fmt(a[8],a[9])
        key=f"{a[0]}:{int(a[1])+1}:{info['END']}"
        o.write("\t".join([
            key,s.get("GT","."),s.get("MC","."),
            s.get("AL","."),s.get("ALLR","."),
            s.get("SD","."),s.get("AP",".")
        ])+"\n")

print("3/6 Preparing master...")
mfile=TMP/"master.tsv"
with M.open() as f,mfile.open("w") as o:
    next(f)
    for l in f:
        a=l.rstrip().split("\t")
        lid,chrom,start,end,ts,te,k,motif,refbp=a
        key=f"{chrom}:{start}:{end}"
        refcopies=float(refbp)/float(k)
        o.write("\t".join([
            key,lid,chrom,start,end,motif,k,refbp,f"{refcopies:g}"
        ])+"\n")

# External sorting = low RAM
for src,name in [(mfile,"master"),(gfile,"gangstr"),(tfile,"trgt")]:
    dst=TMP/f"{name}.sorted.tsv"
    subprocess.run([
        "sort","-T",str(TMP),"-S","512M","-k1,1",
        str(src),"-o",str(dst)
    ],check=True)

ms=TMP/"master.sorted.tsv"
gs=TMP/"gangstr.sorted.tsv"
ts=TMP/"trgt.sorted.tsv"

print("4/6 Coordinate joining...")

mg=TMP/"master_gangstr.tsv"
with mg.open("w") as out:
    subprocess.run([
        "join","-t","\t","-a1","-e",".",
        "-o","1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.2,2.3,2.4,2.5,2.6",
        str(ms),str(gs)
    ],stdout=out,check=True)

allcalls=O/"all_calls.tsv"
with allcalls.open("w") as out:
    subprocess.run([
        "join","-t","\t","-a1","-e",".",
        "-o","1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,2.2,2.3,2.4,2.5,2.6,2.7",
        str(mg),str(ts)
    ],stdout=out,check=True)

print("5/6 Comparing calls...")

stats=defaultdict(int)
qsum=defaultdict(float)
psum=defaultdict(float)

discord=O/"discordant_loci.tsv"
high=O/"high_confidence_gt2.tsv"
lowg=O/"TRGT_high_GangSTR_low.tsv"
bothhigh=O/"both_high_confidence_gt2.tsv"

header=[
    "locus_id","chrom","start","end","motif","REF_copies",
    "GangSTR","GangSTR_CI","GangSTR_DP","GangSTR_Q",
    "TRGT","TRGT_bp","TRGT_range","TRGT_depth","TRGT_purity",
    "allele1_diff","allele2_diff","max_diff"
]
H="\t".join(header)+"\n"

with allcalls.open() as f, \
     discord.open("w") as d, \
     high.open("w") as h, \
     lowg.open("w") as lg, \
     bothhigh.open("w") as bh:

    d.write(H)
    h.write(H)
    lg.write(H)
    bh.write(H)

    for line in f:
        a=line.rstrip().split("\t")
        stats["MASTER"]+=1

        # master columns 1-9
        # GangSTR columns 10-14:
        # GT, REPCN, REPCI, DP, Q
        # TRGT columns 15-20:
        # GT, MC, AL, ALLR, SD, AP
        ggt=a[9]
        grepcn=a[10]
        gci=a[11]
        gdp=a[12]
        gq=a[13]

        tgt=a[14]
        tmc=a[15]
        tal=a[16]
        tallr=a[17]
        tsd=a[18]
        tap=a[19]

        gc=called(ggt)
        tc=called(tgt)

        if gc:
            stats["GangSTR_called"]+=1
        if tc:
            stats["TRGT_called"]+=1

        if gc and tc:
            stats["both_called"]+=1
        elif gc:
            stats["GangSTR_only"]+=1
        elif tc:
            stats["TRGT_only"]+=1
        else:
            stats["neither_called"]+=1

        if not(gc and tc):
            continue

        # -------------------------------------------------------------
        # AGREEMENT:
        # Only REPCN and MC are required here.
        # AP/SD must NOT affect whether a locus enters agreement counts.
        # -------------------------------------------------------------
        g=pair(grepcn)
        t=pair(tmc)

        if not(g and t):
            stats["both_called_not_comparable"]+=1
            continue

        stats["comparable_both_called"]+=1

        d1=abs(g[0]-t[0])
        d2=abs(g[1]-t[1])
        md=max(d1,d2)

        if md==0:
            group="exact"
        elif md<=1:
            group="within_1"
        elif md<=2:
            group="within_2"
        else:
            group="gt_2"

        stats[group]+=1

        # Write ALL discordant loci regardless of whether AP/SD exist.
        if md>0:
            row="\t".join([
                a[1],a[2],a[3],a[4],a[5],a[8],
                grepcn,gci,gdp,gq,
                tmc,tal,tallr,tsd,tap,
                f"{d1:g}",f"{d2:g}",f"{md:g}"
            ])+"\n"
            d.write(row)

        # -------------------------------------------------------------
        # GANGSTR CONFIDENCE:
        # Q is analyzed whenever Q is available.
        # -------------------------------------------------------------
        try:
            q=float(gq)
            stats[group+"_Q_N"]+=1
            qsum[group]+=q

            if q<0.9:
                stats[group+"_Qlt09"]+=1
            if q<0.5:
                stats[group+"_Qlt05"]+=1
        except:
            q=None

        # -------------------------------------------------------------
        # TRGT CONFIDENCE:
        # AP and SD are analyzed separately from agreement.
        # -------------------------------------------------------------
        p=pair(tap)
        sd=pair(tsd)

        if p:
            minp=min(p)
            stats[group+"_P_N"]+=1
            psum[group]+=minp

            if minp<0.9:
                stats[group+"_Plt09"]+=1
        else:
            minp=None

        if sd:
            mins=min(sd)
            stats[group+"_SD_N"]+=1

            if mins<3:
                stats[group+"_SDlt3"]+=1
        else:
            mins=None

        # -------------------------------------------------------------
        # HIGH-CONFIDENCE TRGT >2-REPEAT DISAGREEMENTS
        # Requires:
        #   difference >2
        #   BOTH TRGT alleles purity >=0.9
        #   BOTH TRGT alleles supported by >=3 HiFi reads
        # -------------------------------------------------------------
        if md>2 and minp is not None and mins is not None \
                and minp>=0.9 and mins>=3:

            stats["high_TRGT_gt2"]+=1

            row="\t".join([
                a[1],a[2],a[3],a[4],a[5],a[8],
                grepcn,gci,gdp,gq,
                tmc,tal,tallr,tsd,tap,
                f"{d1:g}",f"{d2:g}",f"{md:g}"
            ])+"\n"

            h.write(row)

            if q is None:
                stats["high_TRGT_GangSTR_Q_missing"]+=1
            elif q<0.5:
                stats["high_TRGT_GangSTR_low"]+=1
                lg.write(row)
            elif q>=0.9:
                stats["both_high_gt2"]+=1
                bh.write(row)
            else:
                stats["high_TRGT_GangSTR_mid"]+=1

print("6/6 Writing final summary...")

agreement_keys=["exact","within_1","within_2","gt_2"]
den=sum(stats[x] for x in agreement_keys)

lines=[]

def P(x=""):
    lines.append(str(x))

P("=== CATALOG / CALLABILITY ===")
P(f"MASTER              {stats['MASTER']:,}")
P(f"GangSTR called      {stats['GangSTR_called']:,} ({100*stats['GangSTR_called']/stats['MASTER']:.2f}%)")
P(f"TRGT called         {stats['TRGT_called']:,} ({100*stats['TRGT_called']/stats['MASTER']:.2f}%)")
P(f"Both called         {stats['both_called']:,} ({100*stats['both_called']/stats['MASTER']:.2f}%)")
P(f"GangSTR only        {stats['GangSTR_only']:,}")
P(f"TRGT only           {stats['TRGT_only']:,}")
P(f"Neither called      {stats['neither_called']:,}")
P(f"Comparable both     {stats['comparable_both_called']:,}")
P(f"Both not comparable {stats['both_called_not_comparable']:,}")

P()
P("=== AGREEMENT AMONG COMPARABLE BOTH-CALLED LOCI ===")
for key,label in [
    ("exact","Exact"),
    ("within_1","Differ by 1 repeat"),
    ("within_2","Differ by 2 repeats"),
    ("gt_2","Differ by >2 repeats")
]:
    P(
        f"{label:22s} "
        f"{stats[key]:9,d} "
        f"({100*stats[key]/den:.2f}%)"
    )

P()
P("=== CONFIDENCE BY AGREEMENT GROUP ===")
P(
    "Group\tN\tMean_GangSTR_Q\tGangSTR_Q<0.9\tGangSTR_Q<0.5\t"
    "Mean_min_TRGT_purity\tTRGT_purity<0.9\tEither_allele_SD<3"
)

for key,label in [
    ("exact","Exact"),
    ("within_1","Differ_by_1"),
    ("within_2","Differ_by_2"),
    ("gt_2","Differ_by_>2")
]:
    n=stats[key]

    qn=stats[key+"_Q_N"]
    pn=stats[key+"_P_N"]
    sdn=stats[key+"_SD_N"]

    meanq=qsum[key]/qn if qn else float("nan")
    meanp=psum[key]/pn if pn else float("nan")

    q09=stats[key+"_Qlt09"]
    q05=stats[key+"_Qlt05"]
    plt=stats[key+"_Plt09"]
    sdlt=stats[key+"_SDlt3"]

    P(
        f"{label}\t{n}\t"
        f"{meanq:.3f}\t"
        f"{q09} ({100*q09/qn if qn else 0:.1f}%)\t"
        f"{q05} ({100*q05/qn if qn else 0:.1f}%)\t"
        f"{meanp:.3f}\t"
        f"{plt} ({100*plt/pn if pn else 0:.1f}%)\t"
        f"{sdlt} ({100*sdlt/sdn if sdn else 0:.1f}%)"
    )

P()
P("=== HIGH-CONFIDENCE TRGT >2-REPEAT DISAGREEMENTS ===")
P("Definition: min TRGT purity >=0.9 AND >=3 HiFi reads supporting each allele")

n=stats["high_TRGT_gt2"]

P(f"Total                         {n:,}")

if n:
    P(
        f"GangSTR Q <0.5                "
        f"{stats['high_TRGT_GangSTR_low']:,} "
        f"({100*stats['high_TRGT_GangSTR_low']/n:.1f}%)"
    )
    P(
        f"GangSTR 0.5<=Q<0.9            "
        f"{stats['high_TRGT_GangSTR_mid']:,} "
        f"({100*stats['high_TRGT_GangSTR_mid']/n:.1f}%)"
    )
    P(
        f"GangSTR Q >=0.9               "
        f"{stats['both_high_gt2']:,} "
        f"({100*stats['both_high_gt2']/n:.1f}%)"
    )
    if stats["high_TRGT_GangSTR_Q_missing"]:
        P(
            f"GangSTR Q missing              "
            f"{stats['high_TRGT_GangSTR_Q_missing']:,} "
            f"({100*stats['high_TRGT_GangSTR_Q_missing']/n:.1f}%)"
        )

P()
P("=== INTERNAL CONSISTENCY CHECKS ===")
P(
    f"Callability total: "
    f"{stats['both_called']+stats['GangSTR_only']+stats['TRGT_only']+stats['neither_called']:,} "
    f"(expected {stats['MASTER']:,})"
)
P(
    f"Agreement total:   {den:,} "
    f"(expected comparable both-called {stats['comparable_both_called']:,})"
)

text="\n".join(lines)+"\n"
(O/"final_summary.txt").write_text(text)

print(text,end="")

print("Outputs:")
for x in [
    allcalls,
    discord,
    high,
    lowg,
    bothhigh,
    O/"final_summary.txt"
]:
    print(" ",x)
