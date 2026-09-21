#!/usr/bin/env bash
#SBATCH --job-name=analyze_STR_caller_results
#SBATCH --output=/home/zw529/donglab/data/target_ALS/WGS_LR/analyze_STR_caller_results.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/WGS_LR/analyze_STR_caller_results.err
#SBATCH --time=04:00:00
#SBATCH --mem=56G
#SBATCH --cpus-per-task=4

# Existing coverage tables: bash analyze_STR_caller_results.sh --read-support-only
# Existing evidence and tables: bash analyze_STR_caller_results.sh --plots-only
# Supply --trgt-spanning-bam /absolute/path/file.bam if auto-detection is ambiguous.

set -euo pipefail
module load SAMtools
exec "${PYTHON:-$HOME/donglab/pipelines/modules/miniconda3/bin/python}" - "$@" <<'PY'
import argparse, ast, csv, functools, gzip, hashlib, json, math, re, shlex, subprocess, sys, textwrap
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import PercentFormatter
def render_plots(output, include_depth=True):
 D=Path(output)
 plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":180})
 A={"exact_rate":"Exact match: both alleles","le1_rate":"Within 1 repeat on both alleles","le2_rate":"Within 2 repeats on both alleles","gt2_rate":"More than 2 repeats on either allele"}
 Q={"gang_q_lt0.9_rate":"GangSTR confidence Q < 0.9","gang_q_lt0.5_rate":"GangSTR confidence Q < 0.5","trgt_min_ap_lt0.9_rate":"TRGT minimum purity AP < 0.9","trgt_min_ap_lt0.5_rate":"TRGT minimum purity AP < 0.5"}
 C=["#0072B2","#E69F00","#009E73","#CC79A7"]
 NOTE="STR = short tandem repeat; locus = one genomic location. GangSTR uses short reads; TRGT uses long reads."
 AG="Agreement compares the two allele repeat counts after sorting them by size; agreement does not establish accuracy."
 QUAL="Q measures genotype confidence; AP measures how closely an allele follows its repeat motif. Minimum AP is the lower purity of the two alleles. These are different measures. Quality percentages use available values for each caller, so their denominators can differ."
 def read(n): return pd.read_csv(D/("analysis_"+n+".tsv"),sep="\t")
 def frame(n,title,sub,w=14):
  f,ax=plt.subplots(n,1,figsize=(w,3.25*n+3.4),squeeze=False)
  f.suptitle(title,fontsize=19,fontweight="bold",y=.98)
  f.text(.08,.93,sub,fontsize=11,va="top")
  return f,ax[:,0]
 def done(f,name,note):
  f.subplots_adjust(top=.85,bottom=.20,hspace=.55,left=.09,right=.98)
  f.text(.08,.025,textwrap.fill(NOTE+" "+note,145),fontsize=10,va="bottom",linespacing=1.45)
  f.savefig(D/("analysis_"+name+".png")); plt.close(f)
 def lines(ax,z,x,items,percent=True,colors=None,solid=False):
  for old_line in list(ax.lines):
   old_line.remove()
  old_legend=ax.get_legend()
  if old_legend is not None:
   old_legend.remove()
  handles=[]
  for i,(col,label) in enumerate(items.items()):
   line,=ax.plot(x,z[col],marker="o",ms=4,color=(C if colors is None else colors)[i],ls="-" if solid else ("--" if i%2 else "-"),label=label)
   handles.append(line)
  if percent:
   ax.yaxis.set_major_formatter(PercentFormatter(1))
   ax.set_ylim(bottom=0)
  ax.grid(axis="y",alpha=.2)
  ax.legend(handles=handles,labels=list(items.values()),fontsize=9,ncol=2,loc="lower left",bbox_to_anchor=(0,1.01))
 for key,label,definition in [("gang_dp","GangSTR informative-read count (DP)","DP counts reads used as evidence by GangSTR; it is not mean BAM coverage."),("trgt_min_sd","TRGT read support for the less-supported allele (minimum SD)","SD is supporting reads per allele; the smaller of the two values defines each bin."),("short_depth","Short-read mean coverage across each STR (x)","Coverage is the average aligned-read depth over the full reference STR interval, including uncovered bases."),("long_depth","Long-read mean coverage across each STR (x)","Coverage is the average aligned-read depth over the full reference STR interval, including uncovered bases.")]:
  if not include_depth and key in ("short_depth","long_depth"): continue
  z=read("by_"+key); x=np.arange(len(z)); bins=z[key+"_bin"].astype(str).str.replace("-<"," to <",regex=False).str.replace(">=","at least ",regex=False)
  f,ax=frame(2,"Repeat-call agreement and quality versus read evidence",label+" | NEUAD700YFB | N beneath each bin = comparable loci")
  lines(ax[0],z,x,{k:A[k] for k in ["exact_rate","le1_rate","gt2_rate"]},colors=["#78B4DC","#8A2BE2","#8B0000"],solid=True); ax[0].set_ylabel("Comparable loci (%)"); ax[0].set_ylim(0,1.15)
  lines(ax[1],z,x,Q); ax[1].set_ylabel("Calls below threshold (%)")
  for b in ax: b.set_xticks(x); b.set_xticklabels([v+"\nN="+format(int(n),",") for v,n in zip(bins,z.N_comparable)],fontsize=9)
  ax[1].set_xlabel(label)
  done(f,"by_"+key,definition+" "+AG+" Within 1 includes exact matches. "+QUAL+" Small-N bins are less stable; these are observational associations.")
 z=read("reference_length_10bp").sort_values("length_bin_start"); z=z.set_index("length_bin_start").reindex(range(0,int(z.length_bin_start.max())+10,10)).rename_axis("length_bin_start").reset_index(); x=z.length_bin_start+5
 f,ax=frame(2,"How does repeat length relate to caller agreement?","NEUAD700YFB | Reference STR length in 10-base-pair bins | All observed lengths shown")
 lines(ax[0],z,x,{k:A[k] for k in ["exact_rate","le1_rate","le2_rate"]}); ax[0].set_ylabel("Comparable loci (%)"); ax[0].set_ylim(0,1.15)
 for col,color in zip(["exact_rate","le1_rate","le2_rate"],C):
  low=z.N_comparable<100; ax[0].scatter(x[low],z.loc[low,col],s=55,facecolors="white",edgecolors=color,zorder=5)
 bottom=np.zeros(len(z))
 for col,label,color in zip(["exact","0_to_1","1_to_2","gt2"],["Exact match","Difference >0 to 1 repeat","Difference >1 to 2 repeats","Difference >2 repeats"],C):
  y=z[col+"_exclusive_N"].fillna(0); ax[1].bar(z.length_bin_start,y,width=10,align="edge",bottom=bottom,label=label,color=color); bottom+=y
 ax[1].set_yscale("symlog",linthresh=1); ax[1].set_ylabel("Number of comparable loci\n(log scale above 1)"); ax[1].legend(fontsize=9,ncol=2); ax[1].grid(axis="y",alpha=.2)
 for b in ax: b.set_xlabel("Reference STR length (base pairs)")
 done(f,"reference_length_10bp",AG+" Bins are 0-9, 10-19, etc.; points sit at bin centers. Hollow points mark fewer than 100 comparable loci. Top curves are cumulative; bottom categories are exclusive. Missing bins have no connecting line.")
 z=read("motif_length"); x=np.arange(len(z))
 f,ax=frame(3,"Does the size of the repeating unit relate to caller results?","NEUAD700YFB | Motif length is the number of DNA bases in one repeating unit (e.g., CAG = 3 bases).")
 ax[0].bar(x,z.exact_rate,color=C[0]); ax[0].set_ylim(0,1.35); ax[0].set_yticks(np.linspace(0,1,6)); ax[0].yaxis.set_major_formatter(PercentFormatter(1)); ax[0].set_ylabel("Exact agreement (%)")
 for i,row in enumerate(z.itertuples()): ax[0].text(i,row.exact_rate+.025,"N="+format(int(row.N_comparable),","),ha="center",rotation=90 if len(z)>10 else 0,fontsize=8)
 lines(ax[1],z,x,{"reference_bp_median":"Reference repeat length","gang_max_bp_median":"GangSTR: longer allele at each locus","trgt_max_allele_bp_median":"TRGT: longer allele at each locus"},False); ax[1].set_ylabel("Median length (base pairs)")
 lines(ax[2],z,x,Q); ax[2].set_ylabel("Calls below threshold (%)")
 for b in ax: b.set_xticks(x); b.set_xticklabels(z.motif_length_group,rotation=90 if len(z)>10 else 0)
 ax[2].set_xlabel("Repeating-unit length (base pairs)")
 done(f,"motif_length",AG+" N counts comparable loci. Length medians use available calls: GangSTR copies x motif length; TRGT allele length. "+QUAL+" Motif groups may differ in length and depth; this plot does not adjust for those differences.")
 raw=pd.read_csv(D/"analysis_loci.tsv.gz",sep="\t",usecols=["motif","max_diff","reference_bp","gang_max_bp","trgt_max_allele_bp","gang_q","trgt_min_ap"])
 rows=[]
 for motif,g in raw.groupby("motif",sort=True):
  v=g.max_diff.dropna()
  row={"motif":motif,"motif_length":len(motif),"N_master":len(g),"N_comparable":len(v),"exact_rate":(v==0).mean() if len(v) else np.nan}
  for c in ["reference_bp","gang_max_bp","trgt_max_allele_bp"]:
   values=g[c].dropna(); row[c+"_median"]=values.median() if len(values) else np.nan
  for c in ["gang_q","trgt_min_ap"]:
   values=g[c].dropna(); row[c+"_N"]=len(values)
   for cutoff in [.9,.5]:
    label=c+"_lt"+str(cutoff)
    row[label+"_N"]=int((values<cutoff).sum())
    row[label+"_rate"]=(values<cutoff).mean() if len(values) else np.nan
  rows.append(row)
 exact=pd.DataFrame(rows).sort_values(["N_comparable","motif"],ascending=[False,True])
 exact.to_csv(D/"analysis_exact_motifs_all.tsv",sep="\t",index=False,na_rep="NA")
 exact.loc[exact.N_comparable>=100].to_csv(D/"analysis_exact_motifs_Nge100.tsv",sep="\t",index=False,na_rep="NA")
 for table_name,key,plot_name,title,unit,xlabel,group_note in [
  ("motif_families_Nge100","motif_family","motif_families_top30","Repeat-sequence families: agreement, allele length and quality","families","Motif family (representative DNA sequence)","Families combine rotations and reverse complements."),
  ("exact_motifs_Nge100","motif","exact_motifs_top30","Original motif sequences: agreement, allele length and quality","original motifs","Original catalog motif (exact nucleotide order)","Motifs are grouped by the original catalog string: no rotations or reverse complements are combined. Exact motif grouping is separate from exact agreement of allele repeat counts.")
 ]:
  z=read(table_name); z=z[z.N_comparable>=100].sort_values(["N_comparable",key],ascending=[False,True]).head(30); x=np.arange(len(z))
  f,a=plt.subplots(3,1,figsize=(20,15)); f.suptitle(title,fontsize=20,fontweight="bold",y=.98)
  f.text(.08,.935,"NEUAD700YFB | Top 30 "+unit+" by comparable-locus count | Each group has at least 100 comparable loci")
  a[0].bar(x,z.exact_rate,color="#666666",label="GangSTR–TRGT exact agreement"); a[0].set_ylim(0,1.35); a[0].set_yticks(np.linspace(0,1,6)); a[0].set_ylabel("Exact agreement (%)"); a[0].yaxis.set_major_formatter(PercentFormatter(1))
  for i,r in enumerate(z.itertuples()): a[0].text(i,r.exact_rate+.02,f"N={int(r.N_comparable):,}",ha="center",rotation=90,fontsize=8)
  for j,(c,label,color) in enumerate([("reference_bp_median","Reference STR length","#999999"),("gang_max_bp_median","GangSTR: longer allele","#0072B2"),("trgt_max_allele_bp_median","TRGT: longer allele","#AE3038")]): a[1].bar(x+(j-1)*.25,z[c],width=.25,color=color,label=label)
  a[1].set_ylabel("Typical repeat length (base pairs)\nMedian across locations in each group")
  for j,(c,label,color,hatch) in enumerate([("gang_q_lt0.9_rate","GangSTR Q < 0.9","#78B4DC",""),("gang_q_lt0.5_rate","GangSTR Q < 0.5","#175A91",""),("trgt_min_ap_lt0.9_rate","TRGT minimum AP < 0.9","#EF9898",""),("trgt_min_ap_lt0.5_rate","TRGT minimum AP < 0.5","#AE3038","")]): a[2].bar(x+(j-1.5)*.2,z[c],width=.2,color=color,hatch=hatch,edgecolor="black",linewidth=.3,label=label)
  a[2].set_ylabel("Calls below threshold (%)"); a[2].yaxis.set_major_formatter(PercentFormatter(1))
  for ax in a: ax.set_xticks(x); ax.set_xticklabels(z[key],rotation=90); ax.grid(axis="y",alpha=.15); ax.set_axisbelow(True); ax.legend(loc="lower left",bbox_to_anchor=(0,1.01),ncol=4,fontsize=10)
  a[2].set_xlabel(xlabel)
  note="STR = short tandem repeat. GangSTR uses short reads; TRGT uses long reads. N counts comparable loci. Exact agreement requires matching both sorted allele repeat counts. Middle panel: at each location, take the longer of the two alleles called in this sample, then show the median of those lengths across the group. Gray is the median reference-genome repeat length; blue is GangSTR; red is TRGT. Equal medians can hide disagreements at individual locations. Q measures genotype confidence; AP measures repeat purity, using the less-pure allele. These scores are not equivalent. Bottom panel: light/dark blue = GangSTR; light/dark red = TRGT. Dark bars show values <0.5, a subset of the light bars showing values <0.9. Percentages use available values for each caller; denominators can differ. Agreement is not proof of accuracy; length and depth differences between groups are not adjusted for."
  note=group_note+" "+note
  if z.empty: a[0].text(.5,.5,"No groups meet the 100-comparable-locus threshold",transform=a[0].transAxes,ha="center")
  f.text(.08,.025,textwrap.fill(note,175),fontsize=10); f.subplots_adjust(top=.87,bottom=.22,hspace=.65,left=.08,right=.98); f.savefig(D/("analysis_"+plot_name+".png"),dpi=180); plt.close(f)
 if include_depth:
  z=pd.read_csv(D/"analysis_loci.tsv.gz",sep="\t",usecols=["short_depth","long_depth"])
  z=z[np.isfinite(z).all(axis=1)&(z>=0).all(axis=1)]
  x,y=np.log1p(z.short_depth.to_numpy()),np.log1p(z.long_depth.to_numpy())
  assert len(x)>2 and np.ptp(x)>0 and np.ptp(y)>0,"Insufficient variation for regression"
  m,b=np.polyfit(x,y,1); r2=1-np.sum((y-(m*x+b))**2)/np.sum((y-y.mean())**2)
  f,a=plt.subplots(figsize=(14,10)); f.suptitle("Short- and long-read coverage at the same STR loci",fontsize=19,fontweight="bold",y=.98)
  f.text(.09,.92,f"NEUAD700YFB | {len(z):,} loci | Median coverage: short reads {z.short_depth.median():.1f}x; long reads {z.long_depth.median():.1f}x")
  h=a.hexbin(x,y,gridsize=70,bins="log",mincnt=1,cmap="viridis"); lim=max(x.max(),y.max()); xx=np.array([x.min(),x.max()])
  a.plot([0,lim],[0,lim],"k--",label="Equal coverage")
  a.plot(xx,m*xx+b,color="red",ls=":",lw=2.5,label="Linear regression on log(1 + depth)")
  a.text(.98,.02,f"$R^2$ = {r2:.4f}\nFit in log(1 + depth) coordinates",transform=a.transAxes,ha="right",va="bottom",color="red",bbox=dict(facecolor="white",alpha=.9,edgecolor="none"))
  t=np.array([0,1,5,10,20,50,100,1000,10000,100000]); t=t[np.log1p(t)<=lim]
  for axis in [a.xaxis,a.yaxis]: axis.set_ticks(np.log1p(t)); axis.set_ticklabels([f"{v:,}" for v in t])
  a.set(xlabel="Short-read mean depth across each STR (x)",ylabel="Long-read mean depth across each STR (x)"); a.legend(loc="upper left",fontsize=9); f.colorbar(h,ax=a,label="Loci per hexagon (log color scale)")
  note="STR = short tandem repeat. Brighter hexagons contain more loci. Above the black line, long-read coverage is greater. Axes use natural log(1 + depth) spacing, labeled in actual coverage. Red: ordinary least-squares fit to individual loci, including zero depths; R-squared is explained variation in transformed long-read depth. Coverage averages the full reference STR interval, including uncovered bases. Association does not establish equal coverage or caller accuracy."
  f.text(.09,.03,textwrap.fill(note,140),fontsize=10); f.subplots_adjust(top=.85,bottom=.24,right=.96); f.savefig(D/"analysis_paired_BAM_depth.png",dpi=180); plt.close(f)
 print("Regenerated "+str(9 if include_depth else 6)+" labeled plots in "+str(D)+"; no BAMs reread.")
def generate_read_supported_evidence(base, output, bam, samtools):
 BINS = [0,50,100,150,200,300,float('inf')]
 LABELS = ['<50','50–99','100–149','150–199','200–299','≥300']
 CATEGORIES = ['Exact genotype match','Within 1 repeat on both alleles','Longer allele underestimated','Longer allele overestimated','Other genotype disagreement','GangSTR no call','GangSTR count unavailable']
 def pair(s):
  try:
   v = tuple(float(x) for x in s.split(','))
   return v if len(v)==2 and all(math.isfinite(x) and x>=0 for x in v) else None
  except (ValueError, AttributeError):
   return None
 def called(s):
  return bool(re.fullmatch(r'\d+[|/]\d+',s))
 def classify(g, t, gt):
  if not called(gt): return 'GangSTR no call'
  if g is None: return 'GangSTR count unavailable'
  g,t = sorted(g),sorted(t)
  if g==t: return 'Exact genotype match'
  if max(abs(g[i]-t[i]) for i in (0,1))<=1: return 'Within 1 repeat on both alleles'
  if g[1]<t[1]: return 'Longer allele underestimated'
  if g[1]>t[1]: return 'Longer allele overestimated'
  return 'Other genotype disagreement'
 def read_length(fields, locus, args):
  flag,pos,mapq = int(fields[1]),int(fields[3])-1,int(fields[4])
  if flag & (4|256|512|1024|2048) or mapq==255 or mapq<args.mapq: return None
  if fields[2]!=locus['chrom'] or fields[9]=='*': return None
  ops = [(int(n),op) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',fields[5])]
  if not ops or any(op in 'NSHP' for n,op in ops): return None
  end = pos + sum(n for n,op in ops if op in 'MD=X')
  if pos>locus['start0']-args.flank or end<locus['end0']+args.flank: return None
  tags = dict((x[:2],x[5:]) for x in fields[11:] if len(x)>5 and x[2]==':')
  try:
   fl = [int(x) for x in tags['FL'].split(',')[1:]]
   if len(fl)!=2 or min(fl)<args.flank: return None
   if float(tags['rq'])<args.rq or not math.isfinite(float(tags['rq'])): return None
   if int(tags['SO'])>-args.flank or int(tags['EO'])<args.flank: return None
  except (KeyError,ValueError): return None
  if sum(n for n,op in ops if op in 'MI=X')!=len(fields[9]): return None
  size = len(fields[9])-sum(fl)
  return size if size>0 else None
 def support(lengths, locus, args):
  targets = sorted(locus['al']); tol = max(2, locus['motif_length'])
  counts=[0,0]; good=0
  for length in lengths:
   delta=[abs(length-t) for t in targets]
   if min(delta)>tol: continue
   good+=1
   if targets[0]==targets[1]: counts[0]+=1
   elif delta[0]<delta[1]: counts[0]+=1
   elif delta[1]<delta[0]: counts[1]+=1
  enough = counts[0]>=2*args.per_allele if targets[0]==targets[1] else min(counts)>=args.per_allele
  rate=good/len(lengths) if lengths else 0
  return enough and rate>=args.consistency,counts,rate
 args=SimpleNamespace(base=Path(base),mapq=20,rq=.99,flank=50,per_allele=3,consistency=.8)
 comp=Path(output);out=comp/'read_supported_STR_benchmark';out.mkdir(parents=True,exist_ok=True)
 (out/'evidence_SUCCESS').unlink(missing_ok=True)
 subprocess.run([samtools,'quickcheck','-v',str(bam)],check=True)
 inputs=[Path(bam),args.base/'inputs/master.tsv',args.base/'comparison/all_calls.tsv',comp/'analysis_loci.tsv.gz']
 before={str(f.resolve()): (f.stat().st_size,f.stat().st_mtime_ns) for f in inputs}
 master={}
 with (args.base/'inputs/master.tsv').open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   s,e,s0,e0=map(int,[r['gangstr_start_1based'],r['gangstr_end_1based'],r['trgt_start_0based'],r['trgt_end_0based']])
   if s0!=s-1 or e0!=e or e0-s0!=int(r['reference_length_bp']): raise ValueError('Master coordinate mismatch')
   lid=r['locus_id']
   if lid in master: raise ValueError('Duplicate locus ID')
   master[lid]=(r['chrom'],s,e,r['motif'])
 cols='key locus_id chrom start end motif motif_length reference_bp reference_copies gang_gt gang_repcn gang_ci gang_dp gang_q trgt_gt trgt_mc trgt_al trgt_range trgt_sd trgt_ap'.split()
 loci={}; seen=set()
 with (args.base/'comparison/all_calls.tsv').open() as f:
  for r in csv.DictReader(f,fieldnames=cols,delimiter='\t'):
   lid=r['locus_id']; s,e=int(r['start']),int(r['end'])
   if lid in seen or master.get(lid)!=(r['chrom'],s,e,r['motif']): raise ValueError('Calls/master mismatch')
   seen.add(lid)
   if int(r['reference_bp'])!=e-s+1: raise ValueError('Inclusive interval length mismatch')
   al,mc=pair(r['trgt_al']),pair(r['trgt_mc'])
   if not called(r['trgt_gt']) or al is None or mc is None: continue
   ml=int(r['motif_length'])
   if ml!=len(r['motif']) or ml<1: raise ValueError('Motif length mismatch')
   loci[lid]={'chrom':r['chrom'],'start0':s-1,'end0':e,'al':al,'mc':mc,'motif_length':ml,'gang':pair(r['gang_repcn']),'gang_gt':r['gang_gt'],'lengths':{},'key':r['key']}
 if seen!=set(master): raise ValueError('Missing master loci in all_calls.tsv')
 del master
 depths={}
 depthfile=comp/'analysis_loci.tsv.gz'
 if depthfile.exists():
  with gzip.open(depthfile,'rt') as f:
   for r in csv.DictReader(f,delimiter='\t'):
    try: depths[r['key']]=(float(r['short_depth']),float(r['long_depth']))
    except (ValueError,KeyError): pass
 print(f'Scanning supporting reads for {len(loci):,} TRGT-called loci; no callers are rerun.',flush=True)
 scanned=matched=qualified=0
 proc=subprocess.Popen([samtools,'view',str(bam)],stdout=subprocess.PIPE,text=True)
 try:
  for line in proc.stdout:
   scanned+=1; fields=line.rstrip('\n').split('\t')
   if scanned%2000000==0: print(f'Read {scanned:,} BAM records',flush=True)
   tr=next((s[5:] for s in fields[11:] if s.startswith('TR:Z:')),None)
   z=loci.get(tr)
   if z is None: continue
   matched+=1; size=read_length(fields,z,args)
   if size is None: continue
   qualified+=1
   molecule='/'.join(fields[0].split('/')[:2]); old=z['lengths'].get(molecule)
   z['lengths'][molecule]=size if old is None else (old if old==size else -1)
  if proc.wait()!=0: raise RuntimeError('samtools view failed')
 finally:
  proc.stdout.close()
  if proc.poll() is None: proc.terminate();proc.wait()
 if not matched: raise ValueError('No TR tags matched catalog IDs; do not infer coordinates from TR labels')
 summaries={name:[dict(length_bin=label,N_screened=0,N_supported=0,**dict.fromkeys(CATEGORIES,0)) for label in LABELS] for name in ['all','short_depth_higher']}
 detail=out/'locus_evidence.tsv.gz.tmp'
 with gzip.open(detail,'wt') as f:
  writer=csv.writer(f,delimiter='\t'); writer.writerow(['locus_id','key','longer_TRGT_allele_bp','TRGT_allele_bp','TRGT_motif_counts','GangSTR_motif_counts','N_distinct_qualifying_molecules','support_short_length','support_long_length','consistent_read_fraction','benchmark_pass','GangSTR_outcome','short_depth','long_depth'])
  for lid,z in loci.items():
   lengths=list(z.pop('lengths').values()); valid=[v for v in lengths if v>0]
   passed,counts,rate=support(valid,z,args)
   passed=passed and len(valid)==len(lengths) # conflicting records cannot establish this benchmark
   outcome=classify(z['gang'],z['mc'],z['gang_gt']) if passed else 'Unresolved / insufficient read support'
   size=max(z['al']); idx=next(i for i in range(len(LABELS)) if BINS[i]<=size<BINS[i+1])
   sd,ld=depths.get(z['key'],(float('nan'),float('nan')))
   for name in ['all']+(['short_depth_higher'] if math.isfinite(sd) and math.isfinite(ld) and sd>ld else []):
    row=summaries[name][idx];row['N_screened']+=1
    if passed: row['N_supported']+=1;row[outcome]+=1
   writer.writerow([lid,z['key'],size,z['al'],z['mc'],z['gang'],len(valid),*counts,rate,passed,outcome,sd,ld])
 for name,summary in summaries.items():
  with (out/(name+'_summary.tsv')).open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)
 after={str(f.resolve()): (f.stat().st_size,f.stat().st_mtime_ns) for f in inputs}
 if before!=after: raise RuntimeError('An evidence input changed during the BAM scan; rerun the analysis.')
 detail.replace(out/'locus_evidence.tsv.gz')
 (out/'provenance.json').write_text(json.dumps({'bam':str(bam),'all_calls':str(args.base/'comparison/all_calls.tsv'),'records_scanned':scanned,'matching_TR_records':matched,'qualifying_records_before_dedup':qualified,'criteria':vars(args),'input_size_mtime_ns':before,'coordinate_convention':'SAM POS and GangSTR start are 1-based; normalized target intervals are 0-based half-open; TR tags join exact catalog IDs.'},default=str,indent=2)+'\n')
 (out/'evidence_SUCCESS').write_text('complete\n')
 print('Read-support evidence saved in '+str(out),flush=True)
def resolve_spanning_bam(base, requested, samtools):
 if requested is not None:
  bam=requested
 else:
  candidates=sorted((base/'trgt').rglob('*spanning*.bam'))
  if len(candidates)!=1:
   raise SystemExit(f'Expected one TRGT spanning BAM, found {len(candidates)}. Supply --trgt-spanning-bam /absolute/path/to/file.bam. Candidates: {candidates}')
  bam=candidates[0]
 if not bam.is_file(): raise SystemExit('Missing supporting-read BAM: '+str(bam))
 subprocess.run([samtools,'quickcheck','-v',str(bam)],check=True)
 return bam
def render_read_supported_plots(output, sample="NEUAD700YFB"):
 with plt.rc_context(matplotlib.rcParamsDefault):
  P=Path(output)
  d=pd.read_csv(P/"read_supported_STR_benchmark/locus_evidence.tsv.gz",sep="\t")
  E=list(range(0,201,20))+[250,300,350,400,450,500,np.inf]
  L=[f"{a}-{int(b)-1}" if np.isfinite(b) else f">={a}" for a,b in zip(E,E[1:])]
  C=["Exact genotype match","Within 1 repeat on both alleles","Longer allele underestimated","Longer allele overestimated","Other genotype disagreement","GangSTR no call","GangSTR count unavailable"]
  T=["Exact genotype match","Within 1 repeat on both alleles (excluding exact)","GangSTR longer motif count lower","GangSTR longer motif count higher","Other genotype disagreement","GangSTR no call","GangSTR count unavailable"]
  K=["#246B45","#8BC898","#B52B38","#8B65BA","#E6A23C","#A7A7A7","#454545"]
  for c in ["longer_TRGT_allele_bp","short_depth","long_depth"]: d[c]=pd.to_numeric(d[c],errors="coerce")
  assert d.longer_TRGT_allele_bp.notna().all() and d.longer_TRGT_allele_bp.ge(0).all(),"Invalid allele length"
  d["bin"]=pd.cut(d.longer_TRGT_allele_bp,E,labels=False,right=False)
  d["passed"]=d.benchmark_pass.astype(str).str.lower().eq("true")
  assert d.loc[d.passed,"GangSTR_outcome"].isin(C).all(),"Unrecognized outcome"
  def pair(s):
   try:
    v=ast.literal_eval(str(s))
    return sorted(map(float,v)) if isinstance(v,(tuple,list)) and len(v)==2 and all(np.isfinite(float(x)) and float(x)>=0 for x in v) else None
   except (ValueError,TypeError,SyntaxError): return None
  def fit(x,y):
   X=np.log1p(x);Y=np.log1p(y)
   if len(X)<3 or np.ptp(X)==0 or np.ptp(Y)==0: return None
   b,a=np.polyfit(X,Y,1);pred=a+b*X
   return a,b,1-np.sum((Y-pred)**2)/np.sum((Y-Y.mean())**2)
  def save(fig,stem):
   fig.savefig(P/(stem+".png"),dpi=180);plt.close(fig)
  for name,mask in [("all",pd.Series(True,index=d.index)),("short_depth_higher",np.isfinite(d.short_depth)&np.isfinite(d.long_depth)&d.short_depth.gt(d.long_depth))]:
   z=d.loc[mask];q=z.loc[z.passed].copy()
   if z.empty: continue
   scope="Full locus set" if name=="all" else "Loci with short-read coverage exceeding long-read coverage"
   screened=z.groupby("bin").size().reindex(range(len(L)),fill_value=0).to_numpy()
   passed=q.groupby("bin").size().reindex(range(len(L)),fill_value=0).to_numpy()
   counts=pd.crosstab(q["bin"],q.GangSTR_outcome).reindex(index=range(len(L)),columns=C,fill_value=0)
   fig,ax=plt.subplots(figsize=(20,9));fig.subplots_adjust(left=.06,right=.99,bottom=.26,top=.74)
   fig.suptitle("STR genotype concordance by allele length",fontsize=22,y=.98)
   fig.text(.06,.93,sample+" | "+scope,fontsize=12)
   fig.text(.06,.90,"N = read-supported loci / screened loci; percentages are calculated among read-supported loci.",fontsize=11)
   bottom=np.zeros(len(L));x=np.arange(len(L))
   for cat,label,color in zip(C,T,K):
    h=np.divide(counts[cat].to_numpy(),passed,out=np.zeros(len(L)),where=passed>0)
    ax.bar(x,h,bottom=bottom,color=color,label=label);bottom+=h
   for i,(a,b) in enumerate(zip(passed,screened)): ax.text(i,1.02,f"N={a:,}/\n{b:,}",ha="center",va="bottom",fontsize=8)
   ax.set_ylim(0,1.17);ax.set_yticks(np.linspace(0,1,6));ax.yaxis.set_major_formatter(PercentFormatter(1))
   ax.set_xticks(x,L,rotation=45,ha="right");ax.set_ylabel("Read-supported loci (%)");ax.set_xlabel("Longer TRGT allele length (bp)")
   ax.spines[["top","right"]].set_visible(False);ax.set_axisbelow(True);ax.grid(axis="y",alpha=.15)
   fig.legend(*ax.get_legend_handles_labels(),loc="upper center",bbox_to_anchor=(.52,.88),ncol=3,fontsize=10)
   fig.text(.06,.11,"Read support: >=3 molecules per distinct allele length (>=6 for equal lengths), MAPQ >=20, read accuracy >=99%, and >=50-bp flanks.\nAt least 80% of qualifying reads must match a TRGT allele length within max(2 bp, one motif unit). No Q/AP filter is applied.\nMutually exclusive outcomes compare sorted motif counts; within-one-repeat matches exclude exact matches. Unsupported loci are excluded.\nTRGT genotyping reads define the benchmark; concordance does not establish independent accuracy. Sparse bins have limited precision.",fontsize=10)
   save(fig,f"analysis_read_supported_STR_singlebar_{name}")
   q["g"]=q.GangSTR_motif_counts.map(pair);q["t"]=q.TRGT_motif_counts.map(pair)
   v=q.loc[q.g.notna()&q.t.notna()&~q.GangSTR_outcome.isin(C[-2:])].copy()
   g=np.array(v.g.tolist(),dtype=float).reshape(-1,2);t=np.array(v.t.tolist(),dtype=float).reshape(-1,2)
   fig,axs=plt.subplots(1,2,figsize=(16,9));fig.subplots_adjust(left=.07,right=.97,top=.80,bottom=.22,wspace=.32)
   fig.suptitle("STR motif-count concordance between GangSTR and TRGT",fontsize=20,y=.98)
   fig.text(.07,.93,sample+" | "+scope,fontsize=11)
   fig.text(.07,.89,f"Read-supported loci with paired counts: {len(v):,} | Excluded for unavailable counts: {len(q)-len(v):,}",fontsize=11)
   for j,ax in enumerate(axs):
    xx,yy=t[:,j],g[:,j];hi=max(1,xx.max(),yy.max()) if len(xx) else 1
    ax.set_xscale("function",functions=(np.log1p,np.expm1));ax.set_yscale("function",functions=(np.log1p,np.expm1))
    if len(xx):
     xy,n=np.unique(np.column_stack([xx,yy]),axis=0,return_counts=True);o=np.argsort(n)
     sc=ax.scatter(xy[o,0],xy[o,1],c=n[o],s=16,cmap="viridis",norm=LogNorm(vmin=1,vmax=max(2,n.max())),linewidths=0)
     fig.colorbar(sc,ax=ax,pad=.025).set_label("Number of loci (log color scale)")
    ax.plot([0,hi],[0,hi],color="black",linestyle="--",linewidth=1,label="Equal motif counts")
    result=fit(xx,yy)
    if result:
     a,b,r2=result;grid=np.linspace(0,np.log1p(hi),200)
     ax.plot(np.expm1(grid),np.expm1(a+b*grid),color="red",linestyle=":",linewidth=1.7,label="Linear fit in log(1 + count)")
     text=f"R² = {r2:.4f}\nFit in log(1 + count) coordinates"
    else: text="R² unavailable: insufficient variation"
    ax.text(.97,.035,text,transform=ax.transAxes,ha="right",va="bottom",color="red",fontsize=10,bbox=dict(facecolor="white",alpha=.85,edgecolor="none"))
    ticks=[k for k in [0,1,10,100,1000,10000,100000] if k<=hi*1.05]
    ax.set_xticks(ticks,[f"{k:,}" for k in ticks]);ax.set_yticks(ticks,[f"{k:,}" for k in ticks])
    ax.set_xlim(0,hi*1.05);ax.set_ylim(0,hi*1.05);ax.set_aspect("equal",adjustable="box")
    ax.set_xlabel("TRGT motif count (repeat units)");ax.set_ylabel("GangSTR motif count (repeat units)")
    ax.set_title("Shorter allele" if j==0 else "Longer allele");ax.legend(loc="upper left",fontsize=8);ax.grid(alpha=.15)
   fig.text(.07,.11,"Alleles are paired by ascending motif count at each locus. Color represents loci with identical count pairs.\nAxes and ordinary least-squares regression use log(1 + count); tick labels show repeat counts. Each locus has equal regression weight.\nR² measures association, not agreement or accuracy. The benchmark reuses TRGT genotyping reads; missing calls are excluded.",fontsize=10)
   save(fig,f"analysis_read_supported_STR_scatter_{name}")
   v["delta"]=g[:,1]-t[:,1]
   stats=v.groupby("bin").delta.agg(N="size",median="median",q25=lambda s:s.quantile(.25),q75=lambda s:s.quantile(.75)).reindex(range(len(L)))
   fig,ax=plt.subplots(figsize=(16,7));fig.subplots_adjust(left=.09,right=.98,top=.78,bottom=.30)
   fig.suptitle("Difference in longer-allele motif count by STR length",fontsize=20,y=.98)
   fig.text(.09,.92,sample+" | "+scope,fontsize=11)
   fig.text(.09,.87,"Points show the median; vertical bars span the middle 50% of loci in each length bin (not confidence intervals).",fontsize=11)
   for i,r in stats.iterrows():
    if pd.isna(r.N): continue
    ax.errorbar(i,r["median"],yerr=[[r["median"]-r.q25],[r.q75-r["median"]]],fmt="o",color="#B52B38",capsize=4)
    ax.text(i,.97,f"N={int(r.N):,}",transform=ax.get_xaxis_transform(),ha="center",va="top",fontsize=8,rotation=65)
   ax.axhline(0,color="black",linestyle="--",linewidth=1);ax.margins(y=.2)
   ax.set_xticks(range(len(L)),L,rotation=45,ha="right");ax.set_xlabel("Longer TRGT allele length (bp)")
   ax.set_ylabel("GangSTR minus TRGT motif count\n(repeat units; longer sorted allele)");ax.grid(axis="y",alpha=.15)
   fig.text(.09,.10,"Negative values indicate fewer repeats reported by GangSTR; positive values indicate more. Zero indicates equal longer-allele counts.\nOnly read-supported loci with counts from both callers are included. Alleles are paired by count rank; TRGT is not independent ground truth.",fontsize=10)
   save(fig,f"analysis_read_supported_STR_difference_{name}")
  print("PNG figures saved in",P)
BASE = Path('/home/zw529/donglab/data/target_ALS/WGS_LR/repeat_comparison_gangSTR_vs_TRGT')
p = argparse.ArgumentParser(description='GangSTR versus TRGT analyses 1-3.')
p.add_argument('--base', type=Path, default=BASE)
p.add_argument('--short-bam', type=Path, default=BASE.parent / 'NEUAD700YFB.SD-029-24-CBLL.2.bam')
p.add_argument('--long-bam', type=Path, default=BASE.parent / 'NEUAD700YFB.SD-029-24-CBLL.3.mapped.bam')
p.add_argument('--out', type=Path, help='Output directory; defaults to BASE/comparison.')
p.add_argument('--trgt-spanning-bam', type=Path, help='TRGT supporting-read BAM; auto-detected under BASE/trgt if unique.')
p.add_argument('--read-support-only', action='store_true', help='Generate supporting-read evidence and six PNGs using existing analysis_loci.tsv.gz; do not recalculate BAM coverage.')
p.add_argument('--sample', default='NEUAD700YFB', help='Sample label for read-supported figures.')
p.add_argument('--plots-only', action='store_true', help='Regenerate plots from saved analysis tables without reading BAMs.')
p.add_argument('--calls-only', action='store_true', help='Run caller summaries only; skip BAM depth and read-support generation.')
p.add_argument('--samtools', default='samtools')
p.add_argument('--allow-sample-id-mismatch', action='store_true')
p.add_argument('--mapq', type=int, default=20)
p.add_argument('--baseq', type=int, default=0)
p.add_argument('--expected-master', type=int, default=784468)
p.add_argument('--expected-groups', default='735319,30151,5729,9804', help='Exact,(0,1],(1,2],>2; use actual expected counts for test data.')
a = p.parse_args()
if a.out is None: a.out=a.base/"comparison"
if a.read_support_only and (a.plots_only or a.calls_only): p.error("--read-support-only cannot be combined with --plots-only or --calls-only")
if min(a.mapq, a.baseq) < 0: p.error('Quality cutoffs must be nonnegative')
if not a.calls_only and (not a.short_bam or not a.long_bam):
 p.error('--short-bam and --long-bam are required unless --calls-only is explicit')
if a.plots_only:
 if not a.calls_only and not (a.out/'read_supported_STR_benchmark/locus_evidence.tsv.gz').is_file():
  p.error('Read-support evidence is missing. Run --read-support-only first, or run the full analysis.')
 render_plots(a.out, include_depth=not a.calls_only)
 if not a.calls_only: render_read_supported_plots(a.out,a.sample)
 sys.exit(0)
spanning_bam=None
if not a.calls_only:
 spanning_bam=resolve_spanning_bam(a.base,a.trgt_spanning_bam,a.samtools)
if a.read_support_only:
 if not (a.out/'analysis_loci.tsv.gz').is_file(): p.error('Run the full analysis first: analysis_loci.tsv.gz is missing.')
 generate_read_supported_evidence(a.base,a.out,spanning_bam,a.samtools)
 render_read_supported_plots(a.out,a.sample)
 sys.exit(0)
calls = a.base / 'comparison/all_calls.tsv'
master = a.base / 'inputs/master.tsv'
for path in (calls, master):
 if not path.is_file(): p.error('Missing input: ' + str(path))
a.out.mkdir(parents=True, exist_ok=True)
(a.out/'analysis_SUCCESS').unlink(missing_ok=True)
def sha(path):
 h = hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda: f.read(1048576), b''): h.update(b)
 return h.hexdigest()
provenance = {'arguments':vars(a).copy(), 'input_sha256':{str(x):sha(x) for x in (calls, master)},
    'python':sys.version, 'pandas':pd.__version__, 'numpy':np.__version__, 'commands':[]}
def save_provenance():
 (a.out/'analysis_provenance.json').write_text(json.dumps(provenance, indent=2, default=str)+'\n')
save_provenance()
def table(df, name):
 df.to_csv(a.out/('analysis_'+name+'.tsv'), sep='\t', index=False, na_rep='NA')
def pair(s):
 try:
  z = [float(v) for v in s.split(',')]
  return sorted(z) if len(z)==2 and all(math.isfinite(v) and v>=0 for v in z) else [np.nan,np.nan]
 except (ValueError, TypeError): return [np.nan,np.nan]
@functools.lru_cache(None)
def canonical(s):
 s=s.upper()
 if not s or set(s)-set('ACGT'): raise ValueError('Non-simple DNA motif: '+s)
 rc=s.translate(str.maketrans('ACGT','TGCA'))[::-1]
 return min(t[i:]+t[:i] for t in (s,rc) for i in range(len(t)))
def called(s):
 return bool(re.fullmatch(r'\d+[|/]\d+',s))
cols='key locus_id chrom start end motif motif_length reference_bp reference_copies gang_gt gang_repcn gang_ci gang_dp gang_q trgt_gt trgt_mc trgt_al trgt_range trgt_sd trgt_ap'.split()
print('Reading canonical calls and validating master coordinates', flush=True)
d=pd.read_csv(calls, sep='\t', header=None, names=cols, dtype=str, keep_default_na=False)
if d.shape[0]!=a.expected_master or d.isna().any().any():
 raise ValueError('Unexpected row count or malformed 20-column canonical table')
if d.key.duplicated().any(): raise ValueError('Duplicate coordinate keys')
for c in ['start','end','motif_length','reference_bp']:
 d[c]=pd.to_numeric(d[c],errors='raise').astype('int64')
m=pd.read_csv(master,sep='\t',dtype=str,keep_default_na=False)
required=['locus_id','chrom','gangstr_start_1based','gangstr_end_1based','trgt_start_0based','trgt_end_0based','motif_length','motif','reference_length_bp']
if list(m.columns)!=required: raise ValueError('Unexpected master schema')
for c in required[2:7]+['reference_length_bp']: m[c]=pd.to_numeric(m[c],errors='raise').astype('int64')
if not ((m.trgt_start_0based==m.gangstr_start_1based-1)&(m.trgt_end_0based==m.gangstr_end_1based)&(m.reference_length_bp==m.trgt_end_0based-m.trgt_start_0based)).all():
 raise ValueError('Master coordinate convention mismatch')
m['key']=m.chrom+':'+m.gangstr_start_1based.astype(str)+':'+m.gangstr_end_1based.astype(str)
if m.key.duplicated().any() or set(m.key)!=set(d.key): raise ValueError('Master/calls coordinate mismatch')
mm=m.set_index('key').loc[d.key]
for dc,mc in [('locus_id','locus_id'),('motif','motif'),('motif_length','motif_length'),('reference_bp','reference_length_bp')]:
 if not np.array_equal(d[dc].values, mm[mc].values): raise ValueError('Master mismatch: '+dc)
if not ((d.start>=1)&(d.end>=d.start)&(d.reference_bp==d.end-d.start+1)&(d.motif_length==d.motif.str.len())).all():
 raise ValueError('Invalid STR intervals or motif lengths')
d['gang_called']=d.gang_gt.map(called); d['trgt_called']=d.trgt_gt.map(called)
d['both_called']=d.gang_called & d.trgt_called
g=np.array(d.gang_repcn.map(pair).tolist()); t=np.array(d.trgt_mc.map(pair).tolist())
d['max_diff']=np.max(np.abs(g-t),axis=1)
d.loc[~d.both_called,'max_diff']=np.nan
d['gang_max_copies']=g[:,1]; d['trgt_max_copies']=t[:,1]
d['gang_max_bp']=g[:,1]*d.motif_length; d['trgt_max_motif_bp']=t[:,1]*d.motif_length
d['trgt_max_allele_bp']=np.array(d.trgt_al.map(pair).tolist())[:,1]
for c in ['gang_q','gang_dp']: d[c]=pd.to_numeric(d[c],errors='coerce')
d['trgt_min_ap']=np.array(d.trgt_ap.map(pair).tolist())[:,0]
d['trgt_min_sd']=np.array(d.trgt_sd.map(pair).tolist())[:,0]
for c in ['gang_q','trgt_min_ap']:
 d.loc[~np.isfinite(d[c])|~d[c].between(0,1),c]=np.nan
for c in ['gang_dp','trgt_min_sd']:
 d.loc[~np.isfinite(d[c])|(d[c]<0),c]=np.nan
for c in ['gang_q','gang_dp','gang_max_bp','gang_max_copies']: d.loc[~d.gang_called,c]=np.nan
for c in ['trgt_min_ap','trgt_min_sd','trgt_max_copies','trgt_max_motif_bp','trgt_max_allele_bp']: d.loc[~d.trgt_called,c]=np.nan
d['length_bin_start']=(d.reference_bp//10)*10
d['motif_length_group']=d.motif_length.map(lambda k:str(k) if k<=6 else '>6')
d['motif_family']=d.motif.map(canonical)
md=d.max_diff
d['agreement_group']=np.select([md==0,(md>0)&(md<=1),(md>1)&(md<=2),md>2],['exact','0_to_1','1_to_2','gt2'],default='not_comparable')
groups=['exact','0_to_1','1_to_2','gt2']
actual=[int((d.agreement_group==v).sum()) for v in groups]
expected=[int(x) for x in a.expected_groups.split(',')]
if actual!=expected or int(d.both_called.sum())!=sum(actual):
 raise ValueError('Canonical agreement check failed: '+str(actual)+' expected '+str(expected))
print('Canonical agreement verified: '+str(actual),flush=True)
table(pd.DataFrame({'group':groups,'N':actual}),'validated_agreement')
def summarize(z):
 out={'N_master':len(z),'N_gang_called':int(z.gang_called.sum()),'N_trgt_called':int(z.trgt_called.sum()),'N_both':int(z.both_called.sum())}
 v=z.max_diff.dropna(); n=len(v); out['N_comparable']=n
 for k,mask in [('exact',v==0),('le1',v<=1),('le2',v<=2),('gt2',v>2)]:
  out[k+'_N']=int(mask.sum()); out[k+'_rate']=mask.mean() if n else np.nan
 for grp in groups: out[grp+'_exclusive_N']=int((z.agreement_group==grp).sum())
 for c in ['reference_bp','gang_max_bp','trgt_max_allele_bp','gang_max_copies','trgt_max_copies','gang_q','trgt_min_ap','trgt_min_sd','gang_dp','short_depth','long_depth']:
  if c not in z: continue
  x=z[c].dropna(); out[c+'_N']=len(x); out[c+'_mean']=x.mean(); out[c+'_median']=x.median()
 for c in ['gang_q','trgt_min_ap']:
  v=z[c].dropna()
  for cut in [.9,.5]:
   key=c+'_lt'+str(cut)
   out[key+'_N']=int((v<cut).sum()); out[key+'_rate']=(v<cut).mean() if len(v) else np.nan
 v=z.trgt_min_sd.dropna(); out['trgt_min_sd_lt3_N']=int((v<3).sum()); out['trgt_min_sd_lt3_rate']=(v<3).mean() if len(v) else np.nan
 return out
def grouped(frame, keys, name):
 rows=[]
 for key,z in frame.groupby(keys,observed=True,sort=True):
  key=key if isinstance(key,tuple) else (key,)
  rows.append(dict(zip(keys,key),**summarize(z)))
 out=pd.DataFrame(rows); table(out,name); return out
grouped(d,['length_bin_start'],'reference_length_10bp')
grouped(d,['motif_length_group'],'motif_length')
families=d.loc[d.max_diff.notna()].motif_family.value_counts()
eligible=families[families>=100].index
grouped(d[d.motif_family.isin(eligible)],['motif_length','motif_family'],'motif_families_Nge100')
grouped(d[d.motif_family.isin(eligible)],['motif_family','length_bin_start'],'motif_family_by_10bp')
grouped(d,['motif_length_group','length_bin_start'],'motif_length_by_10bp')
bins=[0,1,5,10,20,30,50,100,np.inf]
labels=['0-<1','1-<5','5-<10','10-<20','20-<30','30-<50','50-<100','>=100']
def depth_summary(col):
 d[col+'_bin']=pd.cut(d[col],bins,right=False,labels=labels)
 grouped(d,[col+'_bin'],'by_'+col)
depth_summary('gang_dp'); depth_summary('trgt_min_sd')
def bam_header(path):
 subprocess.run([a.samtools,'quickcheck','-v',str(path)],check=True)
 h=subprocess.check_output([a.samtools,'view','-H',str(path)],text=True)
 seq={}; samples=set(); sorted_ok=False
 for line in h.splitlines():
  fields=dict(x.split(':',1) for x in line.split('\t')[1:] if ':' in x)
  if line.startswith('@HD'): sorted_ok=fields.get('SO')=='coordinate'
  if line.startswith('@SQ'): seq[fields['SN']]=fields
  if line.startswith('@RG') and fields.get('SM'): samples.add(fields['SM'])
 if not sorted_ok: raise ValueError('BAM must declare coordinate sorting: '+str(path))
 if len(samples)!=1: raise ValueError('Require one identifiable sample per BAM: '+str(path))
 return seq,samples
def measure_depth(path,label):
 bed=a.out/'analysis_master_intervals.bed'
 cmd=[a.samtools,'depth','-b',str(bed),'-q',str(a.baseq),'-Q',str(a.mapq),'-G','3844','-s',str(path)]
 provenance['commands'].append(shlex.join(cmd)); save_provenance()
 bychrom={c:sorted(zip(z.start,z.end,z.index)) for c,z in d.groupby('chrom',sort=False)}
 sums=np.zeros(len(d),dtype=np.float64); covered=np.zeros(len(d),dtype=np.int64)
 cur=None; active=[]; intervals=[]; j=0; prev=0; seen=set()
 with (a.out/('analysis_'+label+'_samtools.stderr')).open('w') as err:
  proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err,text=True,bufsize=1048576)
  try:
   for line in proc.stdout:
    chrom,pos,val=line.rstrip().split('\t'); pos=int(pos); val=int(val)
    if chrom!=cur:
     if chrom in seen: raise ValueError('Unsorted samtools output')
     seen.add(chrom); cur=chrom; intervals=bychrom.get(chrom,[]); active=[]; j=0; prev=0
    if pos<=prev or val<0: raise ValueError('Invalid samtools depth output')
    prev=pos
    while j<len(intervals) and intervals[j][0]<=pos:
     active.append(intervals[j]); j+=1
    active=[r for r in active if r[1]>=pos]
    for start,end,idx in active:
     sums[idx]+=val; covered[idx]+=int(val>0)
   if proc.wait()!=0: raise RuntimeError('samtools failed; see '+label+'_samtools.stderr')
  except BaseException:
   proc.terminate(); proc.wait(); raise
 d[label+'_depth']=sums/d.reference_bp
 d[label+'_covered_fraction']=covered/d.reference_bp
if not a.calls_only:
 print('Validating BAMs and measuring depth (this can take hours)',flush=True)
 version=subprocess.check_output([a.samtools,'--version'],text=True)
 v=re.search(r'samtools (\d+)\.(\d+)',version)
 if not v or tuple(map(int,v.groups()))<(1,13): raise ValueError('samtools >=1.13 required')
 provenance['samtools_version']=version
 s,ss=bam_header(a.short_bam); l,ls=bam_header(a.long_bam)
 if ss!=ls:
  if not a.allow_sample_id_mismatch and (ss, ls) != ({'NEUAD700YFB_SD029_24_3'}, {'BioSample61'}):
   raise ValueError('BAM sample labels differ. If these are the intended paired BAMs, rerun with --allow-sample-id-mismatch: '+str((ss,ls)))
  print('Using supplied BAM pair despite different sample labels: '+str((ss,ls)), flush=True)
 provenance['bam_sample_ids']={'short':sorted(ss),'long':sorted(ls)}
 provenance['sample_id_mismatch_overridden']=bool(ss!=ls)
 for chrom,z in d.groupby('chrom'):
  if chrom not in s or chrom not in l: raise ValueError('Master contig absent from BAM: '+chrom)
  if s[chrom]['LN']!=l[chrom]['LN'] or z.end.max()>int(s[chrom]['LN']): raise ValueError('Reference lengths disagree: '+chrom)
  if s[chrom].get('M5') and l[chrom].get('M5') and s[chrom]['M5']!=l[chrom]['M5']: raise ValueError('Reference MD5 mismatch: '+chrom)
 provenance['reference_validation']='Contig names/lengths and M5 when present; absent M5 cannot establish sequence identity.'
 for path in [a.short_bam,a.long_bam]:
  st=path.stat(); provenance.setdefault('bam_files',[]).append({'path':str(path.resolve()),'size':st.st_size,'mtime_ns':st.st_mtime_ns})
 d.assign(start0=d.start-1)[['chrom','start0','end']].sort_values(['chrom','start0','end']).to_csv(a.out/'analysis_master_intervals.bed',sep='\t',header=False,index=False)
 measure_depth(a.short_bam,'short'); measure_depth(a.long_bam,'long')
 stats=[]
 for c in ['short_depth','long_depth']:
  x=d[c]; stats.append(dict(dataset=c,N=len(x),mean=x.mean(),median=x.median(),p10=x.quantile(.1),p25=x.quantile(.25),p75=x.quantile(.75),p90=x.quantile(.9),fraction_lt5=(x<5).mean(),fraction_lt10=(x<10).mean(),fraction_ge20=(x>=20).mean(),base_weighted_mean=np.average(x,weights=d.reference_bp)))
  depth_summary(c)
 table(pd.DataFrame(stats),'bam_depth_summary')
 d['short_minus_long_depth']=d.short_depth-d.long_depth
 d['short_to_long_depth_ratio']=d.short_depth/d.long_depth.replace(0,np.nan)
 table(d[['short_depth','long_depth','short_minus_long_depth','short_to_long_depth_ratio']].describe(percentiles=[.1,.25,.5,.75,.9]).reset_index(),'paired_depth_summary')
 grouped(d,['length_bin_start','motif_length_group','short_depth_bin'],'depth_length_motif_strata')
 grouped(d[d.motif_family.isin(eligible)],['motif_family','length_bin_start','short_depth_bin'],'motif_family_length_depth_strata')
 high=d[(d.trgt_min_ap>=.9)&(d.trgt_min_sd>=3)]
 grouped(high,['short_depth_bin'],'short_depth_high_TRGT_evidence')
 corrcols=['reference_bp','motif_length','short_depth','long_depth','gang_q','trgt_min_ap','trgt_min_sd','max_diff']
 table(d[corrcols].corr(method='spearman').reset_index(),'spearman_correlations')
keep=['key','locus_id','chrom','start','end','motif','motif_length','motif_family','reference_bp','length_bin_start','gang_called','trgt_called','agreement_group','max_diff','gang_q','gang_dp','trgt_min_ap','trgt_min_sd','gang_max_copies','trgt_max_copies','gang_max_bp','trgt_max_motif_bp','trgt_max_allele_bp']
keep += [c for c in ['short_depth','long_depth','short_covered_fraction','long_covered_fraction','short_minus_long_depth','short_to_long_depth_ratio'] if c in d]
d[keep].to_csv(a.out/'analysis_loci.tsv.gz',sep='\t',index=False,na_rep='NA',compression='gzip')
render_plots(a.out, include_depth=not a.calls_only)
if not a.calls_only:
 generate_read_supported_evidence(a.base,a.out,spanning_bam,a.samtools)
 render_read_supported_plots(a.out,a.sample)
for path in (calls,master):
 if sha(path)!=provenance['input_sha256'][str(path)]: raise RuntimeError('Input changed during analysis: '+str(path))
provenance['status']='complete_calls_only' if a.calls_only else 'complete'
save_provenance()
(a.out/'analysis_SUCCESS').write_text(provenance['status']+'\n')
print('Finished: '+str(a.out),flush=True)
PY
