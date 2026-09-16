#!/usr/bin/env bash
#SBATCH --job-name=downsample_seq_depth_for_STR_callers
#SBATCH --output=/home/zw529/donglab/data/target_ALS/WGS_LR/downsample_STR_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/WGS_LR/downsample_STR_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --partition=week
#SBATCH --mem=56G
#SBATCH --cpus-per-task=4

set -euo pipefail
module load SAMtools
export STR_DOWNSAMPLE_SCRIPT="$(readlink -f "$0")"
exec "${PYTHON:-$HOME/donglab/pipelines/modules/miniconda3/bin/python}" - "$@" <<'PY'
import argparse, collections, gzip, hashlib, json, math, os, re, shlex
import shutil, subprocess, sys, tempfile, textwrap
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT=Path.home()/'donglab'
BASE=ROOT/'data/target_ALS/WGS_LR/repeat_comparison_gangSTR_vs_TRGT'
p=argparse.ArgumentParser(description='Repeated read downsampling against fixed full-depth TRGT calls.')
p.add_argument('--base',type=Path,default=BASE)
p.add_argument('--short-bam',type=Path,default=BASE.parent/'NEUAD700YFB.SD-029-24-CBLL.2.bam')
p.add_argument('--long-bam',type=Path,default=BASE.parent/'NEUAD700YFB.SD-029-24-CBLL.3.mapped.bam')
p.add_argument('--reference',type=Path,default=ROOT/'references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa')
p.add_argument('--gangstr',default=str(ROOT/'pipelines/modules/gangstr-env/bin/GangSTR'))
p.add_argument('--trgt',default=str(ROOT/'pipelines/modules/trgt-5.1.0/trgt'))
p.add_argument('--samtools',default='samtools')
p.add_argument('--out',type=Path)
p.add_argument('--replicates',type=int,default=5)
p.add_argument('--per-stratum',type=int,default=20)
p.add_argument('--fractions',default='0.8,0.6,0.4,0.3,0.2,0.1,0.05')
p.add_argument('--seed',type=int,default=20260911)
p.add_argument('--threads',type=int,default=int(os.environ.get('SLURM_CPUS_PER_TASK',4)))
p.add_argument('--scratch',type=Path,default=Path(os.environ.get('SLURM_TMPDIR',os.environ.get('TMPDIR',str(BASE.parent)))))
p.add_argument('--arm',choices=['both','short','long'],default='both')
p.add_argument('--prepare-only',action='store_true',help='Validate inputs and select the panel without running callers.')
a=p.parse_args()
fractions=sorted(set(float(x) for x in a.fractions.split(',')),reverse=True)
if not fractions or not all(0<f<1 for f in fractions): p.error('Fractions must be between 0 and 1, excluding endpoints.')
if min(a.replicates,a.per_stratum,a.threads)<1 or a.seed<1: p.error('Counts, threads and seed must be positive.')
if not a.out: a.out=a.base/'comparison'/('downsampling_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
if a.out.exists(): p.error('Output already exists; use a fresh --out directory.')
arms=['short','long'] if a.arm=='both' else [a.arm]

# Keep the whole WGS BAM: cropping reads around STRs changes GangSTR evidence/calibration.
def run(cmd,log=None,capture=False):
    cmd=list(map(str,cmd))
    print(datetime.now().isoformat(timespec='seconds')+' '+shlex.join(cmd),flush=True)
    if log:
        with Path(log).open('w') as f: subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
        return ''
    return subprocess.check_output(cmd,text=True) if capture else subprocess.run(cmd,check=True)

def num(x):
    try:
        v=float(x); return v if math.isfinite(v) else np.nan
    except (ValueError,TypeError): return np.nan

def pair(x):
    z=[num(v) for v in str(x).split(',')]
    return tuple(sorted(z)) if len(z)==2 and all(math.isfinite(v) and v>=0 for v in z) else None

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

def write(frame,name): frame.to_csv(a.out/(name+'.tsv'),sep='\t',index=False,na_rep='NA')

def executable(x):
    z=shutil.which(x)
    if not z: raise RuntimeError('Executable unavailable: '+x)
    return z

a.samtools=executable(a.samtools); a.gangstr=executable(a.gangstr); a.trgt=executable(a.trgt)
inputs=[a.base/'comparison/analysis_loci.tsv.gz',a.base/'comparison/all_calls.tsv',a.base/'inputs/master.tsv',a.base/'inputs/master.gangstr.bed',a.base/'inputs/master.trgt.bed',a.reference,Path(str(a.reference)+'.fai'),a.short_bam,a.long_bam]
for path in inputs:
    if not path.is_file() or path.stat().st_size==0: raise FileNotFoundError(path)
version=run([a.samtools,'--version'],capture=True)
v=re.search(r'samtools (\d+)\.(\d+)',version)
if not v or tuple(map(int,v.groups()))<(1,15): raise RuntimeError('samtools >=1.15 is required.')
a.out.mkdir(parents=True); (a.out/'runs').mkdir(); (a.out/'inputs').mkdir()
manifest={'parameters':vars(a),'samtools':version,'python':sys.version,'commands_policy':'Full WGS template downsampling; same seed across fractions within replicate; each fraction drawn directly from original BAM.'}
for tool in ['gangstr','trgt']: manifest[tool+'_version']=run([getattr(a,tool),'--version'],capture=True).strip()
manifest['input_hashes']={str(x):sha(x) for x in inputs[:5]}
if os.environ.get('STR_DOWNSAMPLE_SCRIPT'): manifest['script_sha256']=sha(os.environ['STR_DOWNSAMPLE_SCRIPT'])
manifest['bam_files']=[{'path':str(x),'bytes':x.stat().st_size,'mtime_ns':x.stat().st_mtime_ns} for x in [a.short_bam,a.long_bam]]

# Use existing indexes through private symlinks; never rename or index the source BAMs.
def link_bam(path,label):
    candidates=[Path(str(path)+'.bai'),path.with_suffix('.bai'),Path(str(path)+'.csi')]
    original='wgs_bam_NEUAD700YFB_NEUAD700YFB.SD-029-24-CBLL.2.bam.bai' if label=='short' else 'lr-wgs_bam_NEUAD700YFB_NEUAD700YFB.SD-029-24-CBLL.3.mapped.bam.bai'
    candidates.append(path.parent/original)
    index=next((x for x in candidates if x.is_file()),None)
    if not index: raise FileNotFoundError('No index found for '+str(path))
    target=a.out/'inputs'/(label+'.bam'); target.symlink_to(path.resolve())
    Path(str(target)+('.csi' if index.suffix=='.csi' else '.bai')).symlink_to(index.resolve())
    run([a.samtools,'quickcheck','-v',target])
    h=run([a.samtools,'view','-H',target],capture=True)
    seq={}; samples=set()
    for line in h.splitlines():
        tags=dict(x.split(':',1) for x in line.split('\t')[1:] if ':' in x)
        if line.startswith('@SQ'): seq[tags['SN']]=tags
        if line.startswith('@RG') and tags.get('SM'): samples.add(tags['SM'])
    if len(samples)!=1: raise ValueError('Each BAM must contain one sample; found '+str(samples))
    manifest[label+'_sample_labels']=sorted(samples)
    return target,seq
bam={}; headers={}
for arm,path in [('short',a.short_bam),('long',a.long_bam)]: bam[arm],headers[arm]=link_bam(path,arm)
fai=pd.read_csv(str(a.reference)+'.fai',sep='\t',header=None,index_col=0)

print('Selecting baseline Q=1 / concordant / strong-TRGT loci',flush=True)
d=pd.read_csv(inputs[0],sep='\t')
cols='key locus_id chrom start end motif motif_length reference_bp reference_copies gang_gt gang_repcn gang_ci gang_dp gang_q trgt_gt trgt_mc trgt_al trgt_range trgt_sd trgt_ap'.split()
c=pd.read_csv(inputs[1],sep='\t',header=None,names=cols,dtype=str,keep_default_na=False)
if c.key.duplicated().any() or d.key.duplicated().any() or set(c.key)!=set(d.key): raise ValueError('Canonical and analysis locus sets differ.')
c=c.set_index('key'); d=d.set_index('key',drop=False)
selected=d[(d.gang_q==1)&(d.agreement_group=='exact')&(d.trgt_min_ap>=.9)&(d.trgt_min_sd>=3)&(d.short_depth>0)&(d.long_depth>0)].copy()
selected['length_stratum']=pd.cut(selected.reference_bp,[0,25,50,100,150,250,500,np.inf],right=False).astype(str)
selected['motif_stratum']=selected.motif_length.map(lambda v:str(int(v)) if v<=6 else '>6')
selected['depth_stratum']=pd.cut(selected.short_depth,[0,10,20,30,40,60,np.inf],right=False).astype(str)
selected['stratum']=selected.length_stratum+'|'+selected.motif_stratum+'|'+selected.depth_stratum
strata=selected.groupby('stratum').size().rename('N_eligible').reset_index()
selected['selection_order']=[hashlib.sha256((str(a.seed)+'|'+k).encode()).hexdigest() for k in selected.key]
panel=selected.sort_values('selection_order').groupby('stratum',sort=True).head(a.per_stratum).copy()
if panel.empty: raise ValueError('No eligible loci.')
panel['benchmark_mc']=[','.join(map(str,pair(c.loc[k,'trgt_mc']) or ())) for k in panel.key]
panel=panel.sort_values(['chrom','start','end']).reset_index(drop=True)
master=pd.read_csv(inputs[2],sep='\t'); master['key']=master.chrom+':'+master.gangstr_start_1based.astype(str)+':'+master.gangstr_end_1based.astype(str)
if master.key.duplicated().any(): raise ValueError('Duplicate master coordinate.')
master=master.set_index('key')
for r in panel.itertuples():
    m=master.loc[r.key]
    if r.reference_bp!=r.end-r.start+1 or m.trgt_start_0based!=r.start-1 or m.trgt_end_0based!=r.end or m.motif.upper()!=r.motif.upper(): raise ValueError('Master mismatch: '+r.key)
    if pair(c.loc[r.key,'gang_repcn'])!=pair(r.benchmark_mc) or num(c.loc[r.key,'gang_q'])!=1 or pair(r.benchmark_mc) is None: raise ValueError('Canonical eligibility mismatch: '+r.key)
    ap=pair(c.loc[r.key,'trgt_ap']); sd=pair(c.loc[r.key,'trgt_sd'])
    if ap is None or sd is None or min(ap)<.9 or min(sd)<3: raise ValueError('Canonical TRGT evidence mismatch: '+r.key)
    for arm in ['short','long']:
        h=headers[arm].get(r.chrom,{})
        if int(h.get('LN',-1))!=int(fai.loc[r.chrom,1]): raise ValueError('BAM/reference mismatch: '+r.chrom)
    h1=headers['short'][r.chrom]; h2=headers['long'][r.chrom]
    if h1.get('M5') and h2.get('M5') and h1['M5']!=h2['M5']: raise ValueError('BAM sequence MD5 mismatch.')
write(panel,'selected_loci')
write(strata.merge(panel.groupby('stratum').size().rename('N_selected'),on='stratum',how='left').fillna(0),'selection_strata')
keys=set(panel.key)
def subset_catalog(src,dst,trgt=False):
    seen=set()
    with src.open() as f,dst.open('w') as o:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            z=line.rstrip().split('\t'); key=z[0]+':'+str(int(z[1])+(1 if trgt else 0))+':'+z[2]
            if key in keys:
                if key in seen: raise ValueError('Duplicate catalog coordinate: '+key)
                seen.add(key); o.write(line if line.endswith('\n') else line+'\n')
    if seen!=keys: raise ValueError('Catalog does not cover selected panel: '+str(src))
gbed=a.out/'inputs/selected.gangstr.bed'; tbed=a.out/'inputs/selected.trgt.bed'; depthbed=a.out/'inputs/selected.depth.bed'
subset_catalog(inputs[3],gbed); subset_catalog(inputs[4],tbed,True)
panel.assign(start0=panel.start-1)[['chrom','start0','end']].to_csv(depthbed,sep='\t',header=False,index=False)
manifest['selected_N']=len(panel)
(a.out/'manifest.json').write_text(json.dumps(manifest,indent=2,default=str)+'\n')
coverage_path=a.base/'comparison/analysis_bam_depth_summary.tsv'
if coverage_path.exists(): shutil.copyfile(coverage_path,a.out/'original_STR_depth_summary.tsv')
if a.prepare_only:
    print('Prepared '+str(len(panel))+' loci: '+str(a.out)); sys.exit(0)
a.scratch.mkdir(parents=True,exist_ok=True)
needed=max(x.stat().st_size for x in [a.short_bam,a.long_bam])*1.2
if shutil.disk_usage(a.scratch).free<needed: raise RuntimeError('Need approximately '+str(round(needed/1e9))+' GB free temporary space. Set --scratch to suitable storage.')

# Read only selected VCF records, retaining missing calls in the denominator.
def read_calls(path,arm):
    out={}; opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            z=line.rstrip().split('\t')
            if len(z)!=10: raise ValueError('Expected single-sample VCF: '+str(path))
            info=dict(x.split('=',1) for x in z[7].split(';') if '=' in x)
            key=z[0]+':'+str(int(z[1])+(1 if arm=='long' else 0))+':'+info['END']
            if key not in keys: continue
            if key in out: raise ValueError('Duplicate VCF coordinate: '+key)
            val=dict(zip(z[8].split(':'),z[9].split(':')))
            gt=val.get('GT','.')
            alleles=pair(val.get('REPCN' if arm=='short' else 'MC','.')) if re.fullmatch(r'\d+[/|]\d+',gt) else None
            ap=pair(val.get('AP','.')); sd=pair(val.get('SD','.'))
            out[key]={'alleles':alleles,'q':num(val.get('Q','.')),'dp':num(val.get('DP','.')),'ap':min(ap) if ap else np.nan,'sd':min(sd) if sd else np.nan}
    return out

def depth(b):
    sums=np.zeros(len(panel)); intervals={ch:sorted(zip(z.start,z.end,z.index)) for ch,z in panel.groupby('chrom')}
    cur=None; active=[]; todo=[]; j=0; prev=0; seen=set()
    with tempfile.TemporaryFile(mode='w+') as err:
        proc=subprocess.Popen([a.samtools,'depth','-b',str(depthbed),'-q','0','-Q','20','-G','3844','-s',str(b)],stdout=subprocess.PIPE,stderr=err,text=True)
        try:
            for line in proc.stdout:
                ch,pos,val=line.split(); pos=int(pos); val=int(val)
                if ch!=cur:
                    if ch in seen: raise ValueError('Unsorted depth output.')
                    seen.add(ch);cur=ch;todo=intervals.get(ch,[]);active=[];j=0;prev=0
                if pos<=prev: raise ValueError('Duplicate/unsorted depth position.')
                prev=pos
                while j<len(todo) and todo[j][0]<=pos: active.append(todo[j]); j+=1
                active=[r for r in active if r[1]>=pos]
                for start,end,idx in active: sums[idx]+=val
            if proc.wait():
                err.seek(0); raise RuntimeError(err.read())
        except BaseException:
            proc.terminate();proc.wait();raise
    return sums/panel.reference_bp.to_numpy()

rows=[]
def experiment(arm,b,frac,rep,folder):
    folder.mkdir(parents=True,exist_ok=True); prefix=folder/'calls'
    if arm=='short': cmd=[a.gangstr,'--bam',b,'--ref',a.reference,'--regions',gbed,'--out',prefix,'--seed',a.seed]
    else: cmd=[a.trgt,'genotype','--genome',a.reference,'--reads',b,'--repeats',tbed,'--output-prefix',prefix,'--preset','wgs','--threads',a.threads]
    run(cmd,folder/'caller.log')
    vcfs=[Path(str(prefix)+'.vcf'),Path(str(prefix)+'.vcf.gz')]
    vcf=next((v for v in vcfs if v.is_file()),None)
    if not vcf: raise RuntimeError('Caller produced no VCF: '+str(folder))
    calls=read_calls(vcf,arm); actual=depth(b); result=[]
    for i,r in enumerate(panel.itertuples()):
        v=calls.get(r.key,{}); al=v.get('alleles'); benchmark=pair(r.benchmark_mc)
        diff=max(abs(al[k]-benchmark[k]) for k in (0,1)) if al else np.nan
        q=v.get('q',np.nan); ap=v.get('ap',np.nan); sd=v.get('sd',np.nan)
        result.append(dict(key=r.key,stratum=r.stratum,length_stratum=r.length_stratum,motif_stratum=r.motif_stratum,arm=arm,fraction=frac,replicate=rep,depth=actual[i],called=al is not None,q=q,dp=v.get('dp',np.nan),ap=ap,sd=sd,max_diff=diff,exact=bool(al and diff==0),within1=bool(al and diff<=1),q1=bool(al and q==1),quality_pass=bool(al and (q>=.9 if arm=='short' else ap>=.9 and sd>=3)),alleles=','.join(map(str,al)) if al else '.'))
    out=pd.DataFrame(result); out['match_and_quality']=out.exact & out.quality_pass
    writepath=folder/'metrics.tsv';out.to_csv(writepath,sep='\t',index=False,na_rep='NA')
    rows.append(out);return out

baseline={}
for arm in arms: baseline[arm]=experiment(arm,bam[arm],1.0,0,a.out/'runs'/arm/'baseline')
stable=set(panel.key)
for arm,z in baseline.items(): stable &= set(z.loc[z.exact & (z.q1 if arm=='short' else z.quality_pass),'key'])
audit=panel[['key','stratum']].copy(); audit['baseline_stable']=audit.key.isin(stable);write(audit,'baseline_reproduction')
if not stable: raise RuntimeError('No loci reproduced the baseline; inspect baseline caller logs before downsampling.')
print('Baseline-stable panel: '+str(len(stable))+' / '+str(len(panel)),flush=True)
for rep in range(1,a.replicates+1):
    for frac in fractions:
        for arm in arms:
            folder=a.out/'runs'/arm/('rep'+str(rep)+'_f'+format(frac,'.8g'))
            with tempfile.TemporaryDirectory(prefix='STR_downsample_',dir=a.scratch) as tmp:
                sub=Path(tmp)/'sample.bam'
                run([a.samtools,'view','-@',max(0,a.threads-1),'-b','--subsample',frac,'--subsample-seed',a.seed+rep,'-o',sub,bam[arm]])
                run([a.samtools,'index','-@',max(0,a.threads-1),sub])
                experiment(arm,sub,frac,rep,folder)
    pd.concat(rows,ignore_index=True).to_csv(a.out/'all_trials.partial.tsv.gz',sep='\t',index=False,na_rep='NA')
r=pd.concat(rows,ignore_index=True);r['baseline_stable']=r.key.isin(stable)
r.to_csv(a.out/'all_trials.tsv.gz',sep='\t',index=False,na_rep='NA')
primary=r[r.baseline_stable].copy()

# Bootstrap entire loci, keeping all repeated draws for each locus together.
def summarize(z):
    med=lambda s:s.median() if s.notna().any() else np.nan
    out={'N_loci':z.key.nunique(),'N_trials':len(z),'median_depth':z.depth.median(),'mean_depth':z.depth.mean(),'N_called':int(z.called.sum()),'N_q':int(z.q.notna().sum()),'median_q':med(z.q),'median_ap':med(z.ap),'median_sd':med(z.sd)}
    for metric in ['called','exact','within1','q1','quality_pass','match_and_quality']: out['p_'+metric]=z[metric].mean()
    out['p_no_call']=1-out['p_called']
    x=z.groupby('key').match_and_quality.mean().to_numpy();rng=np.random.default_rng(a.seed)
    boot=np.array([rng.choice(x,len(x),replace=True).mean() for _ in range(500)])
    out['match_quality_locus_bootstrap_low']=np.quantile(boot,.025);out['match_quality_locus_bootstrap_high']=np.quantile(boot,.975)
    return out

def grouped(z,by,name):
    output=[]
    for k,t in z.groupby(by,observed=True):
        if not isinstance(k,tuple): k=(k,)
        output.append(dict(zip(by,k),**summarize(t)))
    result=pd.DataFrame(output);write(result,name);return result
summary=grouped(primary,['arm','fraction'],'depth_response')
grouped(r,['arm','fraction'],'depth_response_all_selected')
grouped(primary,['arm','length_stratum','fraction'],'depth_response_by_length')
primary['depth_bin']=pd.cut(primary.depth,[0,2,5,10,15,20,25,30,40,60,100,np.inf],right=False).astype(str)
grouped(primary,['arm','depth_bin'],'response_by_measured_depth')
rep_summary=grouped(primary,['arm','replicate','fraction'],'response_by_replicate')
thresholds=[]
for arm,z in summary.groupby('arm'):
    z=z.sort_values('fraction',ascending=False)
    sustained=True; candidates=[]
    for v in z.itertuples():
        sustained &= v.p_match_and_quality>=.95 and v.N_loci>=100
        if sustained: candidates.append(v)
    v=candidates[-1] if candidates else None
    thresholds.append({'arm':arm,'criterion':'At least 95% observed match+quality at this and every higher tested fraction; >=100 distinct loci','lowest_tested_fraction':v.fraction if v else np.nan,'median_measured_depth_at_fraction':v.median_depth if v else np.nan,'observed_probability':v.p_match_and_quality if v else np.nan,'status':'panel-specific observed threshold' if v else 'not established'})
write(pd.DataFrame(thresholds),'minimum_tested_depth')
loss=[]
for (arm,key,rep),z in primary[primary.fraction<1].groupby(['arm','key','replicate']):
    z=z.sort_values('fraction',ascending=False)
    for metric in (['q1','exact','match_and_quality'] if arm=='short' else ['exact','match_and_quality']):
        failed=z[~z[metric]]; f=failed.iloc[0] if len(failed) else None
        loss.append(dict(arm=arm,key=key,replicate=rep,metric=metric,first_failure_fraction=f.fraction if f is not None else np.nan,first_failure_depth=f.depth if f is not None else np.nan,failure_type=('no_call' if not f.called else 'metric_failed') if f is not None else 'no_failure_in_tested_range'))
write(pd.DataFrame(loss),'first_loss_by_locus')

plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
for arm,z in summary.groupby('arm'):
    z=z.sort_values('fraction');name='GangSTR (short reads)' if arm=='short' else 'TRGT (long reads)'
    fig,ax=plt.subplots(2,1,figsize=(13,10))
    fig.suptitle('How much read depth preserves STR calls? — '+name,fontsize=17,fontweight='bold',y=.98)
    fig.text(.08,.93,'NEUAD700YFB | '+str(len(stable))+' baseline-stable loci | '+str(a.replicates)+' random downsampling replicates\nFull-depth TRGT is the fixed benchmark; lines follow tested read-retention fractions.',va='top')
    for col,label in [('p_exact','Exact agreement with full-depth TRGT'),('p_match_and_quality','Exact agreement + quality threshold'),('p_no_call','No genotype returned')]:ax[0].plot(z.median_depth,z[col],'o-',label=label)
    if arm=='short':ax[0].plot(z.median_depth,z.p_q1,'o--',label='GangSTR Q still equals 1')
    ax[0].fill_between(z.median_depth,z.match_quality_locus_bootstrap_low,z.match_quality_locus_bootstrap_high,alpha=.15,label='95% locus-bootstrap interval: match + quality')
    ax[0].axhline(.95,color='gray',ls=':',label='95% observed success target');ax[0].set_ylabel('Fraction of selected loci / trials');ax[0].yaxis.set_major_formatter(PercentFormatter(1));ax[0].set_ylim(0,1.05);ax[0].legend(fontsize=9)
    metric='q' if arm=='short' else 'ap';label='GangSTR genotype confidence (Q)' if arm=='short' else 'TRGT minimum allele purity (AP)'
    for rep,t in primary[(primary.arm==arm)&(primary.fraction<1)].groupby('replicate'):
        t=t.groupby('fraction')[['depth',metric]].median().sort_index();ax[1].plot(t.depth,t[metric],alpha=.25,color='#0072B2')
    ax[1].plot(z.median_depth,z['median_'+metric],'o-',color='#0072B2',label='Median across available values');ax[1].set_ylabel(label);ax[1].set_ylim(0,1.05);ax[1].legend()
    for b in ax:b.set_xlabel('Measured BAM depth at STRs: median across panel (x)');b.grid(alpha=.2)
    quality='Q >= 0.9' if arm=='short' else 'minimum AP >= 0.9 and minimum allele support SD >= 3'
    note='STR = short tandem repeat. Match requires both sorted allele repeat counts to equal full-depth TRGT. Quality threshold: '+quality+'. No-calls count as failures above; unavailable quality values are omitted from the lower panel. Faint lines show replicate medians. Loci were selected from baseline Q=1, concordant calls, with strong TRGT evidence, across length/motif/depth strata. Results are conditional on this panel; they are not a genome-wide minimum-depth guarantee. TRGT self-agreement measures stability, not independent accuracy. Q and AP are not equivalent confidence scores.'
    fig.text(.08,.025,textwrap.fill(note,145),fontsize=9,va='bottom');fig.subplots_adjust(top=.82,bottom=.23,hspace=.35,left=.09,right=.98)
    fig.savefig(a.out/(arm+'_depth_response.png'),dpi=180);plt.close(fig)
text='''Full-depth TRGT remains fixed in both arms. Long-read downsampling measures stability relative to that callset, not independently validated accuracy.
Selection: canonical GangSTR Q=1, exact REPCN/MC agreement, TRGT minimum AP>=0.9 and SD>=3; balanced random selection across reference length, motif length and starting short-read depth. Selection counts and baseline exclusions are recorded.
Each fraction is sampled directly from the original whole-genome BAM. A fixed template-hash seed within a replicate gives nested read sets; different replicates use different seeds. Mates are retained or discarded together. There is no interval cropping or upsampling.
Caller arguments retain the original GangSTR seed and TRGT WGS preset. Library estimates are recomputed from each full-WGS downsample by GangSTR. Baseline reruns identify loci stable under the selected-catalog rerun; results for all originally selected loci are also retained.
Depth uses the full reference STR span, includes zeros, MAPQ>=20/BQ>=0, excludes flags 3844, excludes deletions, and suppresses overlapping-mate double counts. DP and SD remain distinct caller-specific evidence measures.
The minimum is the lowest TESTED retained fraction sustaining >=95% observed exact agreement plus quality at that and every higher fraction, with >=100 distinct loci. Its median depth is a panel descriptor, not a per-locus causal cutoff. Read the first-loss and measured-depth tables for locus heterogeneity. A return to Q=1 after an earlier failure is possible; first loss is not an irreversible threshold.
Locus-bootstrap intervals keep repeated measurements together but do not account for genomic spatial dependence; limited replicate count and selection of initially successful loci limit generalization. Q=1 refers to the reported VCF value and may reflect rounding. AP is purity, not genotype confidence. Nonlinear empirical probabilities are reported instead of an inappropriate linear probability model.
Temporary downsampled BAMs are removed after each run; VCFs, logs, per-locus results and manifests are retained. No source BAM, reference, catalog or canonical result is changed.
Sources: https://www.htslib.org/doc/samtools-view.html ; https://www.htslib.org/doc/samtools-depth.html ; https://github.com/gymreklab/GangSTR ; https://github.com/PacificBiosciences/trgt
'''
(a.out/'METHODS.txt').write_text(text)
for path,expected in manifest['input_hashes'].items():
    if sha(path)!=expected: raise RuntimeError('An input changed during the run: '+path)
(a.out/'SUCCESS').write_text('complete\n')
print('Finished: '+str(a.out),flush=True)
PY
