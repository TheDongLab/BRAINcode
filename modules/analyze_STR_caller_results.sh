#!/usr/bin/env bash
#SBATCH --job-name=analyze_STR_caller_results
#SBATCH --output=/home/zw529/donglab/data/target_ALS/WGS_LR/analyze_STR_caller_results.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/WGS_LR/analyze_STR_caller_results.err
#SBATCH --time=04:00:00
#SBATCH --mem=56G
#SBATCH --cpus-per-task=4

set -euo pipefail
module load SAMtools

exec "${PYTHON:-python3}" - "$@" <<'PY'
import argparse, functools, hashlib, json, math
from pathlib import Path
import re, shlex, subprocess, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path('/home/zw529/donglab/data/target_ALS/WGS_LR/repeat_comparison_gangSTR_vs_TRGT')
p = argparse.ArgumentParser(description='GangSTR versus TRGT analyses 1-3.')
p.add_argument('--base', type=Path, default=BASE)
p.add_argument('--short-bam', type=Path, default=BASE.parent / 'NEUAD700YFB.SD-029-24-CBLL.2.bam')
p.add_argument('--long-bam', type=Path, default=BASE.parent / 'NEUAD700YFB.SD-029-24-CBLL.3.mapped.bam')
p.add_argument('--out', type=Path, default=BASE / 'comparison')
p.add_argument('--calls-only', action='store_true', help='Explicitly run analyses 1-2 and caller-provided DP/SD summaries; skip BAM depth.')
p.add_argument('--samtools', default='samtools')
p.add_argument('--allow-sample-id-mismatch', action='store_true')
p.add_argument('--mapq', type=int, default=20)
p.add_argument('--baseq', type=int, default=0)
p.add_argument('--expected-master', type=int, default=784468)
p.add_argument('--expected-groups', default='735319,30151,5729,9804', help='Exact,(0,1],(1,2],>2; use actual expected counts for test data.')
a = p.parse_args()
if min(a.mapq, a.baseq) < 0: p.error('Quality cutoffs must be nonnegative')
if not a.calls_only and (not a.short_bam or not a.long_bam):
    p.error('--short-bam and --long-bam are required unless --calls-only is explicit')
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

length=grouped(d,['length_bin_start'],'reference_length_10bp')
motlen=grouped(d,['motif_length_group'],'motif_length')
families=d.loc[d.max_diff.notna()].motif_family.value_counts()
eligible=families[families>=100].index
fam=grouped(d[d.motif_family.isin(eligible)],['motif_length','motif_family'],'motif_families_Nge100')
grouped(d[d.motif_family.isin(eligible)],['motif_family','length_bin_start'],'motif_family_by_10bp')
grouped(d,['motif_length_group','length_bin_start'],'motif_length_by_10bp')

def savefig(name):
    plt.tight_layout(); plt.savefig(a.out/('analysis_'+name+'.png'),dpi=180); plt.close()

fig,ax=plt.subplots(2,1,figsize=(11,8),sharex=True)
for c in ['exact_rate','le1_rate','le2_rate']: ax[0].plot(length.length_bin_start+5,length[c],'.-',label=c)
ax[0].set(ylabel='Agreement fraction',ylim=(0,1.02)); ax[0].legend()
bottom=np.zeros(len(length))
for grp in groups:
    y=length[grp+'_exclusive_N']; ax[1].bar(length.length_bin_start,y,width=10,align='edge',bottom=bottom,label=grp); bottom+=y
ax[1].set(xlabel='Reference STR length (bp); bins [0,10), [10,20), ...',ylabel='Comparable loci'); ax[1].legend()
savefig('reference_length_10bp')

def motifplot(frame,key,name):
    if frame.empty: return
    f=frame.sort_values('N_comparable',ascending=False).head(30) if key=='motif_family' else frame
    fig,ax=plt.subplots(3,1,figsize=(max(9,len(f)*.32),10),sharex=True)
    xx=np.arange(len(f))
    ax[0].bar(xx,f.exact_rate); ax[0].set(ylabel='Exact agreement',ylim=(0,1))
    for c in ['reference_bp_median','gang_max_bp_median','trgt_max_allele_bp_median']: ax[1].plot(xx,f[c],'.-',label=c)
    ax[1].set(ylabel='Median length (bp)'); ax[1].legend(fontsize=8)
    for c in ['gang_q_lt0.9_rate','gang_q_lt0.5_rate','trgt_min_ap_lt0.9_rate','trgt_min_ap_lt0.5_rate']: ax[2].plot(xx,f[c],'.-',label=c)
    ax[2].set(ylabel='Fraction below cutoff',xticks=xx,xticklabels=f[key]); ax[2].tick_params(axis='x',rotation=90); ax[2].legend(fontsize=8)
    fig.suptitle('Q = genotype confidence; AP = allele purity (not equivalent accuracy)')
    savefig(name)
motifplot(motlen,'motif_length_group','motif_length')
motifplot(fam,'motif_family','motif_families_top30')

bins=[0,1,5,10,20,30,50,100,np.inf]
labels=['0-<1','1-<5','5-<10','10-<20','20-<30','30-<50','50-<100','>=100']
def depth_summary(col):
    d[col+'_bin']=pd.cut(d[col],bins,right=False,labels=labels)
    out=grouped(d,[col+'_bin'],'by_'+col)
    fig,ax=plt.subplots(2,1,figsize=(10,7),sharex=True)
    for c in ['exact_rate','le1_rate','gt2_rate']: ax[0].plot(out[col+'_bin'],out[c],'.-',label=c)
    ax[0].set(ylabel='Agreement fraction',ylim=(0,1)); ax[0].legend()
    for c in ['gang_q_lt0.9_rate','gang_q_lt0.5_rate','trgt_min_ap_lt0.9_rate','trgt_min_ap_lt0.5_rate']: ax[1].plot(out[col+'_bin'],out[c],'.-',label=c)
    ax[1].set(xlabel=col,ylabel='Fraction below cutoff'); ax[1].legend(fontsize=8)
    savefig('by_'+col)
depth_summary('gang_dp'); depth_summary('trgt_min_sd')

def bam_header(path):
    subprocess.run([a.samtools,'quickcheck','-v',str(path)],check=True)
    h=subprocess.check_output([a.samtools,'view','-H',str(path)],text=True)
    seq={}; samples=set(); order=[]; sorted_ok=False
    for line in h.splitlines():
        fields=dict(x.split(':',1) for x in line.split('\t')[1:] if ':' in x)
        if line.startswith('@HD'): sorted_ok=fields.get('SO')=='coordinate'
        if line.startswith('@SQ'): seq[fields['SN']]=fields; order.append(fields['SN'])
        if line.startswith('@RG') and fields.get('SM'): samples.add(fields['SM'])
    if not sorted_ok: raise ValueError('BAM must declare coordinate sorting: '+str(path))
    if len(samples)!=1: raise ValueError('Require one identifiable sample per BAM: '+str(path))
    return seq,samples,order

def measure_depth(path,label,seq):
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
    s,ss,so=bam_header(a.short_bam); l,ls,lo=bam_header(a.long_bam)
    if ss!=ls:
        if not a.allow_sample_id_mismatch:
            raise ValueError('BAM sample labels differ. If these are the intended paired BAMs, rerun with --allow-sample-id-mismatch: '+str((ss,ls)))
        print('Using supplied BAM pair despite different sample labels: '+str((ss,ls)), flush=True)
    provenance['bam_sample_ids']={'short':sorted(ss),'long':sorted(ls)}
    provenance['sample_id_mismatch_overridden']=bool(ss!=ls and a.allow_sample_id_mismatch)
    for chrom,z in d.groupby('chrom'):
        if chrom not in s or chrom not in l: raise ValueError('Master contig absent from BAM: '+chrom)
        if s[chrom]['LN']!=l[chrom]['LN'] or z.end.max()>int(s[chrom]['LN']): raise ValueError('Reference lengths disagree: '+chrom)
        if s[chrom].get('M5') and l[chrom].get('M5') and s[chrom]['M5']!=l[chrom]['M5']: raise ValueError('Reference MD5 mismatch: '+chrom)
    provenance['reference_validation']='Contig names/lengths and M5 when present; absent M5 cannot establish sequence identity.'
    for path in [a.short_bam,a.long_bam]:
        st=path.stat(); provenance.setdefault('bam_files',[]).append({'path':str(path.resolve()),'size':st.st_size,'mtime_ns':st.st_mtime_ns})
    d.assign(start0=d.start-1)[['chrom','start0','end']].sort_values(['chrom','start0','end']).to_csv(a.out/'analysis_master_intervals.bed',sep='\t',header=False,index=False)
    measure_depth(a.short_bam,'short',s); measure_depth(a.long_bam,'long',l)
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
    plt.figure(figsize=(7,6)); plt.hexbin(np.log1p(d.short_depth),np.log1p(d.long_depth),gridsize=70,bins='log',mincnt=1)
    lim=max(np.log1p(d.short_depth).max(),np.log1p(d.long_depth).max()); plt.plot([0,lim],[0,lim],'k--')
    plt.xlabel('log(1 + short-read mean depth)'); plt.ylabel('log(1 + long-read mean depth)'); plt.colorbar(label='Locus count'); savefig('paired_BAM_depth')
    corrcols=['reference_bp','motif_length','short_depth','long_depth','gang_q','trgt_min_ap','trgt_min_sd','max_diff']
    table(d[corrcols].corr(method='spearman').reset_index(),'spearman_correlations')

keep=['key','locus_id','chrom','start','end','motif','motif_length','motif_family','reference_bp','length_bin_start','gang_called','trgt_called','agreement_group','max_diff','gang_q','gang_dp','trgt_min_ap','trgt_min_sd','gang_max_copies','trgt_max_copies','gang_max_bp','trgt_max_motif_bp','trgt_max_allele_bp']
keep += [c for c in ['short_depth','long_depth','short_covered_fraction','long_covered_fraction','short_minus_long_depth','short_to_long_depth_ratio'] if c in d]
d[keep].to_csv(a.out/'analysis_loci.tsv.gz',sep='\t',index=False,na_rep='NA',compression='gzip')
for path in (calls,master):
    if sha(path)!=provenance['input_sha256'][str(path)]: raise RuntimeError('Input changed during analysis: '+str(path))
provenance['status']='complete_calls_only' if a.calls_only else 'complete'
save_provenance()
(a.out/'analysis_SUCCESS').write_text(provenance['status']+'\n')
print('Finished: '+str(a.out),flush=True)
PY
