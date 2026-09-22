#!/usr/bin/env bash
#SBATCH --job-name=STR_flanks_TiTv
#SBATCH --output=/home/zw529/donglab/data/target_ALS/WGS_LR/repeat_comparison_gangSTR_vs_TRGT/logs/STR_flanks_TiTv_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/WGS_LR/repeat_comparison_gangSTR_vs_TRGT/logs/STR_flanks_TiTv_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail
# For sbatch, create this directory BEFORE submitting: Slurm opens logs before this script runs.
STR_LOG_DIR="${STR_LOG_DIR:-$HOME/donglab/data/target_ALS/WGS_LR/repeat_comparison_gangSTR_vs_TRGT/logs}"
mkdir -p "$STR_LOG_DIR"
STR_RUN_ID="${SLURM_JOB_ID:-$(date +%Y%m%d_%H%M%S)_$$}"
STR_OUT_LOG="$STR_LOG_DIR/STR_flanks_TiTv_${STR_RUN_ID}.out"
STR_ERR_LOG="$STR_LOG_DIR/STR_flanks_TiTv_${STR_RUN_ID}.err"
printf 'Output log: %s\nError log: %s\n' "$STR_OUT_LOG" "$STR_ERR_LOG"
exec >>"$STR_OUT_LOG" 2>>"$STR_ERR_LOG"
finish() {
    str_exit_status=$?
    trap - EXIT
    if [ "$str_exit_status" -eq 0 ]; then
        printf '[%s] COMPLETED successfully\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    else
        printf '[%s] FAILED, exit status %s; inspect %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$str_exit_status" "$STR_ERR_LOG"
        printf '[%s] FAILED, exit status %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$str_exit_status" >&2
    fi
    exit "$str_exit_status"
}
trap finish EXIT
printf '[%s] STARTED host=%s job=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(hostname)" "${SLURM_JOB_ID:-interactive}"
printf 'Arguments:'
printf ' %q' "$@"
printf '\n'
"${PYTHON:-$HOME/donglab/pipelines/modules/miniconda3/bin/python}" -u - "$@" <<'PYTHON_STR_FLANKS'
#!/usr/bin/env python3
"""Reference flank composition and gene-feature overlaps at read-supported STR loci."""
import argparse
import bisect
import collections
import gzip
import json
import re
from pathlib import Path

class IndexedFasta:
    """Read an uncompressed FASTA using its existing samtools .fai index."""
    def __init__(self, path):
        self.path = Path(path)
        if self.path.suffix == '.gz':
            raise ValueError('Use the uncompressed hg38 genome.fa configured in the caller pipeline.')
        self.index = {}
        with Path(str(path) + '.fai').open() as f:
            for line in f:
                chrom, length, offset, bases, width, *_ = line.rstrip().split('\t')
                self.index[chrom] = tuple(map(int, (length, offset, bases, width)))
        self.f = self.path.open('rb')

    def fetch(self, chrom, start, end):
        length, offset, bases, width = self.index[chrom]
        if not 0 <= start <= end <= length:
            raise ValueError(f'Out-of-range interval: {chrom}:{start}-{end}')
        if start == end:
            return ''
        first = offset + start // bases * width + start % bases
        last = offset + (end - 1) // bases * width + (end - 1) % bases + 1
        self.f.seek(first)
        seq = self.f.read(last-first).replace(b'\n', b'').replace(b'\r', b'').decode().upper()
        if len(seq) != end-start:
            raise ValueError('FASTA/index mismatch')
        return seq

    def flanks(self, chrom, start, end, n):
        length = self.index[chrom][0]
        left = self.fetch(chrom, max(0, start-n), start).rjust(n, 'N')
        right = self.fetch(chrom, end, min(length, end+n)).ljust(n, 'N')
        return left, right


def merged(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1][1] = max(end, result[-1][1])
        else:
            result.append([start, end])
    return result


class FeatureIndex:
    def __init__(self, intervals):
        self.spans = merged(intervals)
        self.starts = [x[0] for x in self.spans]

    def overlaps(self, start, end):
        i = bisect.bisect_left(self.starts, end)-1
        return i >= 0 and self.spans[i][1] > start


def load_gtf(path):
    """GTF coordinates are 1-based inclusive; internal intervals are half-open."""
    features = collections.defaultdict(lambda: collections.defaultdict(list))
    transcripts = collections.defaultdict(list)
    chroms = set()
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError('Expected a nine-column GTF file.')
            chrom, _, kind, s, e, _, strand, _, attributes = fields
            start, end = int(s)-1, int(e)
            if start < 0 or end <= start:
                raise ValueError('Invalid GTF interval')
            chroms.add(chrom)
            attrs = dict(re.findall(r'(\S+)\s+"([^"]*)"', attributes))
            if kind in ('gene', 'transcript'):
                features[chrom]['genic'].append((start, end))
            if kind == 'exon':
                features[chrom]['exon'].append((start, end))
                if 'transcript_id' not in attrs:
                    raise ValueError('Exons require transcript_id attributes; supply GTF, not GFF3.')
                transcripts[(chrom, strand, attrs['transcript_id'])].append((start, end))
    if not transcripts:
        raise ValueError('No transcript-associated exons found in annotation.')
    for (chrom, strand, tid), exons in transcripts.items():
        exons = merged(exons)
        features[chrom]['genic'].append((exons[0][0], exons[-1][1]))
        for left, right in zip(exons, exons[1:]):
            if left[1] < right[0]:
                features[chrom]['intron'].append((left[1], right[0]))
    return {c: {k: FeatureIndex(v) for k, v in fs.items()} for c, fs in features.items()}, chroms


def classify_feature(index, chroms, chrom, start, end):
    if chrom not in chroms:
        return 'Unannotated contig', False, False
    fs = index.get(chrom, {})
    hits = {k: v.overlaps(start, end) for k, v in fs.items()}
    exon, intron = hits.get('exon', False), hits.get('intron', False)
    # A region can be exonic in one isoform and intronic in another.
    if exon and intron:
        label = 'Exonic + intronic'
    elif exon:
        label = 'Exonic'
    elif intron:
        label = 'Intronic'
    elif hits.get('genic', False):
        label = 'Other genic'
    else:
        label = 'Intergenic'
    return label, exon, intron


def phase_extension(reference, motif, left, right):
    """Report in-phase periodic bases adjacent to catalog boundaries; do not shift loci."""
    motif = motif.upper()
    if not motif or any(b not in 'ACGT' for b in motif):
        return None, None, False
    aligned = all(b == motif[i % len(motif)] for i, b in enumerate(reference))
    if not aligned:
        return None, None, False
    a = b = 0
    for j, base in enumerate(reversed(left), 1):
        if base != motif[-j % len(motif)]:
            break
        a += 1
    for j, base in enumerate(right):
        if base != motif[(len(reference)+j) % len(motif)]:
            break
        b += 1
    return a, b, True


# Sample substitutions and literal TRGT-allele verification.
"""Sample-VCF substitutions in reference flanks; no inference of mutation direction or causality."""
import bisect
import collections
import gzip
import re


def substitution_class(ref, alt):
    if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT' or ref==alt:
        return None
    return 'Transition' if frozenset((ref,alt)) in (frozenset('AG'),frozenset('CT')) else 'Transversion'


def alternating_at(sequence):
    return bool(sequence) and all(b in 'AT' for b in sequence) and all(a!=b for a,b in zip(sequence,sequence[1:]))


def vcf_records(path, sample=None):
    opener=gzip.open if str(path).endswith('.gz') else open
    selected=None
    with opener(path,'rt') as f:
        for line in f:
            if line.startswith('##'): continue
            if line.startswith('#CHROM'):
                header=line.rstrip().split('\t');names=header[9:]
                if sample:
                    if sample not in names: raise ValueError(f'Sample {sample} not found in {path}: {names}')
                    selected=9+names.index(sample)
                elif len(names)==1: selected=9
                else: raise ValueError(f'Select one sample with --vcf-sample for {path}: {names}')
                continue
            if line.startswith('#'): continue
            if selected is None: raise ValueError('VCF sample header missing')
            z=line.rstrip().split('\t')
            if len(z)<=selected: raise ValueError('Truncated VCF record')
            values=dict(zip(z[8].split(':'),z[selected].split(':')))
            info=dict(v.split('=',1) for v in z[7].split(';') if '=' in v)
            yield dict(chrom=z[0],pos0=int(z[1])-1,ref=z[3].upper(),alts=z[4].upper().split(','),filter=z[6],format=values,info=info)


class FlankLookup:
    """Small fixed-width windows queried by their starts; overlapping loci are retained."""
    def __init__(self, loci, fasta, n):
        windows=collections.defaultdict(list)
        self.n=n
        for r in loci.itertuples():
            for side,s,e in [('left',max(0,r.start0-n),r.start0),('right',r.end0,min(fasta.index[r.chrom][0],r.end0+n))]:
                if e>s: windows[r.chrom].append((s,e,r.locus_id,side,r.start0,r.end0))
        self.windows={c:sorted(w) for c,w in windows.items()}
        self.starts={c:[w[0] for w in v] for c,v in self.windows.items()}

    def find(self, chrom, pos0):
        spans=self.windows.get(chrom,[]);starts=self.starts.get(chrom,[])
        a=bisect.bisect_left(starts,pos0-self.n+1);b=bisect.bisect_right(starts,pos0)
        return [w for w in spans[a:b] if w[0]<=pos0<w[1]]


def read_flank_variants(path, sample, loci, fasta, flank):
    lookup=FlankLookup(loci,fasta,flank)
    rows=[];seen=set();audit=collections.Counter()
    import time
    last_progress=time.monotonic()
    for rec in vcf_records(path,sample):
        audit['VCF_records_scanned']+=1
        if audit['VCF_records_scanned']%10000==0 and time.monotonic()-last_progress>=30:
            print('Scanned '+format(audit['VCF_records_scanned'],',')+' genomic VCF records; qualifying SNV alleles: '+str(audit['unique_qualifying_SNV_alleles']),flush=True)
            last_progress=time.monotonic()
        hits=lookup.find(rec['chrom'],rec['pos0'])
        if not hits: continue
        audit['records_overlapping_AT_TA_flanks']+=1
        if rec['filter']!='PASS':
            audit['excluded_non_PASS']+=1;continue
        fmt=rec['format'];gt=fmt.get('GT','.')
        if not re.fullmatch(r'\d+[/|]\d+',gt):
            audit['excluded_missing_or_non_diploid_GT']+=1;continue
        if fmt.get('FT','PASS') not in ('PASS','.'):
            audit['excluded_sample_filter']+=1;continue
        if len(rec['ref'])!=1 or rec['ref'] not in 'ACGT':
            audit['excluded_non_SNV_REF']+=1;continue
        observed=fasta.fetch(rec['chrom'],rec['pos0'],rec['pos0']+1)
        if observed!=rec['ref']: raise ValueError(f'VCF/hg38 REF mismatch at {rec["chrom"]}:{rec["pos0"]+1}')
        indices=[int(x) for x in re.split(r'[/|]',gt)]
        for i in sorted(set(indices)-{0}):
            if i>len(rec['alts']): raise ValueError('VCF genotype index exceeds ALT count')
            alt=rec['alts'][i-1];kind=substitution_class(rec['ref'],alt)
            if kind is None:
                audit['excluded_non_SNV_ALT']+=1;continue
            key=(rec['chrom'],rec['pos0'],rec['ref'],alt)
            if key in seen: continue
            seen.add(key);audit['unique_qualifying_SNV_alleles']+=1
            for s,e,lid,side,start,end in hits:
                rows.append(dict(locus_id=lid,chrom=rec['chrom'],position_1based=rec['pos0']+1,ref=rec['ref'],alt=alt,substitution=rec['ref']+'>'+alt,substitution_class=kind,side=side,relative_position=rec['pos0']-start if side=='left' else rec['pos0']-end+1,GT=gt,ALT_dosage=indices.count(i),GQ=fmt.get('GQ','.'),DP=fmt.get('DP','.'),PS=fmt.get('PS','.'),variant_key=':'.join(map(str,key))))
    return rows,dict(audit)


def trgt_sequence_purity(path, sample, loci, fasta):
    """Strip only validated reference flanks from literal VCF allele sequences."""
    by_id={str(r.locus_id):r for r in loci.itertuples()}
    by_interval={(r.chrom,r.start0,r.end0):r for r in loci.itertuples()}
    result={};seen=set()
    for rec in vcf_records(path,sample):
        info=rec['info'];r=by_id.get(info.get('TRID',''))
        # Original TRGT VCF POS includes the preceding anchor base in this pipeline.
        if r is None and 'END' in info:
            r=by_interval.get((rec['chrom'],rec['pos0']+1,int(info['END'])))
        if r is None: continue
        lid=str(r.locus_id)
        if lid in seen: raise ValueError('Duplicate TRGT records for '+lid)
        seen.add(lid)
        status='unverifiable';pure=None;lengths=None
        result[lid]=(status,pure,lengths)
        start,end=rec['pos0'],rec['pos0']+len(rec['ref'])
        if rec['chrom']!=r.chrom or not start<=r.start0<r.end0<=end: continue
        if fasta.fetch(r.chrom,start,end)!=rec['ref']: raise ValueError('TRGT VCF/hg38 REF mismatch for '+lid)
        gt=rec['format'].get('GT','.')
        if not re.fullmatch(r'\d+[/|]\d+',gt): continue
        indices=[int(v) for v in re.split(r'[/|]',gt)]
        alleles=[rec['ref']]+rec['alts']
        prefix=rec['ref'][:r.start0-start];suffix=rec['ref'][r.end0-start:]
        seqs=[]
        for i in indices:
            if i>=len(alleles): raise ValueError('TRGT GT index exceeds ALT count')
            seq=alleles[i]
            if any(b not in 'ACGT' for b in seq) or not seq.startswith(prefix) or (suffix and not seq.endswith(suffix)): break
            if len(seq)<len(prefix)+len(suffix): break
            seqs.append(seq[len(prefix):len(seq)-len(suffix) if suffix else None])
        if len(seqs)!=2: continue
        try: al=[int(v) for v in rec['format']['AL'].split(',')]
        except (KeyError,ValueError): continue
        lengths=[len(v) for v in seqs]
        if lengths!=al: continue
        # Verify the literal VCF alleles belong to the evidence callset.
        import ast
        expected=sorted(float(v) for v in ast.literal_eval(str(r.TRGT_allele_bp)))
        if sorted(lengths)!=expected: raise ValueError('TRGT VCF allele lengths differ from saved evidence for '+lid)
        pure=all(alternating_at(v) for v in seqs)
        result[lid]=('sequence-verified alternating AT' if pure else 'sequence-verified interrupted/non-AT',pure,lengths)
    return result


def analyze_substitutions(d, args, out):
    import ast
    import json
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    # IndexedFasta is supplied by the enclosing standalone analysis script.
    target=d.loc[d.motif.isin(['AT','TA'])].copy()
    if target.empty:
        raise ValueError('No read-supported exact AT or TA motifs in the selected data.')
    fasta=IndexedFasta(args.reference)
    print('Scanning sample genomic VCF for AT/TA flank substitutions: '+str(args.snv_vcf),flush=True)
    rows,audit=read_flank_variants(args.snv_vcf,args.vcf_sample,target,fasta,args.flank_bp)
    print('Genomic VCF scan finished: '+json.dumps(audit),flush=True)
    (out/'SNV_scan_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    columns=['locus_id','chrom','position_1based','ref','alt','substitution','substitution_class','side','relative_position','GT','ALT_dosage','GQ','DP','PS','variant_key']
    events=pd.DataFrame(rows,columns=columns)
    events=events.merge(target[['locus_id','motif','length_group']],on='locus_id',validate='many_to_one')
    events.to_csv(out/'AT_TA_flank_SNV_associations.tsv.gz',sep='\t',index=False)
    unique=events.drop_duplicates('variant_key')
    unique.to_csv(out/'AT_TA_unique_flank_SNVs.tsv',sep='\t',index=False)
    counts=events.groupby(['locus_id','substitution_class']).variant_key.nunique().unstack(fill_value=0)
    for kind in ['Transition','Transversion']:
        target[kind+'_N']=target.locus_id.map(counts[kind] if kind in counts else {}).fillna(0).astype(int)
    target['SNV_group']=np.select([(target.Transition_N>0)&(target.Transversion_N>0),target.Transition_N>0,target.Transversion_N>0],['Both types','Transition only','Transversion only'],default='No qualifying SNV recorded')
    target['longer_allele_minus_reference_bp']=target.longer_TRGT_allele_bp-target.reference_bp
    target['shorter_allele_minus_reference_bp']=[min(ast.literal_eval(str(v)))-ref for v,ref in zip(target.TRGT_allele_bp,target.reference_bp)]
    target['max_absolute_allele_length_difference_bp']=target[['shorter_allele_minus_reference_bp','longer_allele_minus_reference_bp']].abs().max(axis=1)
    catalog= pd.read_csv(args.comparison/'all_calls.tsv',sep='\t',header=None,names='key locus_id chrom start end motif motif_length reference_bp reference_copies gang_gt gang_repcn gang_ci gang_dp gang_q trgt_gt trgt_mc trgt_al trgt_range trgt_sd trgt_ap'.split(),dtype=str,keep_default_na=False)
    if catalog.locus_id.duplicated().any(): raise ValueError('Duplicate callset IDs')
    catalog=catalog.set_index('locus_id')
    def reported_pure(lid):
        if lid not in catalog.index: return None
        try:
            ap=[float(x) for x in catalog.loc[lid,'trgt_ap'].split(',')]
            lengths=sorted(float(x) for x in catalog.loc[lid,'trgt_al'].split(','))
            expected=sorted(float(x) for x in ast.literal_eval(str(target_by_id[lid])))
            if lengths!=expected: raise RuntimeError('Original calls differ from saved evidence at '+str(lid))
            return len(ap)==2 and all(v==1 for v in ap)
        except (ValueError,TypeError): return None
    target_by_id=dict(zip(target.locus_id,target.TRGT_allele_bp))
    target['TRGT_reported_AP1_both']=target.locus_id.map(reported_pure)
    vcf=args.trgt_vcf
    if vcf is None:
        candidates=sorted(set((args.comparison.parent/'trgt').rglob('*.vcf'))|set((args.comparison.parent/'trgt').rglob('*.vcf.gz')))
        if len(candidates)==1: vcf=candidates[0]
    if vcf is None: print('No unique original TRGT VCF found; AP-only examples will be labeled unverified. Supply --trgt-vcf to verify sequences.',flush=True)
    if vcf: print('Checking original TRGT allele sequences: '+str(vcf),flush=True)
    sequence=trgt_sequence_purity(vcf,args.trgt_sample,target,fasta) if vcf else {}
    print('TRGT sequence inspection finished; writing locus and example tables.',flush=True)
    target['sequence_status']=[sequence.get(str(lid),('not verified',None,None))[0] for lid in target.locus_id]
    target['both_alleles_alternating_AT']=[sequence.get(str(lid),('not verified',None,None))[1] for lid in target.locus_id]
    target['reference_alternating_AT']=target.reference_sequence.map(alternating_at)
    descriptions=events.groupby('locus_id')[['chrom','position_1based','substitution','side','relative_position']].apply(lambda z:'; '.join(f'{r.chrom}:{r.position_1based} {r.substitution} ({r.side}, {r.relative_position:+d})' for r in z.itertuples())) if len(events) else {}
    target['flank_substitutions']=target.locus_id.map(descriptions).fillna('')
    target.to_csv(out/'AT_TA_loci_substitutions_and_lengths.tsv.gz',sep='\t',index=False)
    seq_pure=target.both_alleles_alternating_AT.eq(True)
    selected=target.loc[(target.Transition_N+target.Transversion_N>0)&target.max_absolute_allele_length_difference_bp.gt(0)&seq_pure&target.reference_alternating_AT].copy()
    selected['absolute_length_difference_bp']=selected.max_absolute_allele_length_difference_bp
    selected=selected.sort_values('absolute_length_difference_bp',ascending=False)
    example_cols=['locus_id','chrom','start','end','motif','reference_bp','TRGT_allele_bp','shorter_allele_minus_reference_bp','longer_allele_minus_reference_bp','sequence_status','TRGT_reported_AP1_both','Transition_N','Transversion_N','flank_substitutions','genomic_feature','N_distinct_qualifying_molecules','consistent_read_fraction']
    selected[example_cols].to_csv(out/'AT_TA_sequence_verified_length_change_examples.tsv',sep='\t',index=False)
    selected[example_cols].head(20).to_csv(out/'AT_TA_top20_sequence_verified_examples.tsv',sep='\t',index=False)
    ap_only=target.loc[(target.Transition_N+target.Transversion_N>0)&target.max_absolute_allele_length_difference_bp.gt(0)&target.TRGT_reported_AP1_both.eq(True)&target.both_alleles_alternating_AT.isna()].copy()
    ap_only[example_cols].to_csv(out/'AT_TA_AP1_examples_not_sequence_verified.tsv',sep='\t',index=False)
    print('Generating transition/transversion figures...',flush=True)
    # Count unique observed REF/ALT SNV alleles, not allele dosage or length changes.
    spectrum=[a+'>'+b for a in 'ACGT' for b in 'ACGT' if a!=b]
    fig,axs=plt.subplots(1,3,figsize=(18,7));fig.subplots_adjust(top=.80,bottom=.27,wspace=.28)
    report=[]
    for ax,motif in zip(axs,['All AT/TA','AT','TA']):
        z=unique if motif=='All AT/TA' else events.loc[events.motif.eq(motif)].drop_duplicates('variant_key')
        n=z.substitution.value_counts().reindex(spectrum,fill_value=0)
        colors=['#0072B2' if substitution_class(v[0],v[2])=='Transition' else '#D55E00' for v in spectrum]
        ax.bar(spectrum,n,color=colors);ax.tick_params(axis='x',rotation=60);ax.set_ylabel('Unique observed SNV alleles')
        ti=int(z.substitution_class.eq('Transition').sum());tv=int(z.substitution_class.eq('Transversion').sum())
        ratio=f'{ti/tv:.2f}' if tv else 'undefined (no Tv)'
        ax.set_title(f'{motif}: Ti={ti:,}, Tv={tv:,}\nTi/Tv={ratio}')
        report.append(dict(motif=motif,transitions=ti,transversions=tv,Ti_Tv_ratio=ti/tv if tv else np.nan))
    fig.suptitle(f'Substitutions in {args.flank_bp}-bp reference flanks of AT and TA repeats',fontsize=17)
    fig.text(.08,.865,'Blue = transition (A↔G, C↔T); orange = transversion (including A↔T).',fontsize=11)
    fig.text(.08,.06,'PASS, non-reference diploid sample genotypes; REF checked against hg38. Symbolic alleles and indels are excluded.\nLabels show hg38 REF → sample ALT, not ancestral mutation direction or de novo mutations. Counts are not mutation rates.\nA site overlapping several loci is counted once per panel; a site may appear in both the AT and TA panels. No callable-base normalization.',fontsize=10)
    fig.savefig(out/'AT_TA_flank_transition_transversion_spectrum.png',dpi=180);plt.close(fig)
    pd.DataFrame(report).to_csv(out/'AT_TA_transition_transversion_summary.tsv',sep='\t',index=False)
    # Stratify by allele length; preserve the identity of each locus in denominators.
    groups=list(target.length_group.cat.categories)
    summary=[]
    fig,axs=plt.subplots(1,2,figsize=(14,7));fig.subplots_adjust(top=.83,bottom=.30,wspace=.25)
    for ax,motif in zip(axs,['AT','TA']):
        n_loci=[];x=np.arange(len(groups));z=target.loc[target.motif.eq(motif)]
        for j,kind in enumerate(['Transition','Transversion']):
            fractions=[]
            for group in groups:
                g=z.loc[z.length_group.eq(group)];n=len(g);positive=int(g[kind+'_N'].gt(0).sum())
                fractions.append(positive/n if n else np.nan)
                summary.append(dict(motif=motif,length_group=group,type=kind,N_loci=n,N_with_SNV=positive,fraction_with_recorded_SNV=positive/n if n else np.nan))
            ax.bar(x+(j-.5)*.36,fractions,width=.36,label=kind,color='#0072B2' if j==0 else '#D55E00')
        n_loci=[int(z.length_group.eq(group).sum()) for group in groups]
        ax.set_xticks(x,[f'{g}\nN={n:,}' for g,n in zip(groups,n_loci)],rotation=40,ha='right')
        from matplotlib.ticker import PercentFormatter
        ax.yaxis.set_major_formatter(PercentFormatter(1));ax.set_ylim(0,1);ax.set_title(motif);ax.legend();ax.set_ylabel('Loci with a recorded flank SNV (%)')
    fig.suptitle('Recorded flank substitutions by longer STR allele length',fontsize=17)
    fig.text(.08,.06,'A locus can have both substitution types, so the two percentages do not add to 100%.\nDenominator: all read-supported exact-motif loci in the length group, not only callable SNV sites.\nAbsence of a qualifying SNV record does not establish a reference genotype. Association does not establish causation.',fontsize=10)
    fig.savefig(out/'AT_TA_flank_substitutions_by_STR_length.png',dpi=180);plt.close(fig)
    pd.DataFrame(summary).to_csv(out/'AT_TA_substitutions_by_length.tsv',sep='\t',index=False)
    # Sequence-verified and AP-only subsets stay visibly distinct.
    categories=['No qualifying SNV recorded','Transition only','Transversion only','Both types']
    for name,mask,title in [('sequence_verified',seq_pure&target.reference_alternating_AT,'Sequence-verified alternating AT alleles'),('AP1_not_sequence_verified',target.TRGT_reported_AP1_both.eq(True)&target.both_alleles_alternating_AT.isna(),'TRGT-reported AP=1; allele sequences not verified')]:
        fig,grid=plt.subplots(2,2,figsize=(14,11),squeeze=False);axs=grid.ravel();fig.subplots_adjust(top=.84,bottom=.24,wspace=.3,hspace=.65)
        for ax,(motif,allele) in zip(axs,[('AT','shorter'),('TA','shorter'),('AT','longer'),('TA','longer')]):
            z=target.loc[mask&target.motif.eq(motif)]
            for i,cat in enumerate(categories):
                values=z.loc[z.SNV_group.eq(cat),allele+'_allele_minus_reference_bp'].to_numpy()
                if len(values):
                    ax.boxplot([values],positions=[i],widths=.55,showfliers=True)
                ax.text(i,.97,f'N={len(values):,}',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=9)
            ax.set_xlim(-.6,3.6);ax.set_xticks(range(4),categories,rotation=25,ha='right');ax.axhline(0,color='gray',ls='--');ax.set_title(motif)
            ax.set_ylabel(allele.capitalize()+' TRGT allele − hg38 repeat length (bp)')
            if z.empty: ax.text(.5,.5,'No eligible loci',transform=ax.transAxes,ha='center')
        fig.suptitle('STR length differences and nearby substitutions',fontsize=17)
        fig.text(.08,.87,title,fontsize=12)
        fig.text(.08,.06,'Box: middle 50%; line: median; whiskers: 1.5×IQR; points: outliers. One observation per locus.\nSubstitutions are in reference flanks, not counted from repeat-length changes. AT and TA retain separate catalog labels.\nNearby SNVs are not phased to the STR alleles; this does not show that a substitution caused a length change.\nNo qualifying SNV recorded includes potentially uncalled or filtered sites.',fontsize=10)
        fig.savefig(out/f'AT_TA_length_change_by_flank_SNV_{name}.png',dpi=180);plt.close(fig)
    audit.update(snv_vcf=str(args.snv_vcf),trgt_vcf=str(vcf) if vcf else None,AT_TA_loci=len(target),sequence_verified_pure_loci=int(seq_pure.sum()),sequence_verified_examples=len(selected),AP1_unverified_examples=len(ap_only),filter_policy='VCF FILTER=PASS; complete diploid nonreference GT; FT PASS or missing; no additional GQ/DP cutoff; REF checked against hg38',interpretation='Reference/sample differences, not ancestral direction, de novo status, causal effects, or mutation rates. No SNV/STR haplotype linkage asserted.')
    (out/'transition_transversion_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    fasta.f.close()
    print(f'AT/TA: {len(unique):,} unique qualifying flank SNV alleles; {len(selected):,} sequence-verified length-change examples; {len(ap_only):,} AP=1 examples without sequence verification.',flush=True)


def main():
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter
    root = Path.home()/'donglab'
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison', type=Path, default=root/'data/target_ALS/WGS_LR/repeat_comparison_gangSTR_vs_TRGT/comparison')
    parser.add_argument('--reference', type=Path, default=root/'references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa')
    parser.add_argument('--gtf', type=Path, default=root/'references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf', help='Matching hg38 GTF, optionally gzip-compressed; exact contig names must match.')
    parser.add_argument('--out', type=Path)
    parser.add_argument('--snv-vcf', type=Path, default=root/'data/target_ALS/WGS_LR/lr-wgs_vcf_genomic_NEUAD700YFB_NEUAD700YFB.SD-029-24-CBLL.3.g.vcf', help='Sample SNV VCF for the optional transition/transversion analysis.')
    parser.add_argument('--substitutions-only', action='store_true', help='Resume from saved locus_flanks_and_locations.tsv.gz; skip completed flank plots and GTF processing.')
    parser.add_argument('--skip-substitutions', action='store_true', help='Run only reference-flank and location analyses.')
    parser.add_argument('--trgt-sample', help='Exact TRGT VCF sample ID if multisample.')
    parser.add_argument('--vcf-sample', help='Exact VCF sample ID; required when there is more than one sample.')
    parser.add_argument('--trgt-vcf', type=Path, help='Optional original TRGT VCF to verify literal allele sequences.')
    parser.add_argument('--length-bins', default='0,20,50,100,200,500', help='Increasing lower boundaries in bp; final bin is open-ended.')
    parser.add_argument('--flank-bp', type=int, default=20)
    parser.add_argument('--top-motifs', type=int, default=10)
    parser.add_argument('--min-motif-loci', type=int, default=10)
    parser.add_argument('--scope', choices=['all-supported', 'gangstr-lower'], default='all-supported')
    a = parser.parse_args()
    if min(a.flank_bp, a.top_motifs, a.min_motif_loci) < 1:
        parser.error('Numeric settings must be positive.')
    if a.substitutions_only and a.skip_substitutions: parser.error('--substitutions-only and --skip-substitutions cannot be combined.')
    for needed in [a.reference,Path(str(a.reference)+'.fai')]+([] if a.substitutions_only else [a.gtf])+([] if a.skip_substitutions else [a.snv_vcf]):
        if not needed.is_file(): parser.error('Missing input: '+str(needed))
    edges=[int(v) for v in a.length_bins.split(',')]
    if not edges or edges[0]!=0 or any(y<=x for x,y in zip(edges,edges[1:])):
        parser.error('Length bins must start at zero and increase.')
    groups=[f'{x}–{y-1} bp' for x,y in zip(edges,edges[1:])]+[f'≥{edges[-1]} bp']
    out = a.out or a.comparison/'STR_flanking_transition_transversion_analysis'
    out.mkdir(parents=True, exist_ok=True)
    (out/'SUCCESS').unlink(missing_ok=True)
    print('Output directory: '+str(out),flush=True)
    if a.substitutions_only:
        saved=out/'locus_flanks_and_locations.tsv.gz'
        if not saved.is_file(): parser.error('Cannot resume: missing '+str(saved))
        print('Resuming substitution analysis from '+str(saved),flush=True)
        d=pd.read_csv(saved,sep='\t')
        if d.locus_id.duplicated().any(): raise ValueError('Duplicate saved locus IDs')
        if not (d.left_flank.str.len().eq(a.flank_bp)&d.right_flank.str.len().eq(a.flank_bp)).all():
            parser.error('Saved flank widths differ from --flank-bp; use the original width or run the full analysis.')
        d['length_group']=pd.cut(d.longer_TRGT_allele_bp,edges+[np.inf],labels=groups,right=False)
        analyze_substitutions(d,a,out)
        (out/'parameters_substitutions.json').write_text(json.dumps(vars(a),default=str,indent=2)+'\n')
        (out/'SUCCESS').write_text('complete\n')
        print('Substitution analysis completed: '+str(out),flush=True)
        return
    analysis = pd.read_csv(a.comparison/'analysis_loci.tsv.gz', sep='\t', usecols=['locus_id', 'chrom', 'start', 'end', 'motif', 'reference_bp'])
    evidence = pd.read_csv(a.comparison/'read_supported_STR_benchmark/locus_evidence.tsv.gz', sep='\t')
    if analysis.locus_id.duplicated().any() or evidence.locus_id.duplicated().any():
        raise ValueError('Duplicate locus IDs')
    passed = evidence.benchmark_pass.astype(str).str.lower().eq('true')
    evidence = evidence.loc[passed]
    if a.scope == 'gangstr-lower':
        evidence = evidence.loc[evidence.GangSTR_outcome.eq('Longer allele underestimated')]
    d = evidence.merge(analysis, on='locus_id', how='left', validate='one_to_one', indicator=True)
    if not d['_merge'].eq('both').all() or d.empty:
        raise ValueError('No selected loci or evidence IDs missing from analysis table.')
    d = d.drop(columns='_merge')
    expected_keys=d.chrom.astype(str)+':'+d.start.astype(str)+':'+d.end.astype(str)
    if not d.key.eq(expected_keys).all():
        raise ValueError('Evidence/analysis coordinate keys differ.')
    if not ((d.end-d.start+1).eq(d.reference_bp)).all():
        raise ValueError('Expected analysis coordinates to be 1-based inclusive.')
    d['start0'] = d.start.astype(int)-1
    d['end0'] = d.end.astype(int)
    d['longer_TRGT_allele_bp']=pd.to_numeric(d.longer_TRGT_allele_bp,errors='raise')
    if not (np.isfinite(d.longer_TRGT_allele_bp)&d.longer_TRGT_allele_bp.ge(0)).all():
        raise ValueError('Invalid TRGT allele length')
    d['length_group'] = pd.cut(d.longer_TRGT_allele_bp, edges+[np.inf], labels=groups, right=False)
    print('Reading reference flanks and GENCODE features for '+format(len(d),',')+' loci...',flush=True)
    fasta = IndexedFasta(a.reference)
    missing = set(d.chrom)-set(fasta.index)
    if missing:
        raise ValueError('Reference lacks contigs: '+str(sorted(missing)))
    annotations = load_gtf(a.gtf) if a.gtf else None
    if annotations and not set(d.chrom).intersection(annotations[1]):
        raise ValueError('No matching annotation contigs; verify hg38 build and chr naming.')
    results = []
    for r in d.itertuples():
        left, right = fasta.flanks(r.chrom, r.start0, r.end0, a.flank_bp)
        reference = fasta.fetch(r.chrom, r.start0, r.end0)
        before, after, aligned = phase_extension(reference, r.motif, left, right)
        feature, exonic, intronic = classify_feature(*annotations, r.chrom, r.start0, r.end0) if annotations else ('Not annotated', False, False)
        results.append((left, right, reference, before, after, aligned, feature, exonic, intronic))
    fasta.f.close()
    names=['left_flank','right_flank','reference_sequence','left_in_phase_bases','right_in_phase_bases','reference_matches_catalog_phase','genomic_feature','overlaps_exon','overlaps_intron']
    d[names] = pd.DataFrame(results, index=d.index)
    d['base_before_repeat']=d.left_flank.str[-1]
    d['base_after_repeat']=d.right_flank.str[0]
    d.to_csv(out/'locus_flanks_and_locations.tsv.gz', sep='\t', index=False)
    counts = pd.crosstab(d.motif, d.length_group).reindex(columns=groups, fill_value=0)
    counts.to_csv(out/'exact_motif_counts.tsv', sep='\t')
    totals=counts.sum(axis=1)
    motifs = totals.loc[totals>=a.min_motif_loci].sort_values(ascending=False).head(a.top_motifs).index.tolist()
    colors = {'A':'#2CA02C', 'C':'#0072B2', 'G':'#E69F00', 'T':'#B52B38'}
    records = []
    def composition(ax, seqs, positions, title, motif, group, side):
        bottom = np.zeros(len(positions))
        total = len(seqs)
        for base, color in colors.items():
            hits = np.array([sum(s[j] == base for s in seqs) for j in range(len(positions))])
            valid = np.array([sum(s[j] in 'ACGT' for s in seqs) for j in range(len(positions))])
            freq = np.divide(hits, valid, out=np.zeros(len(positions)), where=valid > 0)
            ax.bar(positions, freq, bottom=bottom, color=color, width=.9, label=base)
            bottom += freq
            records.extend(dict(motif=motif,group=group,side=side,position=int(pos),base=base,N_base=int(n),N_ACGT=int(v),fraction=float(f) if v else float('nan')) for pos,n,v,f in zip(positions,hits,valid,freq))
        ax.set_title(f'{title} | N={total:,}', fontsize=11)
        ax.set_ylim(0,1);ax.yaxis.set_major_formatter(PercentFormatter(1));ax.set_ylabel('Base frequency')
        if not total:
            ax.text(.5,.5,'No loci',transform=ax.transAxes,ha='center')
        ax.set_xticks(np.unique(np.r_[positions[::max(1, a.flank_bp//10)], positions[0], positions[-1]]))
    for rank, motif in enumerate([None]+motifs):
        z = d if motif is None else d.loc[d.motif.eq(motif)]
        label='All exact motifs pooled' if motif is None else motif
        fig, axs = plt.subplots(len(groups),2,figsize=(14,3*len(groups)+3),sharey=True,squeeze=False)
        fig.subplots_adjust(top=.87,bottom=.10,hspace=.60,wspace=.18)
        fig.suptitle(f'Reference flank composition: {label}',fontsize=18,y=.97)
        fig.text(.08,.945,f'TRGT read-supported loci | Groups use longer-allele length | Scope: {a.scope}',fontsize=10)
        fig.text(.08,.923,'Left flank → [catalog repeat] → right flank | hg38 forward strand; motif strings are not rotated or reverse-complemented.',fontsize=10)
        for i, group in enumerate(groups):
            g = z.loc[z.length_group.eq(group)]
            composition(axs[i,0],g.left_flank.tolist(),np.arange(-a.flank_bp,0),group+' STRs: left flank',label,group,'left')
            composition(axs[i,1],g.right_flank.tolist(),np.arange(1,a.flank_bp+1),group+' STRs: right flank',label,group,'right')
        axs[-1,0].set_xlabel('Position before catalog start (−1 = adjacent base)')
        axs[-1,1].set_xlabel('Position after catalog end (+1 = adjacent base)')
        fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper right',ncol=4,bbox_to_anchor=(.97,.98))
        fig.text(.08,.055,'Each locus contributes once. Ambiguous bases and contig-edge padding are excluded position by position.\nFlanks use reference catalog boundaries, not inferred sample insertion boundaries; left/right are not gene-strand upstream/downstream.\nLength groups are not matched for GC, coverage, or genomic features; pooled profiles also mix motif frequencies. These plots describe composition, not enrichment.',fontsize=9)
        fig.savefig(out/('flank_composition_all_motifs.png' if motif is None else f'flank_composition_motif_{rank:02d}.png'),dpi=180);plt.close(fig)
    pd.DataFrame(records).to_csv(out/'flank_base_frequencies.tsv',sep='\t',index=False)
    if motifs:
        selected = counts.loc[motifs]
        fig, ax = plt.subplots(figsize=(max(9,len(motifs)*.7),6))
        x=np.arange(len(motifs))
        width=.8/len(groups)
        for j,group in enumerate(groups):
            ax.bar(x+(j-(len(groups)-1)/2)*width,selected[group],width=width,label=group)
        ax.set_xticks(x,motifs,rotation=45,ha='right');ax.set_ylabel('Read-supported loci (log scale)');ax.set_yscale('log')
        ax.set_title(f'Exact motif counts | Top {len(motifs)} by total locus count');ax.legend()
        fig.text(.08,.02,'Length groups use the longer TRGT allele. Exact catalog strings; no rotations or reverse complements grouped.',fontsize=9)
        fig.tight_layout(rect=(0,.08,1,1));fig.savefig(out/'exact_motif_counts.png',dpi=180);plt.close(fig)
    if annotations:
        order=['Exonic','Intronic','Exonic + intronic','Other genic','Intergenic','Unannotated contig']
        tab=pd.crosstab(d.genomic_feature,d.length_group).reindex(index=order,columns=groups,fill_value=0)
        tab.to_csv(out/'genomic_feature_counts.tsv',sep='\t')
        nrows=(len(groups)+1)//2
        fig,grid=plt.subplots(nrows,2,figsize=(15,nrows*4+2),squeeze=False);axs=grid.ravel()
        fig.subplots_adjust(bottom=.17,top=.91,wspace=.3,hspace=.9)
        for i,group in enumerate(groups):
            n=int(tab[group].sum());fraction=tab[group]/n if n else tab[group]*0
            bars=axs[i].bar(order,fraction,color=plt.get_cmap('tab10')(i%10))
            axs[i].set_title(f'{group} STRs | N={n:,}');axs[i].tick_params(axis='x',rotation=45)
            axs[i].set_ylim(0,1);axs[i].yaxis.set_major_formatter(PercentFormatter(1))
            axs[i].set_ylabel('Fraction of loci')
            for bar,count in zip(bars,tab[group]):
                axs[i].text(bar.get_x()+bar.get_width()/2,bar.get_height()+.015,f'{count:,}',ha='center',fontsize=9)
        for unused in axs[len(groups):]: unused.set_visible(False)
        fig.suptitle('Genomic locations of read-supported STR loci',fontsize=18)
        fig.text(.08,.045,'Any overlap with the reference STR interval; transcript isoforms are combined. Exonic includes coding and noncoding exons.\nExonic + intronic means both feature types overlap (potentially across different transcripts); other genic overlaps a gene/transcript only.\nIntergenic requires an annotated contig. Locations describe hg38 catalog intervals, not the span of expanded sample alleles.\nThese are observed proportions in the selected catalog, not enrichment relative to the genome.',fontsize=9)
        fig.savefig(out/'genomic_feature_distribution.png',dpi=180);plt.close(fig)
    print('Reference flank and genomic-location outputs finished.',flush=True)
    if not a.skip_substitutions:
        analyze_substitutions(d,a,out)
    (out/'METHODS.txt').write_text('Reference flanks in forward genomic orientation, anchored to unmodified shared catalog intervals. No sample-specific insertion boundary inference. Exact motif strings remain separate. Flank tables retain positions and denominators. In-phase adjacent bases flag potentially extendable pure-reference repeats within the inspected window only; they do not establish shifted sample boundaries. Gene annotation uses any interval overlap and all transcript models. Sample substitutions use PASS SNV ALT alleles carried in a complete diploid GT from the supplied gVCF; reference blocks, indels and symbolic alleles are excluded. No added GQ/DP threshold; those fields are saved. Counts are observed reference/sample differences, not de novo mutation rates or ancestral directions. Allele dosage does not weight counts. No recorded SNV is not proof of absence without callable-site analysis. Sequence-verified alternating AT/TA alleles are analyzed separately from unverified AP=1 calls; no sample haplotype linkage or causal effects of flank substitutions are claimed. TRGT length support reuses genotyping reads and is not independent accuracy validation.\n')
    (out/'parameters.json').write_text(json.dumps(vars(a),default=str,indent=2)+'\n')
    (out/'SUCCESS').write_text('complete\n' if annotations else 'flanks_complete; gene_annotation_not_run\n')
    print('Saved:',out)
    if not annotations:
        print('Gene-location plots require --gtf /path/to/matching_hg38_annotation.gtf.gz')
    if not motifs:
        print('No motifs meet --min-motif-loci; locus and motif-count tables were saved.')


if __name__ == '__main__':
    main()

PYTHON_STR_FLANKS
