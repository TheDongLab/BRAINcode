"""CircRNA-specific preparation shared by the optional eQTL circ mode."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm


def prepare_circ(rpm, snps, cov, outdir, tissue, provenance_file, donor_to_sample):
    """Keep variable RPM rows; no BSJ-read or prevalence threshold is applied."""
    out = Path(outdir)
    if rpm.index.has_duplicates or rpm.columns.has_duplicates:
        raise ValueError('Duplicated circRNA or donor IDs')
    cov = cov.apply(pd.to_numeric, errors='coerce')
    donors = [d for d in rpm.columns if d in cov.columns and np.isfinite(cov[d]).all()]
    if len(donors) < 3:
        raise ValueError('Fewer than three donors with complete covariates')
    rpm, snps, cov = rpm.loc[:, donors], snps.loc[:, donors], cov.loc[:, donors]
    values = rpm.to_numpy(float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('RPM values must be finite and nonnegative')
    # Exact variability, not an abundance/coverage filter.
    keep = rpm.max(axis=1) > rpm.min(axis=1)
    print(f'Circ mode: {len(donors)} complete donors; {int(keep.sum())}/{len(keep)} variable circRNAs; no read/prevalence filter')
    rpm = rpm.loc[keep].copy()
    if rpm.empty:
        raise ValueError('No variable circRNAs in the final donor cohort')
    transformed = pd.DataFrame(norm.ppf((rpm.rank(axis=1, method='average') - 0.5) / len(donors)), index=rpm.index, columns=donors)
    # Covariate intersection can change the cohort: repeat variant QC on those donors.
    snps = snps.apply(pd.to_numeric, errors='coerce')
    called = snps.notna().sum(axis=1)
    af = snps.sum(axis=1) / (2 * called)
    groups = snps.apply(lambda r: len(r.value_counts()) >= 2 and r.value_counts().min() >= 5, axis=1)
    snps = snps.loc[(called > 0) & (np.minimum(af, 1-af) >= 0.05) & groups]
    if snps.empty:
        raise ValueError('No SNPs survive QC in the final cohort')
    # CIRCexplorer2 BED start is zero-based; SNP positions are one-based.
    coords = rpm.index.to_series().str.extract(r'^(?P<chr>[^:]+):(?P<start>\d+)-(?P<end>\d+):(?P<strand>[+-])$')
    if coords.isna().any().any():
        raise ValueError('Invalid circRNA coordinate ID')
    left = coords.start.astype('int64') + 1
    right = coords.end.astype('int64')
    if (left > right).any():
        raise ValueError('Invalid circRNA interval')
    loc = pd.DataFrame({'geneid': rpm.index, 'chr': coords.chr.to_numpy(), 'left': left.to_numpy(), 'right': right.to_numpy()})
    loc.to_csv(out/'gene_location.txt', sep='\t', index=False)
    # Host annotations come from the same per-sample CIRCexplorer2 files used by the builder.
    provenance = pd.read_csv(provenance_file, sep='\t', index_col=0)
    genes, transcripts = {}, {}
    selected = set(rpm.index)
    for donor in donors:
        sample = donor_to_sample[donor]
        path = provenance.loc[sample, 'circ_file']
        tab = pd.read_csv(path, sep='\t', usecols=['chrom','start','end','strand','geneName','isoformName'])
        ids = tab.chrom.astype(str)+':'+tab.start.astype(str)+'-'+tab.end.astype(str)+':'+tab.strand.astype(str)
        for cid, gene, transcript in zip(ids, tab.geneName, tab.isoformName):
            if cid not in selected:
                continue
            for value, target in [(gene, genes), (transcript, transcripts)]:
                if pd.notna(value) and str(value) not in ('', '.', 'NA'):
                    target.setdefault(cid, set()).add(str(value))
    pd.DataFrame({'geneid': rpm.index,
                  'host_gene': [';'.join(sorted(genes.get(c, []))) for c in rpm.index],
                  'host_transcript': [';'.join(sorted(transcripts.get(c, []))) for c in rpm.index]}
                 ).to_csv(out/'circ_annotation.tsv', sep='\t', index=False)
    rpm.to_csv(out/f'circ_rpm_{tissue}.txt', sep='\t', index_label='circ_id')
    pd.DataFrame({'donor_id': donors, 'sample_id': [donor_to_sample[d] for d in donors]}).to_csv(out/'selected_samples.tsv', sep='\t', index=False)
    return transformed, snps, cov
