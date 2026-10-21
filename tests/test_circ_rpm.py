"""Synthetic checks for the Python embedded in the circ matrix builder."""
import ast
import contextlib
import io
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from collections import defaultdict

SCRIPT = Path(__file__).resolve().parents[1] / 'modules/build_circ_matrix_for_QTL.sh'
CODE = SCRIPT.read_text().split("python3 - <<'EOF'\n", 1)[1].split('\nEOF', 1)[0]
NS = dict(defaultdict=defaultdict)
for node in ast.parse(CODE).body:
    if isinstance(node, ast.FunctionDef) and node.name in ('mapped_fragment_stats', 'audit_star_junctions'):
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(SCRIPT), 'exec'), NS)


def mapping(path, unique=15000000, multi=5000000):
    (path / 'STAR.Log.final.out').write_text(
        'Number of input reads | 25000000\n'
        f'Uniquely mapped reads number | {unique}\n'
        f'Number of reads mapped to multiple loci | {multi}\n')


class CircRPMTests(unittest.TestCase):
    def test_mapping_denominator(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d); mapping(p)
            stats = NS['mapped_fragment_stats'](p)
            self.assertEqual(stats['mapped_fragments'], 20000000)
            self.assertAlmostEqual(stats['rpm_multiplier'] * 10, .5)
            mapping(p, 0, 0)
            with self.assertRaises(ValueError): NS['mapped_fragment_stats'](p)
            (p / 'STAR.Log.final.out').unlink()
            with self.assertRaises(ValueError): NS['mapped_fragment_stats'](p)

    def test_junction_audit(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            row = 'chr1\t201\t+\tchr1\t100\t+\t1\t0\t0\tfrag1\t90\t10M10S\t190\t10S10M\n'
            (p / 'STAR.Chimeric.out.junction').write_text(row)
            (p / 'back_spliced_junction.bed').write_text('chr1\t100\t200\tFUSIONJUNC_0/1\t0\t+\n')
            self.assertEqual(NS['audit_star_junctions'](p)['junction_audit'],
                             'unique_fragment_support_matches_parsed_bed')
            (p / 'STAR.Chimeric.out.junction').write_text(row + row)
            with self.assertRaisesRegex(ValueError, 'repeated fragment'): NS['audit_star_junctions'](p)
            (p / 'STAR.Chimeric.out.junction').write_text(row.replace('201', '301'))
            with self.assertRaisesRegex(ValueError, 'disagree'): NS['audit_star_junctions'](p)

    def test_full_builder_exports_unfiltered_counts_and_rpm(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            # Reverse metadata order to exercise name-aligned normalization.
            pd.DataFrame(dict(externalsampleid=['B', 'A'], externalsubjectid=['P2', 'P1'],
                              tissue=['Motor Cortex'] * 2)).to_csv(root / 'targetALS_rnaseq_metadata.csv', index=False)
            for sid, unique, counts, fractions in [('A', 15000000, [10, 2], [.1, .001]),
                                                   ('B', 35000000, [20, 0], [.2, 0])]:
                p = root / 'Motor_Cortex/RNAseq/Processed' / sid
                p.mkdir(parents=True); mapping(p, unique, 5000000)
                # Correct total input for sample B.
                f = p / 'STAR.Log.final.out'
                f.write_text(f.read_text().replace('25000000', '50000000'))
                pd.DataFrame(dict(chrom=['chr1'] * 2, start=[100, 300], end=[200, 400],
                                  strand=['+'] * 2, readNumber=counts, circ_reads=counts,
                                  circ_percent=fractions)).to_csv(p / 'circularRNA_known_circ_percentage.txt', sep='\t', index=False)
            code = CODE.replace('BASE = Path("/home/zw529/donglab/data/target_ALS")', f'BASE = Path({str(root)!r})')
            with contextlib.redirect_stdout(io.StringIO()): exec(compile(code, str(SCRIPT), 'exec'), {})
            count = pd.read_csv(root / 'QTL/circ_counts_matrix.txt', sep='\t', index_col=0)
            rpm = pd.read_csv(root / 'QTL/circ_rpm_matrix.txt', sep='\t', index_col=0)
            pct = pd.read_csv(root / 'QTL/circ_matrix.txt', sep='\t', index_col=0)
            meta = pd.read_csv(root / 'QTL/circ_rpm_sample_metadata.tsv', sep='\t', index_col=0)
            self.assertEqual(list(count.columns), ['B', 'A'])
            self.assertEqual(count.shape, (2, 2))
            self.assertEqual(pct.shape, (1, 2))  # low fraction SD does not remove counts/RPM
            self.assertEqual(count.loc['chr1:100-200:+', 'A'], 10)
            self.assertAlmostEqual(rpm.loc['chr1:100-200:+', 'A'], .5)
            self.assertAlmostEqual(rpm.loc['chr1:100-200:+', 'B'], .5)
            self.assertEqual(rpm.loc['chr1:300-400:+', 'B'], 0)
            self.assertEqual(meta.loc['B', 'mapped_fragments'], 40000000)
            self.assertEqual(meta.loc['A', 'junction_audit'], 'unavailable_no_chimeric_file')

if __name__ == '__main__': unittest.main()
