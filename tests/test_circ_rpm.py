"""Synthetic checks for the Python embedded in the circ matrix builder."""
import ast
import contextlib
import io
from pathlib import Path
import tempfile
import re
import unittest
import numpy as np
import pandas as pd
from collections import defaultdict

SCRIPT = Path(__file__).resolve().parents[1] / 'modules/build_circ_matrix_for_QTL.sh'
CODE = SCRIPT.read_text().split("python3 - <<'EOF'\n", 1)[1].split('\nEOF', 1)[0]
NS = dict(defaultdict=defaultdict, re=re)
for node in ast.parse(CODE).body:
    if isinstance(node, ast.FunctionDef) and node.name in ('mapped_read_stats', 'audit_star_junctions', 'norm_id', 'find_circ'):
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(SCRIPT), 'exec'), NS)


def mapping(path, passed=15000000, failed=5000000):
    (path / 'bam.stat').write_text(
        '99999999\n'
        '99999999 + 0 in total (QC-passed reads + QC-failed reads)\n'
        '60000000 + 5000000 primary\n'
        '90000000 + 0 mapped (90.00% : N/A)\n'
        f'{passed} + {failed} primary mapped (80.00% : N/A)\n')


class CircRPMTests(unittest.TestCase):
    def test_first_matching_path(self):
        first = Path('/data/Cervical_Spinal_Cord/RNAseq/Processed/CGND_HRA_03699/circularRNA_known_circ_percentage.txt')
        second = Path('/data/Cervical_spinal_cord/RNAseq/Processed/CGND_HRA_03699/circularRNA_known_circ_percentage.txt')
        NS['all_circ'] = [first, second]
        self.assertEqual(NS['find_circ']('CGND-HRA-03699'), first)
        NS['all_circ'] = [second, first]
        self.assertEqual(NS['find_circ']('CGND-HRA-03699'), second)
        self.assertIsNone(NS['find_circ']('absent'))

    def test_mapping_denominator(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d); mapping(p)
            stats = NS['mapped_read_stats'](p)
            self.assertEqual(stats['primary_mapped_reads'], 20000000)
            self.assertAlmostEqual(stats['rpm_multiplier'] * 10, .5)
            self.assertEqual(stats['primary_mapped_qc_fail'], 5000000)
            # Never accept the generic mapped field or a duplicate denominator line.
            original = (p / 'bam.stat').read_text()
            for broken in [original.replace('primary mapped', 'other mapped'),
                           original + '1 + 0 primary mapped (1% : N/A)\n',
                           original.replace('60000000 + 5000000 primary', '1 + 0 primary')]:
                (p / 'bam.stat').write_text(broken)
                with self.assertRaises(ValueError): NS['mapped_read_stats'](p)
            mapping(p, 0, 0)
            with self.assertRaises(ValueError): NS['mapped_read_stats'](p)
            (p / 'bam.stat').unlink()
            with self.assertRaises(ValueError): NS['mapped_read_stats'](p)

    def test_junction_audit(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            row = 'chr1\t201\t+\tchr1\t100\t+\t1\t0\t0\tfrag1\t90\t10M10S\t190\t10S10M\n'
            (p / 'STAR.Chimeric.out.junction').write_text(row)
            (p / 'back_spliced_junction.bed').write_text('chr1\t100\t200\tFUSIONJUNC_0/1\t0\t+\n')
            self.assertEqual(NS['audit_star_junctions'](p)['junction_audit'],
                             'junction_records_match_parsed_bed')
            (p / 'STAR.Chimeric.out.junction').write_text(row + row)
            (p / 'back_spliced_junction.bed').write_text('chr1\t100\t200\tFUSIONJUNC_0/2\t0\t+\n')
            self.assertEqual(NS['audit_star_junctions'](p)['repeated_read_name_junction_records'], 1)
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
            self.assertEqual(meta.loc['B', 'primary_mapped_reads'], 40000000)
            self.assertEqual(meta.loc['A', 'junction_audit'], 'unavailable_no_chimeric_file')

if __name__ == '__main__': unittest.main()
