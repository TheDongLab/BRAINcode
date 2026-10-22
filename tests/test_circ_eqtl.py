import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'modules'))
from circ_eqtl import prepare_circ

class CircPreparation(unittest.TestCase):
    def test_final_cohort_no_support_filter(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            donors = [f'D{i}' for i in range(16)]
            ids = ['chr1:100-200:+', 'chr1:300-400:-', 'chr1:500-600:+', 'chr1:700-800:+']
            # Keep a singleton below five reads; drop zero, constant and a row
            # whose only variation is in the donor with missing covariates.
            rpm = pd.DataFrame([[0.01]+[0]*15, [0]*16, [1]*16, [0]*15+[1]], index=ids, columns=donors)
            snps = pd.DataFrame([[0]*5+[1]*5+[2]*6], index=['rs1'], columns=donors)
            cov = pd.DataFrame([range(16)], index=['age'], columns=donors, dtype=float)
            cov.loc['age','D15'] = np.nan
            source = p/'circ.txt'
            pd.DataFrame({'chrom':['chr1'], 'start':[100], 'end':[200], 'strand':['+'], 'geneName':['HOST'], 'isoformName':['TX']}).to_csv(source,sep='\t',index=False)
            pd.DataFrame({'sample_id':donors,'circ_file':[str(source)]*16}).set_index('sample_id').to_csv(p/'provenance.tsv',sep='\t')
            z,g,c = prepare_circ(rpm,snps,cov,p,'Test',p/'provenance.tsv',dict(zip(donors,donors)))
            self.assertEqual(list(z.index), ids[:1])
            self.assertEqual(list(z.columns), donors[:15])
            self.assertEqual(list(z.columns), list(g.columns))
            self.assertEqual(list(z.columns), list(c.columns))
            self.assertTrue(np.isfinite(z.to_numpy()).all())
            self.assertGreater(z.iloc[0,0], z.iloc[0,1])
            self.assertEqual(z.iloc[0,1], z.iloc[0,2])
            loc=pd.read_csv(p/'gene_location.txt',sep='\t')
            self.assertEqual(loc.left.iloc[0],101)
            self.assertEqual(loc.right.iloc[0],200)
            ann=pd.read_csv(p/'circ_annotation.tsv',sep='\t')
            self.assertEqual(ann.host_gene.iloc[0],'HOST')
            self.assertEqual(pd.read_csv(p/'circ_rpm_Test.txt',sep='\t').iloc[0,1],0.01)

if __name__ == '__main__': unittest.main()
