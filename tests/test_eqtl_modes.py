"""Check legacy argument compatibility without running cluster jobs."""
import subprocess
import unittest
from pathlib import Path

MODULES = Path(__file__).resolve().parents[1]/'modules'
class ModeArguments(unittest.TestCase):
    def test_arguments(self):
        s=(MODULES/'run_eQTL.sh').read_text()
        parser=s[s.index('TISSUE="$1"'):s.index('TISSUE_DIR=')]
        for args,expected in [(['Cerebellum'],'gene standard  eQTL'),(['Cerebellum','circ'],'circ standard  circ_eQTL'),(['Cerebellum','interaction'],'gene interaction  eQTL'),(['Cerebellum','interaction','Male'],'gene interaction Male eQTL'),(['Cerebellum','circ','interaction'],'circ interaction  circ_eQTL'),(['Cerebellum','interaction','Male','circ'],'circ interaction Male circ_eQTL')]:
            result=subprocess.check_output(['bash','-c','set -eu\n'+parser+'\nprintf "%s %s %s %s" "$MODE" "$RUN_TYPE" "$SUB_DIR" "$ANALYSIS"','test',*args],text=True)
            self.assertEqual(result,expected)
    def test_syntax(self):
        for n in ('prep_eQTL.sh','run_eQTL.sh'):
            subprocess.run(['bash','-n',str(MODULES/n)],check=True)
        s=(MODULES/'prep_eQTL.sh').read_text().split('python3 << EOF\n',1)[1].split('\nEOF',1)[0]
        compile(s,'embedded preparation Python','exec')

if __name__=='__main__': unittest.main()
