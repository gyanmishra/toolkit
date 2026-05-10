import sys
import tempfile
from pathlib import Path
import types
import unittest

# Stub pandas for import-time dependency
sys.modules.setdefault("pandas", types.SimpleNamespace(DataFrame=None))

from GenerateQualimap_report import parse_report


class TestQualimapParsing(unittest.TestCase):
    def test_parse_report_key_values(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "rnaseq_qc_results.txt"
            p.write_text("A = 1\nnotkv\nB=2\n", encoding="utf-8")
            result = parse_report(p)
            self.assertEqual(result["A"], "1")
            self.assertEqual(result["B"], "2")
            self.assertNotIn("notkv", result)


if __name__ == "__main__":
    unittest.main()
