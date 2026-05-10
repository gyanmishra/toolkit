import sys
import tempfile
from pathlib import Path
import types
import unittest
from unittest.mock import patch

# Stub optional third-party modules for import-time dependency
sys.modules.setdefault("GEOparse", types.SimpleNamespace(GEO=types.SimpleNamespace(GSE=object)))
sys.modules.setdefault("pandas", types.SimpleNamespace(DataFrame=None))
sys.modules.setdefault("Bio", types.SimpleNamespace(Entrez=types.SimpleNamespace()))

from download_GEOmetadata import maybe_gzip, download_sra


class TestDownloadGeoMetadata(unittest.TestCase):
    @patch("download_GEOmetadata.subprocess.run")
    def test_maybe_gzip_runs_when_plain_file_exists(self, mock_run):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "x.fastq"
            p.write_text("abc", encoding="utf-8")
            gz = maybe_gzip(p)
            self.assertEqual(gz.name, "x.fastq.gz")
            mock_run.assert_called_once()

    @patch("download_GEOmetadata.subprocess.run")
    def test_download_sra_skips_when_gz_exists(self, mock_run):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            (d / "SRR1.fastq.gz").write_text("gz", encoding="utf-8")
            download_sra("SRR1", d)
            mock_run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
