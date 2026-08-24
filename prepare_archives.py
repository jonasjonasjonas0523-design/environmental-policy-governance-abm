"""Prepare archived results for local reproduction without altering the immutable ZIP files."""
from pathlib import Path
import zipfile, shutil

ROOT = Path(__file__).resolve().parent
R = ROOT / "results"

with zipfile.ZipFile(R / "v12_results_bundle.zip") as z:
    for member in z.namelist():
        if member.startswith("results_v12/"):
            z.extract(member, ROOT)

with zipfile.ZipFile(R / "v11_stage1_archive.zip") as z:
    for member in z.namelist():
        if member.startswith("project/results/") and not member.endswith("/"):
            rel = Path(member).relative_to("project")
            out = ROOT / rel
            out.parent.mkdir(parents=True, exist_ok=True)
            with z.open(member) as src, open(out, "wb") as dst:
                shutil.copyfileobj(src, dst)

print("Prepared working directories: results/ and results_v12/")
