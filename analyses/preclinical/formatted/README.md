# Formatted cell line features and metadata
The contents of this directory are generated a few scripts: `01.organize-samples.ipynb` and `02.process-cell-line-molecular-features.ipynb`. After running these scripts, the following files should be present in this directory:
- `ccle.copy-numbers.txt`
- `ccle.variants.txt`
- `cell-line-names.raw.txt`
- `cell-line-names.formatted.txt`
- `cell-lines.summary.txt`
- `sanger-fusions.txt`
- `sanger.gdsc.pairwise-sensitive.txt`
- `sanger.gdsc.txt`

`cell-line-names.formatted.txt` is a file that underwent some manual checking from `cell-line-names.raw.txt`. The formatted version is committed in this Github repo, please regenerate the other files to keep this repository (relatively) lean.
