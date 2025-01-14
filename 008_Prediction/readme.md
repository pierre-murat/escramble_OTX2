# Readme

The python script `borzoi_helpers.py` and other borzoi-related (`targets_gtex.txt`, `params_pred.json`) files were downloaded from https://github.com/calico/borzoi, which was distributed under Apache 2.0 (https://www.apache.org/licenses/LICENSE-2.0).

Before running `predict_deletions.py`, make sure to 1) follow borzoi installation instructions and 2) download borzoi model folds and place in `examples/saved_models/f0/model0_best.h5` etc.
Further, in `predict_deletions.py`, change the line 
```
borzoi_dir =  # SET THIS TO DIRECTORY CONTAINING 'examples/saved_models/f0/model0_best.h5' ETC
```
to point to the right directory.

Furthermore, download https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz and use `gunzip` to unzip, `predict_deletions.py` uses this to get the sequence.

After running `predict_deletions.py`, run the analysis `process_prediction_data.Rmd` to generate plots.
