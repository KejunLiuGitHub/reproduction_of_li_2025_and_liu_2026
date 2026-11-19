# Li 2025: interfacial C–S formation

Simulation inputs, selected trajectories and analysis code for the SO₃²⁻ + HCHO pathway in [Li et al., *Nature Communications* **16**, 5187 (2025)](https://doi.org/10.1038/s41467-025-59712-3).

We obtained short C–S states in interfacial metadynamics and followed two of them for 5 ps after removing the bias. **Unbiased interfacial C–S formation starting from separated reactants has not been reproduced.** This repository records the calculations and provides a starting point for continuing the study.

## Results so far

| Calculation | Coverage | Main observation |
|---|---|---|
| Original interfacial metaD chain | 50 ps | Visits short C–S states; returns to longer distances in the final segment |
| Separate reactant-started metaD chain | 20 ps | Samples a second set of short-distance configurations |
| Two selected-state releases | 5 ps each | Mean S–C distances of 1.913 and 1.890 Å after bias removal |
| Four local-water orientation tests | 5 ps each | No C–S formation geometry observed from the long-distance starts |
| Two unrotated controls | 2.930 and 3.100 ps in this snapshot | Still running toward 5 ps |

The selected-state releases begin with short C–S distances. The water-orientation tests begin with separated reactants. Their different starting conditions are part of the comparison. [Progress and interpretation](docs/PROGRESS.md) describes the complete results, including the solvent depth, earlier unsuccessful releases and method differences.

## Get the data and run an analysis

```bash
git clone https://github.com/KejunLiuGitHub/li2025-interface-cs-reproduction.git
cd li2025-interface-cs-reproduction
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python tools/fetch_data.py
python tools/reproduce_case.py --case natural_original40 --out results/natural_original40
```

The [v2026.11.19 data release](https://github.com/KejunLiuGitHub/li2025-interface-cs-reproduction/releases/tag/v2026.11.19) is **122.5 MB compressed** (390.1 MB unpacked). It contains 46 saved trajectories, including short smoke runs, and associated scientific files. `fetch_data.py` checks the archive and every unpacked file against SHA-256 records. [Data layout and selection](docs/DATA.md) explains which files were included.

Analysis uses Python and NumPy; it does not require cluster access. [Reproduction instructions](docs/REPRODUCE.md) also show how to prepare the original-control inputs for a new CP2K run. The archived inputs retain their original settings and file paths; the preparation tool writes a separate portable copy.

## Repository layout

- `campaign/`: original scientific scripts, inputs and parameter records, preserved from the source commit.
- `tools/`: data download, trajectory reanalysis and portable input preparation.
- `data/DATASETS.json`: trajectory names, boxes, bias history and snapshot coverage.
- `metadata/`: source and data checksums, selection and release manifests.
- `docs/`: results, model details and reproduction instructions.

Large data files are attached to the release and are not part of Git history. The two running controls are intermediate snapshots; all other statuses should be read from their run records. The source project continues to retain full verbose logs and omitted restart material.

## Source data and citation

Starting structures are derived from the authors’ [Zenodo dataset, DOI 10.5281/zenodo.15253851](https://zenodo.org/records/15253851). The original workbook and article PDFs are not redistributed here. Cite the paper and dataset when using these structures. See [model and provenance](docs/METHODS.md) and [reuse notes](REUSE.md).

## 中文概述

目前已在界面metaD中获得短C–S态，两条撤偏轨迹随后运行5 ps，S–C均值为1.913和1.890 Å。从分离反应物出发的无偏界面成键尚未复现：四条旋水试验各完成5 ps均未出现目标几何，两条原样对照正运行至2.930和3.100 ps。接下来完成对照，结合水网络、溶质取向和界面深度比较完整轨迹。
