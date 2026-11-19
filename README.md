# Reproduction of Li 2025 & Liu 2026

This repository documents my attempts to reproduce the interfacial chemistry reported by [Li et al. (2025)](https://doi.org/10.1038/s41467-025-59712-3) and [Liu et al. (2026)](https://doi.org/10.5194/acp-26-9357-2026). I have not yet reproduced the target interfacial reaction behavior or established converged reproductions of the reported free-energy barriers. **This outcome does not demonstrate that either paper is incorrect.**

The calculations do not yet provide a fully matched comparison with the original studies. In particular, the Liu reactant starting structures and preparation procedure have not been matched, some implementation details remain unverified against the original inputs, and sampling and free-energy convergence are incomplete. These gaps may affect the comparison, but their individual effects have not been established. The inputs, selected trajectories and analysis scripts are shared here to make those questions explicit and support further work.

| Study | System examined | Outcome of these calculations |
|---|---|---|
| Li 2025 | SO₃²⁻ + HCHO, 96 waters, BLYP-D3(BJ) | Unbiased interfacial formation from separated reactants was not reproduced; no converged reproduction of the reported barrier was established. |
| Liu 2026 | HOSO₂⁻ + HCHO + K⁺, 128 waters, PBE-D3 | Unbiased interfacial formation and the reported interfacial barrier were not reproduced. Sampling is paused. |

Short C–S states were obtained under metadynamics or restraints, and some prepared states retained short distances after bias removal. The progress reports distinguish these observations from formation under unbiased conditions and describe the remaining comparisons:

- [Li 2025: results, model and remaining comparisons](docs/PROGRESS.md)
- [Liu 2026: results, model and reasons the reproduction remained incomplete](docs/LIU2026_PROGRESS.md)
- [Additional results available on request](docs/ADDITIONAL_RESULTS.md): a detailed catalogue of material left out of the downloads.

## Downloads and reanalysis

Install Python 3.12 and NumPy:

```bash
git clone https://github.com/KejunLiuGitHub/reproduction_of_li_2025_and_liu_2026.git
cd reproduction_of_li_2025_and_liu_2026
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

For Li 2025:

```bash
python tools/fetch_data.py --dataset li2025
python tools/reproduce_case.py --case natural_original40 --out results/li_natural_original40
```

For Liu 2026:

```bash
python tools/fetch_data.py --dataset liu2026
python tools/reproduce_liu_case.py --case frame_ic_1 --out results/liu_frame_ic_1
```

The [Li data release](https://github.com/KejunLiuGitHub/reproduction_of_li_2025_and_liu_2026/releases/tag/v2026.11.19) contains 46 saved trajectories, including smoke runs, in a 122.5 MB archive. The [Liu data release](https://github.com/KejunLiuGitHub/reproduction_of_li_2025_and_liu_2026/releases/tag/liu2026-v2026.11.19) is 67.5 MB compressed (199.1 MB unpacked) and contains selected unbiased and guided trajectories, the four-arm proton-location comparison, and three complete archived metadynamics CV/bias histories. Full metadynamics XYZ histories and larger supporting collections are listed in the on-request catalogue.

Both download commands check archive and individual-file SHA-256 values. The two studies have separate atom mappings, data manifests and analyses. See [Li reproduction instructions](docs/REPRODUCE.md) and [Liu reproduction instructions](docs/LIU2026_REPRODUCE.md).

## What is included

`campaign/` holds the scientific analysis code and selected input records. `data/` contains case manifests; large files are downloaded from Releases. `metadata/` records checksums and provenance. `docs/` describes the calculations, results and omitted material.

The Li snapshot includes two running controls through 2.930 and 3.100 ps on 19 November 2026. The Liu record includes five complete 3 ps first segments from an unfinished 5 ps protocol. These coverage limits are retained in the case descriptions. Larger historical files can be requested using the identifiers in [the catalogue](docs/ADDITIONAL_RESULTS.md).

## Citation and reuse

Cite the original papers when discussing their methods or conclusions, and this repository when using these calculations. Li starting coordinates derive from the authors' [Zenodo dataset](https://zenodo.org/records/15253851). Article PDFs and the author workbook are obtained from their original sources. [Reuse notes](REUSE.md) describe attribution.
