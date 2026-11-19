# Model and provenance

The calculations follow the SO₃²⁻ + HCHO pathway of [Li et al. (2025)](https://doi.org/10.1038/s41467-025-59712-3). Coordinates originate from the authors' Figure 6b source workbook, available at [Zenodo, DOI 10.5281/zenodo.15253851](https://zenodo.org/records/15253851).

## Molecular system

The Li system contains 296 atoms: 96 waters, formaldehyde and sulfite. The atom order is H,H,O for atoms 1–288, followed by C289/O290/H291/H292 and S293/O294/O295/O296. The total charge is −2, with singlet spin and no explicit counterions. The interface cell is 15 × 15 × 30 Å; bulk calculations use a cubic 14.2296 Å cell.

## Electronic structure and dynamics

The recorded software is CP2K 9.1 (`cp2k.psmp`) with PLUMED 2.7.3. Inputs use BLYP-D3(BJ), DZVP-MOLOPT-SR-GTH/GTH-BLYP, 300 Ry density and 40 Ry relative cutoffs, and EPS_SCF = 10⁻⁶. NVT dynamics use a 1 fs step and CSVR at 298.15 K with a 100 fs time constant. Production calculations used 64 MPI ranks and one OpenMP thread per rank. Coordinates are written every 10 fs and energies every step.

Ordinary fixed-height metadynamics biases the S293–C289 distance. Hills have width 0.1 Å and are deposited every 30 fs, with heights of 0.03137545 kcal/mol for the interface and 0.0627509 kcal/mol for bulk. The PLUMED inputs specify the 4 Å upper wall and KAPPA = 200. The exact as-run parameters are retained in the input and PLUMED files.

The bias is removed in the separate unbiased runs. Their starting structures can still carry metadynamics preparation history. The use of ordinary metadynamics, wall strength, thermostat time constant, exact basis labels and charged-slab treatment are among the points to compare with the authors' original implementation before interpreting a free-energy difference as a reproduction of the published barrier.

## Original-coordinate and water-orientation controls

The original controls use the actual 1150 and 2200 fs frames of the 5 ps R pre-equilibration (job 8568513), with seeds 2026111603 and 2026111604. The donor/reverse pairs rotate H202/H203 around water O204; the original controls retain all 296 source coordinates. Within each source pair the same velocity protocol and seed are used. Each production run starts afresh from the specified coordinates after a separate 20-step smoke run, so smoke frames are not appended to production.

## Trajectory observables

`analyze_event_unbiased.py` calculates all-frame minimum-image distances in the supplied box. Its original filename is retained even when the wrapper applies its distance calculations to a biased trajectory; the wrapper records the actual bias history.

`analyze_solvation.py` reports carbon geometry, local water coordination, hydrogen bonds, orientation and a signed distance from an operational Gibbs dividing surface. The surface is estimated from the smoothed water-oxygen density profile. Positive signed distances are on the vapour side, negative distances on the liquid side. This is a density-defined surface rather than an instantaneous molecular surface.

`analyze_network.py` uses directed water-hydrogen/acceptor contacts. The strict, standard and loose distance/angle thresholds are 2.4 Å/160°, 2.5 Å/150° and 2.6 Å/140°, respectively. Parent water identities and threshold sensitivity are retained. Geometric contact spans are reported as observations, not fitted hydrogen-bond lifetimes. Numerical geometric-screen definitions are in `campaign/li2025_structural_gate_20261116/thresholds.json`.

## Author workbook

Download `source-data.xlsx` from the Zenodo record if you need to repeat coordinate extraction. Its recorded size is 53,375,254 bytes, MD5 is `53a0f1915b8f15fe1cd5d56742d16820`, and SHA-256 is `f4e2bedc004f0abaf180f55e95bd2aa9f7ddc2aa86ef09a683fd2d5afe0a9992`. The author requests citation of the paper and relevant literature when using the reference geometries. The workbook and article PDFs are not part of this release.

The source code snapshot is commit `c3693948a16bb02a85b48b6bfc8b3b0ad228896b` of the working project. `metadata/SOURCE_FILES.json` identifies the exported files byte for byte. New portable inputs have separate hashes and a recorded path-only diff.
