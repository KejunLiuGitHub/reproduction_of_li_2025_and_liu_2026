# Li 2025 reproduction: current results and continuation

We have obtained short-C–S states in interfacial metadynamics and followed two of these states for 5 ps after removing the bias. We have **not yet reproduced unbiased interfacial C–S formation starting from separated reactants**. The next step is to finish the two running reactant controls and compare them with the completed water-orientation experiments.

Completed-run results and the matched-window comparison below were compiled on **19 November 2026 at 21:06:23 UTC+08:00**. The release includes later running-control snapshots, described in the coverage table. It accompanies a selected reproducibility archive containing complete chosen trajectories and their associated inputs and analysis, including runs that did not produce the sought geometry. The archive selection manifest records its contents and omissions; `metadata/DATA_FILES.json` records the released files and their hashes.

## System and calculation

The target is the SO₃²⁻ + HCHO branch of Li *et al.*, *Nature Communications* **16**, 5187 (2025), [paper](https://doi.org/10.1038/s41467-025-59712-3). Starting coordinates come from the authors’ [public source data](https://doi.org/10.5281/zenodo.15253851). The interfacial cell contains 96 waters, HCHO and SO₃²⁻: 296 atoms, charge −2, singlet, a 15 × 15 × 30 Å box, and no explicit counterions.

The calculations use CP2K with BLYP-D3(BJ), DZVP-MOLOPT-SR-GTH/GTH, 300/40 Ry cutoffs and EPS_SCF 10⁻⁶. Dynamics use a 1 fs timestep and CSVR at 298.15 K with a 100 fs time constant; coordinates are saved every 10 fs. Metadynamics uses S–C distance as the collective variable, fixed-height hills and a 4 Å upper wall. The subsequent unbiased trajectories remove the metadynamics and distance restraints.

The exact basis labels, wall strength, thermostat time constant and choice of ordinary rather than well-tempered metadynamics are implementation choices whose equivalence to the authors’ original inputs remains to be resolved. These runs have not established the paper’s free-energy barrier or interfacial acceleration. The present comparison concerns geometries, solvent environments and their time evolution.

## Metadynamics reached short-distance states

The original interfacial chain reached **50 ps**. Its continuation jobs were 8568711, 8571154 and 8572387. A separate chain, started after 5 ps of reactant pre-equilibration (8568513), reached **20 ps** through jobs 8569493 and 8571156. The two bias histories are retained separately.

In the original chain, the 39–40 ps window had mean S–C/C–O distances of **1.815/1.386 Å**, with S–C below 2 Å in 99 of 100 saved frames. Later sampling returned to longer S–C distances: the 40–50 ps segment had 376 of 1000 frames below 2 Å, mean S–C/C–O of 2.3081/1.3243 Å, and final S–C of 3.1214 Å. The final 49–50 ps window contained no saved S–C distance below 2 Å.

The separate R-started chain also visited short-distance states; its final S–C at 20 ps was 2.3762 Å. These observations motivated releasing selected configurations into unbiased dynamics. The metadynamics trajectories are used to locate and compare states; their elapsed times are not physical reaction times. Both chains stopped at their planned endpoints.

## Two selected states retained short C–S distances for 5 ps

Two naturally sampled metadynamics configurations were used without coordinate editing, with newly initialized velocities and no bias during the ensuing 5 ps:

| Trajectory | Starting configuration | Mean S–C / C–O, Å | Saved frames with S–C <2 Å | Final S–C, Å |
|---|---|---|---:|---:|
| `natural_original40`, 8572454 | Original-chain frame at 40 ps | 1.9128 / 1.3661 | 439/501 | 1.8239 |
| `natural_r19830`, 8572455 | R-chain frame at 19.83 ps | 1.8904 / 1.3826 | 478/501 | 1.8767 |

This is the strongest retention result obtained so far. Both trajectories stayed predominantly in the short-distance region over the full 5 ps, with excursions across the 2 Å threshold. Under the combined distance and local-carbon geometry definition, their final continuous segments span **750 and 1190 fs**. Tightening the S–C threshold to 1.95 Å leaves a terminal segment of at least 100 fs only for the second trajectory. The first result is therefore more sensitive to the geometric definition.

The starting S–C distances were already 1.8523 and 1.9683 Å. These trajectories show what happens to selected short-distance states after removal of the bias; they do not measure formation from separated reactants. Each starting condition has one trajectory.

The two states also occupy different depths in the water film. Mean signed distances of C/S/Oc from the density-defined Gibbs dividing surface are **−2.960/−3.518/−3.514 Å** and **−4.440/−4.777/−4.790 Å**, respectively. All three atoms remain on the liquid side in the saved frames. Carbon P10–P90 ranges are −3.303 to −2.559 Å and −5.049 to −4.022 Å. These positions should be retained when comparing solvent environments: the second configuration lies deeper, and neither is centered exactly on the dividing surface.

## Earlier controls explain why the starting environment matters

Public product-labeled snapshots retained short distances more often than the first short-distance event taken from our own metadynamics. The following table reports distance observations over each complete trajectory.

| Starting condition | Job(s) | Duration / saved frames | Frames with S–C <2 Å | Mean S–C, Å |
|---|---|---|---|---|
| Public interface P | 8568510 | 3 ps / 301 | 227 | 1.9596 |
| Public bulk P | 8568511 | 3 ps / 301 | 291 | 1.8882 |
| Public interface P, new velocities | 8569494 | 5 ps / 501 | 437 | 1.9341 |
| Our interface event at 10.41 ps, three velocity seeds | 8566899, 8568085, 8568086 | 3 ps / 301 each | 2 / 5 / 2 | 3.0470 / 2.5834 / 2.9459 |
| Same event after 2 ps distance-restrained relaxation, inherited velocities | 8568949 | 5 ps / 500 | 7 | 3.0305 |
| Same relaxed endpoint, new velocities | 8569495 | 5 ps / 501 | 26 | 2.6300 |

The 2 ps restrained preparation had mean S–C/C–O of 1.9914/1.3345 Å and final S–C of 1.9963 Å. After release, however, both velocity choices moved predominantly to longer S–C distances; neither had a saved S–C distance below 2 Å in its final 1 ps. Simply bringing S and C close together did not reproduce the retention seen in the later metadynamics-selected states.

The public-P 5 ps trajectory contains five frames where the fixed parent-water connectivity becomes ambiguous. Its distance record is retained, with the affected atom identities and frames documented. The legacy inherited-velocity release lacks its stage-t0 XYZ frame and an explicit return-code field in the archived run record; its 5 ps integration and 500 saved frames are distinguished in the archive. These details matter when replaying or comparing the trajectories.

The differences between public-P, early-event releases and later retained states point to the surrounding configuration as a useful variable. They do not identify water coordination alone as the cause: solute orientation, depth, carbon geometry and solvent organization change together.

## Long-distance reactant tests and the two running controls

To examine water orientation without starting from a short-C–S state, two actual frames were selected from the 5 ps reactant trajectory: 1150 and 2200 fs. Their initial S–C distances are **3.8711 and 4.2418 Å**. For each frame, H202/H203 of water O204 were rotated in opposite directions while the other 294 atoms remained unchanged. All four runs then evolved for 5 ps without bias or restraints.

| Job / condition | Mean / minimum / final S–C, Å |
|---|---|
| 8572523, R1150 donor | 3.4148 / 2.6307 / 3.4646 |
| 8572524, R1150 reverse | 4.4117 / 3.1595 / 6.1397 |
| 8572525, R2200 donor | 3.1311 / 2.3902 / 3.3554 |
| 8572526, R2200 reverse | 3.8503 / 2.6855 / 3.8811 |

None produced a predefined C–S formation geometry in its 501 saved frames, under any of the three geometric definitions. Thus these four particular water-orientation preparations did not produce the sought event within 5 ps. The original reactant pre-equilibration likewise had no such event.

Two unrotated controls now complete the original/donor/reverse comparison. They use the same source frames and pair-specific velocity seeds, with all source coordinates unchanged:

| Job / condition | Coverage in the 21:58 snapshot (UTC+08:00) | Mean / minimum / final S–C, Å |
|---|---|---|
| 8572577, R1150 original | 294 frames, 0–2.930 ps | 3.6160 / 2.7181 / 3.5676 |
| 8572578, R2200 original | 311 frames, 0–3.100 ps | 3.4989 / 2.8017 / 3.6076 |

Both are still running toward 5 ps. Neither prefix contains a saved S–C distance below 2 Å. The earlier matched-window comparison uses the first 2570/2720 fs of the released original-control trajectories:

| Source / common window | Original / donor / reverse mean S–C, Å | H203→Oc contact frames |
|---|---|---|
| R1150, 0–2570 fs | 3.6770 / 3.6748 / 3.9145 | 193 / 129 / 5 of 258 |
| R2200, 0–2720 fs | 3.5113 / 3.2741 / 4.0683 | 170 / 217 / 213 of 273 |

For R1150, original and donor have almost the same mean S–C despite different target-contact counts. For R2200, donor and reverse have similar target-contact counts but different mean S–C. Both rotations also raise the initial potential energy relative to the corresponding original coordinates. The complete comparison therefore needs the actual water contacts, solute depth and orientation together, rather than treating the initial donor/reverse labels as permanent network states.

These are selected frames from one R history, involving the same local water and one trajectory per condition. Adjacent saved frames are correlated. The current prefixes describe their evolution so far; the final comparison awaits the remaining 5 ps coverage.

## Continuing the study

1. **Complete the two existing original controls.** Compare all three arms over the full 5 ps and the fixed 0–0.1, 0.1–1 and 1–5 ps windows, retaining every observed outcome. The first two windows are already covered; the last is still partial.
2. **Compare retained and non-retained configurations.** Use S–C/C–O, carbon geometry, water identities, orientation and signed-GDS depth together. A focused electronic-structure characterization of representative retained states would clarify the nature of the short-C–S interaction. Existing Mulliken/Hirshfeld populations are population-scheme descriptors, not bond orders.
3. **Resolve method differences before extending the free-energy claim.** The charged-slab treatment, finite-size and force convergence, publication-unspecified input choices and bias convergence remain the relevant checks for reproducing the published barrier. More elapsed metadynamics time alone would not settle these questions.
4. **Keep future comparisons interpretable.** Any further trajectories should address a specific difference revealed by the completed controls, with independently prepared starting conditions where replication is needed. The present results give no reason to strengthen the distance restraint or continue a broad rotation search.

The archive preserves full selected trajectories, including the returns to long S–C distances and the four unsuccessful formation trials. It also keeps the incomplete status of the two running controls. This makes the distinction between the results already obtained—biased short-distance states and subsequent unbiased retention—and the remaining objective—formation from separated reactants under unbiased interfacial dynamics—explicit.
