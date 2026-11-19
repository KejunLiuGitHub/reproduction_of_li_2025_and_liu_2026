# Liu 2026 reproduction: results and remaining questions

Our calculations have not yet reproduced unbiased interfacial C–S formation from the intended reactant complex or the reported interfacial free-energy barrier. This does not demonstrate that the published results are incorrect. The paper makes the original data [available on request](https://acp.copernicus.org/articles/26/9357/2026/); the public article and supplement do not include complete input/restart files or numerical starting coordinates for the reported trajectories. The [supplement](https://acp.copernicus.org/articles/26/9357/2026/acp-26-9357-2026-supplement.pdf) describes quadratic metadynamics walls without specifying their positions or force constants. These original records are needed to resolve the run-specific setup choices.

The calculations did produce short C–S distances under applied restraints and retention of short distances after controlled proton placement. These are useful observations about the prepared configurations, and do not complete the reproduction target.

The final sampling record is the stop on 10 November 2026, reviewed on 11 November. The selected archive combines complete saved trajectories for the chosen dynamical tests with the archived collective-variable and bias histories of three stopped metadynamics chains. It contains retained short-distance states, controls without retention, and unfinished experiments. Full metadynamics coordinate trajectories and other unselected collections are described in the additional-results catalogue.

## Target and molecular model

The target is Liu et al., *Atmospheric Chemistry and Physics* **26**, 9357–9371 (2026), [DOI 10.5194/acp-26-9357-2026](https://doi.org/10.5194/acp-26-9357-2026), “Interface-dominated hydroxymethanesulfonate and its isomer formation provides key mechanisms for reconciling the atmospheric sulfur budget gap in polluted and cold environments.” The paper reports interfacial HMS formation within approximately 3 ps in unbiased dynamics starting from an optimized reactant complex, and an interfacial HMS forward free-energy barrier of 0.6 kcal/mol. The corresponding bulk HMS, bulk HMSi and interface HMSi barriers are 7.7, 7.6 and 6.1 kcal/mol.

Our Liu model contains **394 atoms**: 128 waters, HCHO, HOSO₂⁻ and K⁺, with total charge zero and singlet spin. The interfacial cell is 17 × 17 × 47 Å. C385 and S389 define the C–S distance; O386 is the formaldehyde oxygen; O390/H393 identify the initial acidic O–H group; K is atom 394. Atom 393 is a hydrogen identifier, not the system size.

This is chemically distinct from the Li 2025 branch in the same repository: Li uses SO₃²⁻ + HCHO, 96 waters, 296 atoms, charge −2 and a 15 × 15 × 30 Å interfacial cell. Their trajectories, atom mappings, bias histories and free-energy interpretations must remain separate.

The Liu dynamics use CP2K/PBE-D3, DZVP-MOLOPT-SR-GTH/GTH, 300/40 Ry cutoffs, NVT CSVR at 298.15 K, a 1 fs timestep and a 100 fs thermostat time constant. The metadynamics coordinate is the five-distance combination `d1+d2−d3−d4−d5`, with `d4 = r(S389,C385)`. Well-tempered settings are hill width 0.4 Å, height 0.16 kcal/mol, deposition interval 30 fs and bias factor 10. Those hill parameters and the coordinate expression were checked against the paper's supplementary information. The implemented walls at −14 and +10 Å with KAPPA 200 are project choices; their equivalence to the authors' walls has not been established.

## Three metadynamics chains stopped without reaching short C–S distances

The three chains were stopped on 10 November at 00:40 UTC+08:00 after persistently sampling non-target configurations. The values below come from the final archived COLVAR records, including their original segment-boundary records.

| Chain | Job running at stop | Last recorded time | COLVAR rows | Minimum / final S–C, Å | Rows with S–C <2 Å |
|---|---:|---:|---:|---|---:|
| Original interface chain | 8554932 | 60.650 ps | 6069 | 2.6932 / 4.9823 | 0 |
| RC-history frame 340 | 8562431 | 19.290 ps | 1929 | 3.1754 / 7.5164 | 0 |
| RC-history frame 890 | 8562419 | 13.610 ps | 1363 | 2.5159 / 4.3177 | 0 |

These are three histories from the same RC preparation lineage, not three independently constructed equilibrium reactant ensembles. The original chain's median CV was −11.422 Å, with 23.68% of recorded values below −13 Å. The carbonyl-to-designated-water proton distance `d3` had a median of 7.217 Å; 94.48% of its values exceeded 4 Å. The two alternate starts also spent most of their records at `d3 >4 Å` (79.37% and 75.13%).

The component distances explain why elapsed metadynamics time alone was a poor measure of progress: the algebraic coordinate changed substantially while the designated solvent contact became long and S–C remained outside the short-distance region. This identifies a sampling problem worth investigating. It does not establish that a narrower wall or a different coordinate would recover the paper's mechanism.

The planned sampling and converged free-energy comparison were not completed. No barrier from these chains is presented as a reproduction of 0.6 kcal/mol. The public selection retains the complete locally archived COLVAR/HILLS histories and the actual input records, **not the full coordinate trajectories of these three chains**. Segment-boundary overlaps remain in the original files and must be handled explicitly when counting unique times.

## Unbiased trajectories from reactant-like starts

The original prepared RC began at S–C = 3.563 Å. Its uninterrupted physical history was extended from 3 to 5 ps through jobs 8554054 and 8554134. The full saved trajectory contains 501 frames over 0–5000 fs. It contains no S–C distance below 2 Å and ends at 5.427 Å. This is a completed unsuccessful formation test for that prepared RC.

A separate experiment selected five actual configurations at 1, 2, 3, 4 and 5 ps from one orientation/preparation trajectory. Jobs 8554943–8554947 each completed a 3 ps first segment, with 301 saved frames. Their wrappers then failed because they looked for a restart filename different from the one written by CP2K. The restart was present at step 3000. Recovery jobs 8562388–8562392 were subsequently cancelled after only very short second segments or before starting.

The archive retains all five complete 3 ps first segments. They are **unfinished 5 ps experiments**, not five completed negative results. Their source configurations come from different times in one history, whose preparation constraints must be retained; independence or equilibrium sampling cannot be inferred simply from the spacing between selected frames.

Other numerical interruptions also occurred. A 10 ps pre-equilibration attempt (8554980) hit its wall-clock limit: its energy record reached step 3225 and its saved restart step 3200. Metadynamics continuation segments also encountered wall-clock limits. These interruptions are distinct from the physical absence of the sought geometry in a completed trajectory.

## Guided trajectories and controlled proton-location tests

Two guided examples document the response to external forces. In the 4 ps weak-guidance trajectory (8554057, extended by 8554111), the final continuous S–C <2 Å segment begins at 2990 fs and lasts 1010 fs. In the 3 ps strong S–C restraint run (8554102), the corresponding segment begins at 300 fs and lasts 2700 fs; the minimum S–C is 1.685 Å. Both are biased simulations. Their event times are not unbiased reaction times or barrier measurements.

The later four-arm experiment began from one D15-derived configuration. It changed the placement of H393 between control and proton-prepositioned arms, used paired velocity seeds, then applied a 1 ps S–C restraint at 2.0 Å with KAPPA 10 before a 5 ps integration with all PLUMED forces removed. Proton relocation changes the initial chemical connectivity as well as the local solvent environment; it is more than a rotation-only perturbation.

| Arm | Job | Mean S–C in the unbiased segment, Å | Saved frames with S–C <2 Å | Final continuous short-distance span |
|---|---:|---:|---:|---:|
| Proton-prepositioned, seed 1 | 8555023 | 1.8436 | 499/500 | 3.960 ps |
| Proton-prepositioned, seed 2 | 8555025 | 1.8579 | 500/500 | 4.990 ps |
| Control, seed 1 | 8555022 | 3.7964 | 0/500 | None |
| Control, seed 2 | 8555024 | 3.4019 | 0/500 | None |

The proton-prepositioned arms were already at short C–S distance when the unbiased segment began. They provide evidence of **retention after controlled preparation**, with two velocity seeds. They do not demonstrate spontaneous formation from the unreacted RC. The two controls make the difference between these particular prepared states clear, but two seeds of one coordinate construction do not establish an equilibrium yield or isolate a general protonation effect.

Each complete preparation-plus-release file has 601 saved frames from step 0 to 6000. The separate unbiased-segment file contains 500 frames from step 1010 to 6000: a 5 ps integration with a 4.990 ps interval between its first and last saved frames. Both versions are retained to preserve the existing analyses and their hashes; the shared frames must not be counted twice. Relative-time analyses accompany the original records because some earlier labels treated absolute step numbers or the missing phase-boundary frame as a coverage failure.

## Proton geometry, depth and orientation

The relevant observables include C–S/C–O distances, O390–H393, the closest hydrogen to O386, solvent contacts and signed distance from a density-defined interface. A short C–S distance alone does not establish the complete proton-relay pathway. After H393 is relocated or moves away from O390, the vector O390→H393 is also no longer the direction of a bonded O–H group, so its orientation cannot be compared across arms as though it always represented the same hydroxyl.

The public reanalysis uses periodic distances and the retained core interface calculation. Older diagnostic scripts used a different surface binning convention, contained a lower-interface sign inconsistency, or computed proton distances without periodic boundaries. Their unique quantitative orientation and relay assignments are therefore not used here to establish the mechanism or reproduce the paper's product-depth distribution. The original outputs remain available as historical diagnostics, with these limitations recorded.

## What would make a continuation informative

The original reactant-complex coordinates and associated input/restart files are not included in the public article or supplement; the stated route to these data is a request to the authors. Those records would resolve the exact initial states of the unbiased and metadynamics trajectories. The preparation documented here instead used selected explicit-solvent frames, restrained solvent contacts, small-cluster optimization and embedding adjustments.

A useful continuation would first obtain or reconstruct those inputs and resolve the walls, dispersion damping, thermostat and free-energy reconstruction details. Component-coordinate sampling should then be checked before committing to long metadynamics. A modified coordinate or confinement scheme would be a separate methodological experiment with a fresh bias history. Further unbiased tests should use clearly documented RCs and complete their stated observation windows.

The target results have not been reproduced in these calculations. The retained trajectories document the tested initial conditions and the gaps a continuation would need to resolve. Separate bulk, static and machine-learned-potential collections are listed in the on-request catalogue.

## Data and source files

`data/LIU2026_DATASETS.json` distinguishes 16 coordinate views from 12 experiments and three additional metadynamics CV histories. Four phase-B-only views overlap the corresponding full A+B files. `metadata/LIU2026_DATA_FILES.json` records every released file and its source hash. Four PLUMED comment headers were translated into English; their executable directives are unchanged and both hashes are recorded. Small input, energy, seed and checkpoint dependencies were recovered by exact hashes from the original run records.
