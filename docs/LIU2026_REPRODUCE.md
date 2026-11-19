# Recalculate the Liu trajectory descriptors

The Liu calculations use a 394-atom water slab in a 17 × 17 × 47 Å periodic cell. They are a separate system from the 296-atom Li calculations. The selected trajectories and their preparation histories are listed in `data/LIU2026_DATASETS.json`.

After obtaining the Liu data package, list the available trajectory identifiers:

```sh
python3 -c 'import json; print("\n".join(json.load(open("data/LIU2026_DATASETS.json"))["cases"]))'
```

Choose an identifier from that list and use a new output directory:

```sh
python3 tools/reproduce_liu_case.py --case frame_ic_1 --out results/liu_frame_ic_1
```

The command needs Python 3.9 or newer and the standard library. The repository's existing Python environment also works. It reads the released data, checks their hashes, and leaves the source files unchanged. It refuses an existing output directory.

## Outputs

- `FRAME_METRICS.tsv`: one row per unique saved MD step, with S389–C385, C385–O386, O390–H393 and O386–H393 distances; the nearest hydrogen to O386 and its atom index; temperature where an energy record has the same step and time; and signed distances to the water surface for C385, S389, O386 and K394.
- `SUMMARY.json`: distances, temperature, signed surface positions and vapour-side fractions, with separate summaries for declared bias stages and explicit trajectory coverage.
- `CASE_MANIFEST.json`: the case entry used for this calculation, including preparation history and stage boundaries.
- `PROVENANCE.json` and `OUTPUT_SHA256.json`: the command, source code identity, source hashes before and after reading, and output hashes.

The tool checks all 394 atom identities and finite coordinates in every original frame. The archived parser keeps the last occurrence of a repeated MD step and sorts the resulting steps. Repeated endpoint counts, coordinate differences and time conflicts remain in the summary. Incomplete frames cause an error. Missing or inconsistent temperature records are left empty, without interpolation.

## Definitions

Distances use the original minimum-image function, periodic along all three cell axes. Atom numbers are 1-based and follow the original trajectory order. A short distance is reported as a distance; this command does not assign a bond, a reaction event or a successful reproduction.

The surface calculation uses the original `analyze_iface.frame_metrics` function: water-oxygen density in 0.5 Å bins, a three-bin moving average and the half-maximum density boundary. Signed distance is positive on the vapour side and negative inside the water slab. The reported vapour-side fraction uses only frames in which the surface calculation resolves the slab; unresolved frames and their reasons are retained. This is an operational estimate of the surface position.

The nearest-H distance considers all hydrogen atoms and retains the corresponding identity at every frame. It does not track a proton-transfer pathway. The water-mediated distances `d2`, `d3` and `d5` are computed only when a case supplies the original mediator-water atom map. Otherwise those fields remain empty. Water coordination numbers are not added because the four archived core files do not supply a shared coordination definition.

## Stage boundaries and coverage

The four proton-location A+B trajectories combine an initially restrained segment with an unbiased segment. Full A+B entries therefore carry `biased: "mixed"`, and each row retains its own stage. In these entries, steps 0–1000 belong to preparation stage A and steps 1001–6000 to stage B.

The corresponding B-only files contain 500 saved frames from steps 1010 to 6000. Their observed span is 4990 fs, although the planned integrated stage lasted 5000 fs. The tool preserves the original step and time origin. A B-only file and its parent A+B file contain overlapping observations and are not independent trajectories.

Other selected runs include initially prepared frame controls, a previously selected seed, and weakly or strongly guided trajectories. Their bias and preparation flags are copied from the dataset catalog. Temperature and trajectory endpoints are reported separately; the tool does not infer missing coverage from another file.

The four original source files are unchanged. Their commit and hashes are recorded in `metadata/LIU2026_SOURCE_FILES.json`. The public wrapper calls their parsing and descriptive geometry functions, without running their historical decision routines.

To test the wrapper in an isolated copy, using any downloaded trajectory with at least two frames:

```sh
python3 tools/test_liu_tools.py --case frame_ic_1 -v
```

The tests independently recalculate S–C and C–O distances for the full selected trajectory and exercise duplicate endpoints, incomplete frames, atom identities, source hashes, stage boundaries and missing temperature or surface values.

The release test record, `examples/LIU2026_REPLAY_TEST_REPORT.json`, covers the full saved trajectories for `proton_r1`, `proton_r1_full_ab`, `frame_ic_1`, `original_seed_unbiased_5ps` and `weak_guided_4ps`. It records the repeated observations shared by the two proton_r1 views and compares available historical distance values with the replay.
