# Local hydrogen-bond networks and natural-candidate selection

Method recorded on 16 November 2026 for the Li Figure 6b system: SO₃²⁻, HCHO and
96 H₂O molecules, with the fixed 296-atom order. The analysis reuses the atom mapping,
periodic minimum-image convention and full-frame geometrical descriptors in
`analyze_solvation.py`.

## Directed edges and local shells

Each edge is the one-based identity tuple `(donor water O, H, acceptor O)`. Acceptors
include all water oxygens and the solute oxygens Oc and Os1–Os3; the donor's own oxygen
is excluded. The standard definition is H···O ≤ 2.5 Å and O–H···O ≥ 150°.
The strict and loose definitions are 2.4 Å/160° and 2.6 Å/140°, respectively.
These are project definitions for geometrical hydrogen-bond observations.

A parent-water O–H distance above 1.3 Å is flagged using the existing assignment.
All edges and frames are retained, alongside a sensitivity summary excluding flagged
frames. Hydrogen ownership is never silently reassigned.

A separate graph is constructed for each threshold. Waters donating directly to Oc
or Os form the first graph shell. Waters connected to those donors by a water–water
hydrogen bond in either direction form the second shell. The local graph is the
subgraph induced by these two shells, together with water-to-solute edges. These graph
shells differ from the geometrical coordination shell defined by a 3.5 Å cutoff.

Donor/acceptor degree families are counted in the full graph, including bifurcated
hydrogens, and summarized for local waters. Both directed and undirected shortest
water paths are reported between the Oc-donating and Os-donating waters. A path length
of zero denotes a shared donor water. Connectivity is a geometrical description;
a path alone does not establish proton transfer.

## Sampled edge observations

Continuous segments follow each edge's identity across frames. Their endpoint span is
`(last step − first step) × 1 fs`; a single-frame observation therefore has zero span.
Intermittent segments allow at most one intervening sampled frame without that edge
(10 fs). This explicit gap rule is used without fitting a lifetime distribution.
Left and right censoring are recorded for every segment.

An edge is observed only while it satisfies the hydrogen-bond definition and belongs
to the local graph. Leaving the local region ends its observation even if the same
hydrogen bond persists elsewhere. Spans from biased trajectories are reported as
sampled geometrical observations rather than physical lifetimes or reaction rates.
Inputs must have complete, uniformly spaced 10-step frames, recorded source hashes,
and agreement with the prior solute hydrogen-bond counts at every frame and threshold.

## Analysis set and threshold-entry frames

The analysis set contains the public-P fresh-velocity 5 ps run, the original-chain
30–40 ps segment, the independent-R 10–20 ps segment, and the two 5 ps releases after
restrained relaxation. Summaries use all frames and, separately, short-SC entry frames.
An entry is the first frame with SC < 2 Å immediately following a frame with SC ≥ 2 Å.
A trajectory segment beginning at short distance is marked as left-censored.
Candidate selection does not use future retention duration.

## Network groups and matching

The two predefined network groups have an undirected shortest path between Oc- and
Os-donating waters of either at most two edges, including a shared donor, or more than
two edges/no local connection. The labels describe connectivity rather than favorable
or unfavorable reaction conditions. Only entry frames without parent-OH flags enter
the candidate pool.

Matched candidates must have identical standard Oc hydrogen-bond counts and identical
integer Oc water-coordination counts. Additional absolute-difference limits are:

| Descriptor | Maximum difference |
|---|---:|
| SC distance | 0.05 Å |
| CO distance | 0.03 Å |
| CO and CS normal cosines, each | 0.10 |
| C and Oc signed-GDS distances, each | 0.30 Å |
| S signed-GDS distance | 0.50 Å |
| OCS angle | 5° |
| Carbon out-of-plane height | 0.03 Å |
| Temperature | 10 K |

Cross-group matches are sought separately within the original and independent-R
histories. Frames from the same history must be at least 1 ps apart; this removes
nearby frames but does not establish statistical independence. A four-run stratified
exploratory set is proposed only if both preparation histories supply one matched
pair. Network decorrelation and preparation-history effects remain relevant to its
interpretation.

## Fallback retention pilots

If four strictly matched candidates are unavailable, the thresholds stay fixed.
Instead, one latest natural frame with SC < 2 Å and no parent-OH flag is selected from
each complete original-chain 30–40 ps and independent-R 10–20 ps trajectory. All
matching descriptors and their differences are recorded. These two candidates are
explicitly unmatched and have different biased preparation histories.

Each pilot starts a fresh 5 ps unbiased trajectory from one selected configuration.
The comparison concerns geometrical retention after removing the bias. With one
initial condition per history, it is not a matched network-causality experiment or
a replicated equilibrium sample. The predefined seeds are 2026111601 for the original
chain and 2026111602 for the R chain. Provenance retains `prepared_under_bias=true`;
the event10410 configuration is not reused.

Source XYZ files are identified by stable paths and hashes. Candidate records bind
the selected frame step, source job and trajectory SHA to the source run's numerical
checks. This selection is independent of the ongoing original-chain continuation:
that chain ends at a cumulative 50 ps, while the independent R-metaD chain ends at 20 ps.
