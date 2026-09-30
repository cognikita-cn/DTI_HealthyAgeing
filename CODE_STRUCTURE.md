# DTI Healthy Ageing — post-review code structure

## Active post-review graph pipeline

1. `06_Combined_Graph_theory_analysis_setup.m`
   - loads both datasets;
   - uses the final cross-dataset tractography matrices;
   - defines healthy young/older selections;
   - builds the healthy-derived combined edge masks;
   - builds the final permutation caches.

2. `04_Mass_matched_lesion_and_permutations.m`
   - current mass-matched lesion null;
   - current distance-matched lesion null;
   - combination/conversion utilities;
   - current rich-club null launcher and combiner.

3. `06_Combined_Graph_theory_analysis.m`
   - healthy-only primary graph inference;
   - mass-matched null;
   - distance-matched sensitivity null;
   - revised rich-club definition;
   - supplementary tables.

## Active functions

- `functions/load_conn_mat_dem.m`
  - current connectivity loader.
  - Same calculations as the original loader, with explicit file/count checks.

- `functions/make_permutation_cache.m`
  - current cache builder.
  - Keeps all valid-DTI participants and stores `healthy_selection`.

- `functions/get_aal_centroid_distances.m`
  - current distance helper for the distance-matched null.
  - Requires SPM12 on the MATLAB path.

- `functions/perm_les_func_streams.m`
  - current mass-matched lesion null.

- `functions/perm_les_func_streams_distance.m`
  - current five-bin centroid-distance-matched sensitivity null.

- `functions/rich_club_perm_current.m`
  - replacement for the old rich-club permutation function.
  - Uses current caches and current saved masks with reproducible seeds.
  - Retains the historical BCT null model `null_model_und_sign(M3,5,0.1)`.

## Rich-club status

The old rich-club permutations are not valid final post-review inputs because they were
generated from the old threshold/matrix workflow. `run_rich_club_analysis` is therefore
disabled by default in the main analysis. Recompute and combine the current rich-club
nulls from `04_Mass_matched_lesion_and_permutations.m`, then enable the analysis.

## TBSS / cross-dataset hotspot logic

The final manuscript uses **within-dataset TBSS discovery on healthy participants**,
followed by **cross-dataset application of the resulting age-sensitive cluster**.

- TempSeqAges graph analysis uses the AudMemDem-derived cluster.
- AudMemDem graph analysis uses the TempSeqAges-derived cluster.
- The pooled `tbss_combined` branch is not part of the final analysis.

Active TBSS material is under `tbss/`. The supplied `ttest_perm_Final.m` is the
AudMemDem dataset-specific permutation workflow. The matching TempSeqAges-specific
cluster-generation/permutation code was not included in the uploaded function set and
should be added before final archival if full raw-to-result reproducibility is desired.

The historical pooled TBSS scripts are retained in `legacy/` for provenance only.

## Utilities

- `utilities/tbss_size_checker.m`
  - diagnostic only; it checks the size of the combined TBSS 4D file.

## Legacy

The following are retained only for provenance and should not be called by the
post-review manuscript pipeline:

- `legacy/05_Graph_theory_analysis_legacy.m`
- `legacy/04_permutation_workflow_pre_review_legacy.m`
- `legacy/perm_les_func_streams_old.m`
- `legacy/perm_les_func_prop_legacy.m`
- `legacy/rich_club_perm_new_legacy.m`

`perm_les_func_prop.m` is especially unsafe as an active function because its internal
function name is `perm_les_func_streams` and its implementation is another mass-style
streamline permutation rather than a clearly separate proportion-null implementation.

## Remaining provenance checks before repository freeze

- Add or identify the TempSeqAges-specific TBSS permutation code that generated the
  cluster later applied to AudMemDem.
- Confirm that `03_Probabilistic_tractography.m` contains, or is paired with, the exact
  cross-dataset waypoint/tractography commands that generated both final cross-check
  matrix sets.
- Recompute the rich-club nulls with the current masks before reporting final RC results.
- Keep the historical preprocessing script unchanged; participant 0001 was a pilot and
  some preprocessing loops intentionally begin at participant 2.
