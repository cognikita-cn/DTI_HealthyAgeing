# DTI_Healthy Ageing

This repository contains full code used for the "Efficient ageing: Simulated lesion of the structural connectome reveals optimised decline in the healthy ageing brain" work. 

NOTE: code is shared for transparency and is not necessarily presented in a reusable format. In case you have any questions / comments / suggestion on the current code, please feel free to reach out on nikud@clin.au.dk


MATLAB scripts for diffusion-MRI preprocessing, age-related white-matter cluster analysis, probabilistic tractography, and structural-network virtual lesions.

## Files

| File | Role |
| --- | --- |
| `Pipeline.m` | Section-based workflow covering preprocessing, voxelwise statistics, tractography, network analysis and rich-club dependence |
| `ClustperimInfAcrDens.m` | TempSeq analysis script comparing network measures and permutation outcomes across nine density percentiles |
| `Functions/ttest_perm_func.m` | Voxelwise age permutations and descending cluster-size distributions |
| `Functions/load_conn_mat.m` | Loading, symmetrisation and organisation of raw, cluster-constrained and virtual-lesion matrices |
| `Functions/perm_les_func_streams.m` | Multinomial streamline-count lesion permutations |
| `Functions/perm_les_func_prop.m` | Permutation of cluster-dependent proportions across edges |
| `Functions/rich_club_perm_new.m` | Random-network rich-club coefficient curves |

`ClustperimInfAcrDens.m` is a script; the five files under `Functions` are MATLAB functions. Function help comments describe the fields read from their input structures and the files they write.

## Workflow

| Stage | Processing | Principal outputs |
| --- | --- | --- |
| Diffusion preprocessing | AP/PA DICOM conversion, TOPUP, brain extraction, EDDY and tensor fitting | Corrected diffusion data, masks and FA images |
| TBSS and voxelwise statistics | Registration, skeleton projection, age-group comparisons and cluster-size permutations | Skeletonised FA and significant white-matter masks |
| Probabilistic tractography | BEDPOSTX, spatial transforms, AAL90 ROI preparation and tractography with/without a cluster waypoint | Raw and cluster-constrained connectivity |
| Virtual lesions | Cluster-weight subtraction, edge selection, ROI-size normalisation and null-network comparisons | Efficiency, node-strength and clustering summaries |
| Rich-club analysis | Random-network coefficient curves, node membership and edge-class comparisons | Cluster dependence of rich-club, peripheral and connecting edges |

The pipeline includes separate dataset selections and combined selections. Combined voxelwise analyses require participant selections aligned with the FA-volume order; `SubSet` selects young and healthy older participants. Dataset-specific and cross-dataset cluster locations appear in the tractography and network-analysis sections.

## Configuration and execution

Open the scripts in the MATLAB editor and use the numbered `%%` sections. Some sections consume manually loaded workspace variables or saved intermediate files. Cluster jobs are asynchronous: dependent sections require their completed output files. Dataset selections, subject ranges and density selections are specified within individual sections.

Replace every `/path/to/...` root and `YOUR_CLUSTER_PROJECT` with local locations or the cluster account. Keep trailing slashes and subfolder/file suffixes. The pipeline header defines the placeholder roots. The AudMemDem directory also contains shared cross-dataset results. Configure the function search path to include the `Functions` directory; the pipeline contains `addpath` calls to the corresponding function location.

The functions use internal TempSeq paths. `load_conn_mat` ignores its input argument. The permutation functions read the following fields:

| Function | Fields read | Saved result |
| --- | --- | --- |
| `ttest_perm_func` | `S.nperm`, `S.iteration` | `perm_clust_matrix_Age_ttest`: permutation × 3000 descending cluster sizes |
| `perm_les_func_streams` | `S.perc`, `S.nperm`, `S.iteration` | Per-participant `this_part_struct`, with iteration in the filename |
| `perm_les_func_prop` | `S.perc`, `S.nperm`, `S.iteration` | Per-participant `this_part_struct`; iteration is read but omitted from the filename |
| `rich_club_perm_new` | `S.perc`, `S.nperm`, `S.iteration` | `perm_mat`: permutation × 90 coefficient positions × participant |

`S.perc` is a percentile, so a value of `5` means the fifth percentile. Fields such as `S.path`, `S.Roi_sizes`, `S.outdir`, `S.subdir_raw`, `S.subdir_clust` and `S.partStart` do not override these functions' internal data selection or output locations. Network functions contain 78-participant allocations or loops. The pipeline's dataset job structures therefore do not, by themselves, configure these functions for a different dataset.

The paths and output formats in the function bodies must correspond to the inputs read by the selected downstream sections. The pipeline combines streamline files as 40 batches of 250, while its submission blocks request 20 batches of 500. Both total 10,000 permutations but have different batch layouts. Relative output names use the current MATLAB directory. Proportion-permutation output names depend on participant and density, so repeated iterations at the same location overwrite the same filenames.

## Matrix conventions

Raw and cluster-constrained connectivity is symmetrised by averaging reciprocal weights. `order = [1:2:90 90:-2:2]` reorders the AAL regions. `load_conn_mat` returns:

- `Raw_streams` and `Cluster_streams`: reordered connectivity weights.
- `Lesion_mats`: raw minus cluster-constrained weights.
- `Prop_Mats`: cluster-constrained divided by raw weights, elementwise.
- `Thresholding_Mat`: seed-wise `waytotal` normalisation followed by symmetrisation, in the initial AAL order.
- `ROI_sizes`: the array loaded from disk, without reordering in the loader.

Edge consistency is the between-participant standard deviation divided by the mean. Edges below the selected percentile are retained. ROI-size normalisation divides weights by the sum of the two endpoint ROI sizes. The streamline sampler uses `max(round(weight)-1,0)` as its baseline counts and draws `floor(cluster weight)` removals, rejecting negative residuals. Proportion permutations shuffle cluster/raw ratios; total removed weight can vary.

`degree` denotes weighted node strength in the network-change calculations. Its dimensions differ between scripts: ROI × participant in `Pipeline.m`, participant × ROI in `ClustperimInfAcrDens.m`. The latter fits HDI across regions within each participant; the pipeline's HDI expressions index rows of its ROI × participant matrices. In `rich_club_perm_new`, the normalisation expression uses the final `Roi_subj` from the participant loop for all connectivity slices.

The voxelwise permutation function uses a two-sided 0.01 threshold on `abs(t)` and pools both signs into a binary cluster image. The pipeline's observed-statistic block uses a 0.005-derived threshold on positive `t`. These settings describe different cluster definitions and require a consistent choice when preparing matched observed and null results.

## Dependencies and data

Dependencies include MATLAB statistical and image-processing functions, FSL, `dcm2nii`, NIfTI helpers (`load_nii`, `load_untouch_nii`, `save_nii`), Brain Connectivity Toolbox, `fdr_bh`, Leonardo/OSL helpers and the cluster interfaces (`clusterconfig`, `job2cluster`, `submit_to_cluster`). The BCT path references the `2019_03_03_BCT` directory.

Research data, participant metadata, atlas assets, helper toolboxes and saved intermediate results are external inputs. This repository contains section-based analysis code; a complete fresh-workspace run is not established by the files alone.
