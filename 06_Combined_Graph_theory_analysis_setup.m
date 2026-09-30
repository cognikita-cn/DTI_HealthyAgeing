%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% COMBINED GRAPH THEORY ANALYSIS
%%%%%%%%%%%%%%%%% TempSeqAges + AuditMemDement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script sets up both datasets for a combined graph theory analysis.
% Dataset-specific variables are retained until the final combination step.
% The combined variables keep the names used in the original pipeline so
% that the following graph theory sections can be adapted with few changes.

clear all

%%  %%% STEP 1 %%% General setup

addpath('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/Functions')
addpath('/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/Graph_Theory/BCT/2019_03_03_BCT')

% Select which white matter cluster is used for tractography matrices:
% 0 = corresponding/within-dataset cluster
% 1 = cross-dataset cluster (FINAL POST-REVIEW ANALYSIS)
% 2 = combined cluster
MNI_clust = 1;

% AAL90 ordering used throughout the original graph theory pipeline
order = [1:2:90 90:-2:2];

%%  %%% STEP 2 %%% Set up paths for both datasets

% TempSeqAges
S_TempSeq = [];
S_TempSeq.path = '/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/00*';
S_TempSeq.Roi_sizes = '/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/Graph_Theory_Age_Diff/ROI_sizes.mat';
S_TempSeq.subdir_raw = '/Tractography_AgeDiff/NoCluster/Age_ttest/AAL90_5000stream/';

if MNI_clust == 0
    S_TempSeq.subdir_clust = '/Tractography_AgeDiff_005/AAL_Clust_full/Age_ttest/AAL90_5000stream/';
elseif MNI_clust == 1
    S_TempSeq.subdir_clust = '/Tractography_CrossCheckOY_Fin/AAL_Clust_full/Age_ttest/AAL90_5000stream/';
else
    S_TempSeq.subdir_clust = '/Tractography_AgeDiff/WM_Cluster_Combined/Age_ttest/AAL90_5000stream/';
end

% AuditMemDement
S_AudMemDem = [];
S_AudMemDem.path = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/0*';
S_AudMemDem.Roi_sizes = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing/ROI_sizes.mat';
S_AudMemDem.subdir_raw = '/Tractography_AgeDiff/NoCluster/Age_ttest/AAL90_5000stream/';

if MNI_clust == 0
    S_AudMemDem.subdir_clust = '/Tractography_AgeDiff/WM_Cluster_ThisDat/Age_ttest/AAL90_5000stream/';
elseif MNI_clust == 1
    S_AudMemDem.subdir_clust = '/Tractography_AgeDiff/WM_Cluster_CrossCheck/Age_ttest/AAL90_5000stream/';
else
    S_AudMemDem.subdir_clust = '/Tractography_AgeDiff/WM_Cluster_Combined/Age_ttest/AAL90_5000stream/';
end

%%  %%% STEP 3 %%% Load participant information and define groups

% TempSeqAges participant information
participant_infoTempSeq = readtable('/aux/MINDLAB2021_MEG-TempSeqAges/Nikita/TSA2021_Nikita.xlsx');

% Converting 0 to NaN as in the original pipeline
participant_infoTempSeq.WMCombined(participant_infoTempSeq.WMCombined==0) = NaN;

ValidRowsTempSeq = ~isnan(participant_infoTempSeq.WMCombined); % for later WM relations
%Valid_datTempSeq = participant_infoTempSeq(ValidRowsTempSeq,:);
Valid_datTempSeq = participant_infoTempSeq(:,:); % all participants, as in the original graph theory setup

YoungTempSeq = Valid_datTempSeq.Subject(Valid_datTempSeq.Age < 27);
OldTempSeq = Valid_datTempSeq.Subject(Valid_datTempSeq.Age > 55);
Y_SelectionTempSeq = ismember(Valid_datTempSeq.Subject,YoungTempSeq);
O_SelectionTempSeq = ismember(Valid_datTempSeq.Subject,OldTempSeq);


% AuditMemDement participant information
participant_infoAudMemDem = readtable('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/tbss/Descript_Dat/Participant_categories_AuditMemDement_ages.xlsx');

path_outAudMemDem = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/';
path5AudMemDem = dir([path_outAudMemDem '0*']);
DTI_AudMemDem = zeros(length(participant_infoAudMemDem.Group),1);

for gg = 1:length(path5AudMemDem)
    path6AudMemDem = dir([path5AudMemDem(gg).folder '/' path5AudMemDem(gg).name '/20*']);
    if ~isempty(path6AudMemDem)
        DTI_AudMemDem(gg,1) = 1;
    end
end

participant_infoAudMemDem.DTI = DTI_AudMemDem;
Valid_participant_infoAudMemDem = participant_infoAudMemDem(participant_infoAudMemDem.DTI==1,:);

Y_SelectionAudMemDem = ismember(Valid_participant_infoAudMemDem.Group,{'A'});
O_SelectionAudMemDem = ismember(Valid_participant_infoAudMemDem.Group,{'B'});
%O_SelectionAudMemDem = ismember(Valid_participant_infoAudMemDem.Group,{'B','B_MCI'});
D_SelectionAudMemDem = ismember(Valid_participant_infoAudMemDem.Group,{'C'});

%%  %%% STEP 4 %%% Load connectivity matrices and ROI sizes

% The same loading function is used for both datasets in the original code
mat_structTempSeq = load_conn_mat_dem(S_TempSeq);
mat_structAudMemDem = load_conn_mat_dem(S_AudMemDem);


% TempSeqAges matrices
Thresholding_MatTempSeq = mat_structTempSeq.Thresholding_Mat;
Raw_streamsTempSeq = mat_structTempSeq.Raw_streams;
Cluster_streamsTempSeq = mat_structTempSeq.Cluster_streams;
Lesion_matsTempSeq = mat_structTempSeq.Lesion_mats;
Prop_MatsTempSeq = mat_structTempSeq.Prop_Mats;
ROI_sizesTempSeq = mat_structTempSeq.ROI_sizes;
npartTempSeq = size(Raw_streamsTempSeq,3);


% AuditMemDement matrices
Thresholding_MatAudMemDem = mat_structAudMemDem.Thresholding_Mat;
Raw_streamsAudMemDem = mat_structAudMemDem.Raw_streams;
Cluster_streamsAudMemDem = mat_structAudMemDem.Cluster_streams;
Lesion_matsAudMemDem = mat_structAudMemDem.Lesion_mats;
Prop_MatsAudMemDem = mat_structAudMemDem.Prop_Mats;
ROI_sizesAudMemDem = mat_structAudMemDem.ROI_sizes;
npartAudMemDem = size(Raw_streamsAudMemDem,3);

%%  %%% STEP 5 %%% Check participant and matrix alignment

% These checks catch dimension and participant-count mismatches before the
% datasets are combined. Participant ordering should still be confirmed
% against the subject order returned by load_conn_mat_dem.
assert(size(Raw_streamsTempSeq,1)==90 && size(Raw_streamsTempSeq,2)==90, ...
    'TempSeqAges connectivity matrices are not 90 by 90.')
assert(size(Raw_streamsAudMemDem,1)==90 && size(Raw_streamsAudMemDem,2)==90, ...
    'AuditMemDement connectivity matrices are not 90 by 90.')

assert(length(Y_SelectionTempSeq)==npartTempSeq, ...
    'TempSeqAges participant information does not match the connectivity matrices.')
assert(length(Y_SelectionAudMemDem)==npartAudMemDem, ...
    'AuditMemDement participant information does not match the connectivity matrices.')

assert(size(Thresholding_MatTempSeq,3)==npartTempSeq && ...
    size(Cluster_streamsTempSeq,3)==npartTempSeq && ...
    size(Lesion_matsTempSeq,3)==npartTempSeq && ...
    size(Prop_MatsTempSeq,3)==npartTempSeq, ...
    'TempSeqAges matrices do not contain the same number of participants.')
assert(size(Thresholding_MatAudMemDem,3)==npartAudMemDem && ...
    size(Cluster_streamsAudMemDem,3)==npartAudMemDem && ...
    size(Lesion_matsAudMemDem,3)==npartAudMemDem && ...
    size(Prop_MatsAudMemDem,3)==npartAudMemDem, ...
    'AuditMemDement matrices do not contain the same number of participants.')

assert(size(ROI_sizesTempSeq,2)==npartTempSeq, ...
    'TempSeqAges ROI sizes do not match the connectivity matrices.')
assert(size(ROI_sizesAudMemDem,2)==npartAudMemDem, ...
    'AuditMemDement ROI sizes do not match the connectivity matrices.')

%%  %%% STEP 6 %%% Combine both datasets

% Combined matrices. These retain the original generic variable names so
% they can be used directly in the following graph theory sections.
Thresholding_Mat = cat(3,Thresholding_MatTempSeq,Thresholding_MatAudMemDem);
Raw_streams = cat(3,Raw_streamsTempSeq,Raw_streamsAudMemDem);
Cluster_streams = cat(3,Cluster_streamsTempSeq,Cluster_streamsAudMemDem);
Lesion_mats = cat(3,Lesion_matsTempSeq,Lesion_matsAudMemDem);
Prop_Mats = cat(3,Prop_MatsTempSeq,Prop_MatsAudMemDem);
ROI_sizes = cat(2,ROI_sizesTempSeq,ROI_sizesAudMemDem);

npart = size(Raw_streams,3);


% Dataset indicators in the same participant order as the combined matrices
TempSeq_Selection = [true(npartTempSeq,1); false(npartAudMemDem,1)];
AudMemDem_Selection = [false(npartTempSeq,1); true(npartAudMemDem,1)];


% Combined group indicators
Y_Selection = [Y_SelectionTempSeq(:); Y_SelectionAudMemDem(:)];
O_Selection = [O_SelectionTempSeq(:); O_SelectionAudMemDem(:)];
D_Selection = [false(npartTempSeq,1); D_SelectionAudMemDem(:)];

Healthy_Selection = Y_Selection|O_Selection;
All_Selection = true(npart,1);


% Combined age vector for later age-related analyses
Age = [Valid_datTempSeq.Age; Valid_participant_infoAudMemDem.Age];


% Keep a combined structure in the same format as the original loader output
mat_struct = [];
mat_struct.Thresholding_Mat = Thresholding_Mat;
mat_struct.Raw_streams = Raw_streams;
mat_struct.Cluster_streams = Cluster_streams;
mat_struct.Lesion_mats = Lesion_mats;
mat_struct.Prop_Mats = Prop_Mats;
mat_struct.ROI_sizes = ROI_sizes;

%%  %%% STEP 7 %%% Final setup checks and summary

assert(npart==npartTempSeq+npartAudMemDem, ...
    'The number of participants changed while combining datasets.')
assert(length(Age)==npart, ...
    'The combined age vector does not match the combined matrices.')
assert(length(Y_Selection)==npart && length(O_Selection)==npart && length(D_Selection)==npart, ...
    'The combined group vectors do not match the combined matrices.')

disp('Combined graph theory setup complete')
disp(['TempSeqAges participants: ' num2str(npartTempSeq)])
disp(['AuditMemDement participants: ' num2str(npartAudMemDem)])
disp(['Total participants: ' num2str(npart)])
disp(['Young participants: ' num2str(sum(Y_Selection))])
disp(['Older healthy participants: ' num2str(sum(O_Selection))])
disp(['Dementia participants: ' num2str(sum(D_Selection))])


%%  %%% STEP 8 %%% Build combined healthy-sample edge masks
% The edge mask is derived once from healthy participants in both datasets
% and saved for reuse by the post-review permutation and graph analyses.

%% Combined healthy-sample edge masks

order = [1:2:90 90:-2:2];
dens_range = [10 20 30 40 50];

Healthy_SelectionTempSeq = ...
    Y_SelectionTempSeq | O_SelectionTempSeq;

Healthy_SelectionAudMemDem = ...
    Y_SelectionAudMemDem | O_SelectionAudMemDem;

Thresholding_MatTempSeq_Healthy = ...
    Thresholding_MatTempSeq(:,:,Healthy_SelectionTempSeq);

Thresholding_MatAudMemDem_Healthy = ...
    Thresholding_MatAudMemDem(:,:,Healthy_SelectionAudMemDem);

Thresholding_Mat = cat(3, ...
    Thresholding_MatTempSeq_Healthy, ...
    Thresholding_MatAudMemDem_Healthy);

%% Consistency matrix

sd_mat = std(Thresholding_Mat,0,3);
mean_mat = mean(Thresholding_Mat,3);

Consistency_mat = sd_mat ./ mean_mat;
Consistency_mat(~isfinite(Consistency_mat)) = Inf;

upper_ind = triu(true(90),1);
Consistency_values = Consistency_mat(upper_ind);

%% Edge masks

edge_masks = false(90,90,length(dens_range));
consistency_thresholds = nan(length(dens_range),1);
n_edges = zeros(length(dens_range),1);

for ii = 1:length(dens_range)

    Density_thresh = dens_range(ii);

    Consistency_threshold = ...
        prctile(Consistency_values,Density_thresh);

    Edge_Mask = ...
        Consistency_mat < Consistency_threshold;

    Edge_Mask = Edge_Mask | Edge_Mask.';
    Edge_Mask(1:91:end) = false;

    Edge_Mask = Edge_Mask(order,order);

    edge_masks(:,:,ii) = Edge_Mask;
    consistency_thresholds(ii) = Consistency_threshold;
    n_edges(ii) = nnz(triu(Edge_Mask,1));

end

%% Save

outdir = ...
    ['/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/' ...
     'GraphTheory_Ageing_Cross_new'];

if ~exist(outdir,'dir')
    mkdir(outdir);
end

n_healthy_TempSeq = sum(Healthy_SelectionTempSeq);
n_healthy_AudMemDem = sum(Healthy_SelectionAudMemDem);
n_healthy_total = n_healthy_TempSeq + n_healthy_AudMemDem;

save( ...
    fullfile(outdir,'combined_edge_masks.mat'), ...
    'edge_masks', ...
    'dens_range', ...
    'consistency_thresholds', ...
    'n_edges', ...
    'order', ...
    'Healthy_SelectionTempSeq', ...
    'Healthy_SelectionAudMemDem', ...
    'n_healthy_TempSeq', ...
    'n_healthy_AudMemDem', ...
    'n_healthy_total');

%% Summary

fprintf('Healthy TempSeq participants: %d\n',n_healthy_TempSeq);
fprintf('Healthy AudMemDem participants: %d\n',n_healthy_AudMemDem);
fprintf('Combined healthy sample: %d\n\n',n_healthy_total);

disp(table( ...
    dens_range(:), ...
    consistency_thresholds, ...
    n_edges, ...
    'VariableNames', ...
    {'Density','ConsistencyThreshold','NumberOfEdges'}));

%%  %%% STEP 9 %%% Check participant-specific AAL centroid distances
% One-participant visual QC for the distance matrix used by the distance-matched null.

%%
addpath('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/Functions')
masks_file = ...
    ['/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/' ...
     '0001/masks.txt'];

[D,C] = get_aal_centroid_distances(masks_file);

figure
imagesc(D)
axis square
colorbar
title('AAL centroid distance (mm)')

%%  %%% STEP 10 %%% Build final permutation caches
% Uses the participant selections and final cross-dataset cluster paths defined
% above. Run Steps 1-10 in the same MATLAB session.
%
% get_aal_centroid_distances requires SPM12. The explicit check below avoids
% silently relying on a placeholder or an unknown MATLAB path.

if isempty(which('spm_vol')) || isempty(which('spm_read_vols'))
    error(['SPM12 is not on the MATLAB path. Add SPM12 before building ' ...
           'the permutation caches.']);
end

root_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';

Healthy_SelectionTempSeq = ...
    Y_SelectionTempSeq | O_SelectionTempSeq;

Healthy_SelectionAudMemDem = ...
    Y_SelectionAudMemDem | O_SelectionAudMemDem;

assert(length(Healthy_SelectionTempSeq)==npartTempSeq, ...
    'TempSeq healthy selection does not match the valid-DTI matrix order.');

assert(length(Healthy_SelectionAudMemDem)==npartAudMemDem, ...
    'AudMemDem healthy selection does not match the valid-DTI matrix order.');


%% TempSeq cache

S = struct();

S.path = ...
    '/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/00*';

S.Roi_sizes = ...
    ['/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/' ...
     'Graph_Theory_Age_Diff/ROI_sizes.mat'];

S.subdir_raw = ...
    'Tractography_AgeDiff/NoCluster/Age_ttest/AAL90_5000stream';

S.subdir_clust = ...
    'Tractography_CrossCheckOY_Fin/AAL_Clust_full/Age_ttest/AAL90_5000stream';

fprintf('\nBuilding TempSeq cache...\n');

cache = make_permutation_cache( ...
    S, ...
    Healthy_SelectionTempSeq);

tempseq_dir = fullfile(root_dir,'TempSeq');

if ~exist(tempseq_dir,'dir')
    mkdir(tempseq_dir);
end

save( ...
    fullfile(tempseq_dir,'permutation_input.mat'), ...
    'cache', ...
    '-v7.3');

clear cache S


%% AudMemDem cache

S = struct();

S.path = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/0*';

S.Roi_sizes = ...
    ['/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/' ...
     'GraphTheory_Ageing/ROI_sizes.mat'];

S.subdir_raw = ...
    'Tractography_AgeDiff/NoCluster/Age_ttest/AAL90_5000stream';

S.subdir_clust = ...
    'Tractography_AgeDiff/WM_Cluster_CrossCheck/Age_ttest/AAL90_5000stream';

fprintf('\nBuilding AudMemDem cache...\n');

cache = make_permutation_cache( ...
    S, ...
    Healthy_SelectionAudMemDem);

audmem_dir = fullfile(root_dir,'AudMemDem');

if ~exist(audmem_dir,'dir')
    mkdir(audmem_dir);
end

save( ...
    fullfile(audmem_dir,'permutation_input.mat'), ...
    'cache', ...
    '-v7.3');

clear cache S

fprintf('\nFinal permutation caches built.\n');
