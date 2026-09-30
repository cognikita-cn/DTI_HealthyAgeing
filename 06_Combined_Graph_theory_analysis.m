%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% COMBINED GRAPH THEORY ANALYSIS
%%%%%%%%%%%%%%%%% TempSeqAges + AuditMemDement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script runs both datasets as a combined graph theory analysis.
% Dataset-specific variables are retained until the final combination step.
% The combined variables keep the names used in the original pipeline so
% that the following graph theory sections can be adapted with few changes.

% Post-review setup utilities, cache construction and permutation launchers are
% kept in 06_Combined_Graph_theory_analysis_setup.m and
% 04_Mass_matched_lesion_and_permutations.m respectively.


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
    % AudMemDem-derived ageing cluster applied to TempSeqAges
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
    % TempSeqAges-derived ageing cluster applied to AudMemDem
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

%%  %%% STEP 8 %%% Analysis settings and permutation paths

dens_range = [10 20 30 40 50];
dens_tested = length(dens_range);

root_new = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';

% Mass-matched virtual lesion permutations
perm_dirTempSeq = ...
    [root_new '/TempSeq/Permutation_Combined/perm_str_'];

perm_dirAudMemDem = ...
    [root_new '/AudMemDem/Permutation_Combined/perm_str_'];

% Current rich-club nulls. These must be recomputed with the current saved
% edge masks before rich-club inference is enabled.
RC_perm_dirTempSeq = ...
    [root_new '/TempSeq_RichClub/Permutation_Combined/RC_comb_'];

RC_perm_dirAudMemDem = ...
    [root_new '/AudMemDem_RichClub/Permutation_Combined/RC_comb_'];

% Combined post-review output directory
outdir = ...
    [root_new '/Results/Combined/'];

if ~exist(outdir,'dir')
    mkdir(outdir)
end

run_virtual_lesion_analysis = 1;
run_nodewise_summary = 1;

% Keep disabled until the current rich-club permutations have been generated
% and combined using 04_Mass_matched_lesion_and_permutations.m.
run_rich_club_analysis = 0;

make_figures = 1;
save_results = 1;

%%  %%% STEP 9 %%% Investigating cluster impact across connectivity matrix densities

if run_virtual_lesion_analysis == 1
    % Participant-level permutation significance across densit ies
    Eff_acr_dens = nan(npart,dens_tested);
    Eff_acr_dens_for_plot = nan(npart,dens_tested,3);
    Clust_Coef_acr_dens = nan(npart,90,dens_tested,2);
    Deg_dep_acr_dens = nan(npart,90,dens_tested,2);
    HDI_dens = nan(npart,dens_tested);

    % Original and lesioned graph measures retained for later group analyses
    Efficiency_scores_acr_dens = nan(npart,dens_tested);
    Efficiency_scores_les_acr_dens = nan(npart,dens_tested);
    Degree_acr_dens = nan(90,npart,dens_tested);
    Degree_les_acr_dens = nan(90,npart,dens_tested);
    Clustering_coef_acr_dens = nan(90,npart,dens_tested);
    Clustering_coef_les_acr_dens = nan(90,npart,dens_tested);

    % Paired original-versus-lesioned tests in healthy participants
    Or_Les_eff = nan(4,dens_tested);
    FDR_singif_clustCoef = false(90,dens_tested);
    FDR_singif_degree = false(90,dens_tested);
    Edge_Mask_acr_dens = false(90,90,dens_tested);
    mask_data = load( ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new/combined_edge_masks.mat', ...
    'edge_masks','dens_range');

    Analysis_Selection = Healthy_Selection;

    for perc = 1:dens_tested
        Density_thresh = dens_range(perc);


        density_ind = find(mask_data.dens_range == Density_thresh,1);

        Edge_Mask = logical(mask_data.edge_masks(:,:,density_ind));
        Edge_Mask_acr_dens(:,:,perc) = Edge_Mask;


        Edge_Mask_tri = Edge_Mask&triu(true(size(Edge_Mask)),1);
        My_edge_ind = find(Edge_Mask_tri);

        % Setting up participant ROI sizes for normalisation
        Roi_Sub_mat = zeros(90,90,npart);
        for participant = 1:npart
            Roi_subj = ROI_sizes(:,participant);
            Roi_subj = Roi_subj+Roi_subj.';
            Roi_Sub_mat(:,:,participant) = Roi_subj;
        end

        % Pruning connectivity matrices
        Raw_streams_norm = Raw_streams.*Edge_Mask;
        Cluster_streams_norm = Cluster_streams.*Edge_Mask;
        Lesion_mats_norm = Lesion_mats.*Edge_Mask;

        % Normalising by the summed size of the two connected ROIs
        Raw_streams_norm = Raw_streams_norm./Roi_Sub_mat;
        Cluster_streams_norm = Cluster_streams_norm./Roi_Sub_mat;
        Lesion_mats_norm = Lesion_mats_norm./Roi_Sub_mat;

        Raw_streams_norm(~isfinite(Raw_streams_norm)) = 0;
        Cluster_streams_norm(~isfinite(Cluster_streams_norm)) = 0;
        Lesion_mats_norm(~isfinite(Lesion_mats_norm)) = 0;


        %%%%%%
        % Computing graph theory measures
        %%%%%%

        % Original connectivity matrices
        degree = zeros(90,npart);
        Efficiency_scores = zeros(npart,1);
        clustering_coef = zeros(90,npart);

        for ii = 1:npart
            M3 = Raw_streams_norm(:,:,ii);
            degree(:,ii) = sum(M3,2); % weighted degree (node strength)
            Efficiency_scores(ii,1) = efficiency_wei(M3);
            clustering_coef(:,ii) = clustering_coef_wu(M3);
        end
        disp('Measures original computed')


        % Virtual lesion connectivity matrices
        degree_les = zeros(90,npart);
        Efficiency_scores_les = zeros(npart,1);
        clustering_coef_les = zeros(90,npart);

        for ii = 1:npart
            M3 = Lesion_mats_norm(:,:,ii);
            degree_les(:,ii) = sum(M3,2);
            Efficiency_scores_les(ii,1) = efficiency_wei(M3);
            clustering_coef_les(:,ii) = clustering_coef_wu(M3);
        end
        disp('Measures lesioned computed')


        % Effect of the lesion on efficiency, clustering and weighted degree
        Eff_change = Efficiency_scores-Efficiency_scores_les;
        Perc_Eff_Change = Eff_change./Efficiency_scores;
        Perc_Eff_Change(~isfinite(Perc_Eff_Change)) = NaN;

        clust_coef_change = clustering_coef-clustering_coef_les;
        degree_change = degree-degree_les;
        degree_change_prop = degree_change./degree;
        degree_change_prop(~isfinite(degree_change_prop)) = NaN;


        % Hub disruption index: relationship between node strength and
        % proportional node loss, computed separately for each participant
        HDI = nan(npart,1);
        for participant = 1:npart
            k_raw = degree(:,participant);
            delta = degree_change_prop(:,participant);
            ValidNodes = isfinite(k_raw)&isfinite(delta);
            if sum(ValidNodes)>2 && std(k_raw(ValidNodes))>0
                b = polyfit(zscore(k_raw(ValidNodes)),delta(ValidNodes),1);
                HDI(participant,1) = b(1);
            end
        end


        % Saving graph measures across densities
        HDI_dens(:,perc) = HDI;
        Efficiency_scores_acr_dens(:,perc) = Efficiency_scores;
        Efficiency_scores_les_acr_dens(:,perc) = Efficiency_scores_les;
        Degree_acr_dens(:,:,perc) = degree;
        Degree_les_acr_dens(:,:,perc) = degree_les;
        Clustering_coef_acr_dens(:,:,perc) = clustering_coef;
        Clustering_coef_les_acr_dens(:,:,perc) = clustering_coef_les;


        % Paired tests of original versus lesioned matrices
        [~,p,~,stats] = ttest( ...
            Efficiency_scores(Analysis_Selection), ...
            Efficiency_scores_les(Analysis_Selection));
        Or_Les_eff(1,perc) = p;
        Or_Les_eff(2,perc) = stats.tstat;

        [~,p,~,stats] = ttest( ...
            mean(clustering_coef(:,Analysis_Selection),1), ...
            mean(clustering_coef_les(:,Analysis_Selection),1));
        Or_Les_eff(3,perc) = p;
        Or_Les_eff(4,perc) = stats.tstat;

        clust_pvals_curr = nan(90,1);
        degree_pvals_curr = nan(90,1);
        for ii = 1:90
            [~,clust_pvals_curr(ii)] = ttest( ...
                clustering_coef(ii,Analysis_Selection), ...
                clustering_coef_les(ii,Analysis_Selection));
            [~,degree_pvals_curr(ii)] = ttest( ...
                degree(ii,Analysis_Selection), ...
                degree_les(ii,Analysis_Selection));
        end
        clust_pvals_curr(~isfinite(clust_pvals_curr)) = 1;
        degree_pvals_curr(~isfinite(degree_pvals_curr)) = 1;
        FDR_singif_clustCoef(:,perc) = fdr_bh(clust_pvals_curr,0.05,'pdep','yes');
        FDR_singif_degree(:,perc) = fdr_bh(degree_pvals_curr,0.05,'pdep','yes');


        %%%%%%
        % Comparing the observed lesion with mass-matched permutations
        %%%%%%

        perm_fileTempSeq = [perm_dirTempSeq num2str(Density_thresh) '.mat'];
        perm_fileAudMemDem = [perm_dirAudMemDem num2str(Density_thresh) '.mat'];

        tmp_PermTempSeq = load(perm_fileTempSeq);
        tmp_PermAudMemDem = load(perm_fileAudMemDem);
        perm_strTempSeq = tmp_PermTempSeq.perc_perm_str;
        perm_strAudMemDem = tmp_PermAudMemDem.perc_perm_str;

        % Combining permutation results in the same order as Raw_streams
        this_perc_perm_eff = cat(2,perm_strTempSeq.eff,perm_strAudMemDem.eff);
        this_perc_perm_clust = cat(3,perm_strTempSeq.clust,perm_strAudMemDem.clust);
        this_perc_perm_deg = cat(3,perm_strTempSeq.deg,perm_strAudMemDem.deg);

        nperm = size(this_perc_perm_eff,1);
        assert(size(this_perc_perm_eff,2)==npart, ...
            'Efficiency permutations do not match the combined participant order.')
        assert(size(this_perc_perm_clust,3)==npart && size(this_perc_perm_deg,3)==npart, ...
            'Nodewise permutations do not match the combined participant order.')


        % Efficiency permutation significance
        % Only healthy participants have current post-review nulls.
        % The +1 correction prevents zero Monte-Carlo p-values.
        Eff_change_perm = Efficiency_scores-this_perc_perm_eff.';

        pvals = nan(npart,1);
        pvals(Healthy_Selection) = ...
            (sum( ...
                Eff_change_perm(Healthy_Selection,:) >= ...
                Eff_change(Healthy_Selection), ...
                2) + 1) ./ ...
            (nperm + 1);

        Max_perm_eff = max(this_perc_perm_eff,[],1);
        Eff_acr_dens_for_plot(:,perc,1) = Efficiency_scores;
        Eff_acr_dens_for_plot(:,perc,2) = Efficiency_scores_les;
        Eff_acr_dens_for_plot(:,perc,3) = Max_perm_eff.';


        % Clustering coefficient permutation significance
        this_perc_perm_clust = permute(this_perc_perm_clust,[1 3 2]);
        clust_coef_change_perms = clustering_coef-this_perc_perm_clust;

        pvals_clust_coef_change_more = nan(90,npart);
        pvals_clust_coef_change_less = nan(90,npart);

        pvals_clust_coef_change_more(:,Healthy_Selection) = ...
            (sum( ...
                clust_coef_change_perms(:,Healthy_Selection,:) >= ...
                clust_coef_change(:,Healthy_Selection), ...
                3) + 1) ./ ...
            (nperm + 1);

        pvals_clust_coef_change_less(:,Healthy_Selection) = ...
            (sum( ...
                clust_coef_change_perms(:,Healthy_Selection,:) <= ...
                clust_coef_change(:,Healthy_Selection), ...
                3) + 1) ./ ...
            (nperm + 1);


        % Weighted degree permutation significance
        this_perc_perm_deg = permute(this_perc_perm_deg,[1 3 2]);
        degree_dep_perm = degree-this_perc_perm_deg;

        pvals_deg_dep_change_more = nan(90,npart);
        pvals_deg_dep_change_less = nan(90,npart);

        pvals_deg_dep_change_more(:,Healthy_Selection) = ...
            (sum( ...
                degree_dep_perm(:,Healthy_Selection,:) >= ...
                degree_change(:,Healthy_Selection), ...
                3) + 1) ./ ...
            (nperm + 1);

        pvals_deg_dep_change_less(:,Healthy_Selection) = ...
            (sum( ...
                degree_dep_perm(:,Healthy_Selection,:) <= ...
                degree_change(:,Healthy_Selection), ...
                3) + 1) ./ ...
            (nperm + 1);


        % Saving permutation significance across densities
        Eff_acr_dens(:,perc) = pvals;
        Clust_Coef_acr_dens(:,:,perc,1) = pvals_clust_coef_change_more.';
        Clust_Coef_acr_dens(:,:,perc,2) = pvals_clust_coef_change_less.';
        Deg_dep_acr_dens(:,:,perc,1) = pvals_deg_dep_change_more.';
        Deg_dep_acr_dens(:,:,perc,2) = pvals_deg_dep_change_less.';

        disp(['DENSITY ' num2str(Density_thresh) '% COMPLETE'])
    end
end

%%  %%% STEP 10 %%% Nodewise effects shared across participants

if run_nodewise_summary == 1
    % Only healthy young and old participants are eligible for this part.
    % The same analysis is run for:
    % 1) young and old together, 2) young only, and 3) old only.
    Nodewise_Selection = {Healthy_Selection,Y_Selection,O_Selection};
    Nodewise_Group = {'YoungOld','Young','Old'};
    nNodewiseGroups = length(Nodewise_Selection);

    assert(any(Healthy_Selection) && any(Y_Selection) && any(O_Selection), ...
        'Young and/or old participants are missing from the combined dataset.')
    assert(~any(Healthy_Selection&D_Selection), ...
        'Dementia participants overlap with the eligible young/old sample.')

    % Column 1 = more loss than expected; column 2 = less loss than expected.
    sig_acr_dens_degree_all = false(90,2,dens_tested,nNodewiseGroups);
    sig_acr_dens_clust_all = false(90,2,dens_tested,nNodewiseGroups);
    nodewise_degree_pval = ones(90,dens_tested,nNodewiseGroups);
    nodewise_clust_pval = ones(90,dens_tested,nNodewiseGroups);
    nodewise_degree_prop_more = nan(90,dens_tested,nNodewiseGroups);
    nodewise_clust_prop_more = nan(90,dens_tested,nNodewiseGroups);

    for group = 1:nNodewiseGroups
        CurrentSelection = Nodewise_Selection{group};
        disp([Nodewise_Group{group} ' participants: ' num2str(sum(CurrentSelection))])

        for perc = 1:dens_tested
            % Weighted degree
            P_more = squeeze(Deg_dep_acr_dens(:,:,perc,1)).'<0.05;
            P_less = squeeze(Deg_dep_acr_dens(:,:,perc,2)).'<0.05;
            P_more = P_more(:,CurrentSelection);
            P_less = P_less(:,CurrentSelection);

            Eff_present = sum(P_more|P_less,2);
            more_count = sum(P_more,2);
            pval_bin = ones(90,1);

            idx = Eff_present>0;
            k = more_count(idx);
            n = Eff_present(idx);
            pval_bin(idx) = 2.*min(binocdf(k,n,0.5),1-binocdf(k-1,n,0.5));
            pval_bin = min(pval_bin,1);

            h = fdr_bh(pval_bin,0.01,'pdep','yes');
            phat = more_count./max(Eff_present,1);
            sig_acr_dens_degree_all(:,1,perc,group) = h&phat>0.5;
            sig_acr_dens_degree_all(:,2,perc,group) = h&phat<0.5;
            nodewise_degree_pval(:,perc,group) = pval_bin;
            nodewise_degree_prop_more(:,perc,group) = phat;


            % Clustering coefficient
            P_more = squeeze(Clust_Coef_acr_dens(:,:,perc,1)).'<0.05;
            P_less = squeeze(Clust_Coef_acr_dens(:,:,perc,2)).'<0.05;
            P_more = P_more(:,CurrentSelection);
            P_less = P_less(:,CurrentSelection);

            Eff_present = sum(P_more|P_less,2);
            more_count = sum(P_more,2);
            pval_bin = ones(90,1);

            idx = Eff_present>0;
            k = more_count(idx);
            n = Eff_present(idx);
            pval_bin(idx) = 2.*min(binocdf(k,n,0.5),1-binocdf(k-1,n,0.5));
            pval_bin = min(pval_bin,1);

            h = fdr_bh(pval_bin,0.01,'pdep','yes');
            phat = more_count./max(Eff_present,1);
            sig_acr_dens_clust_all(:,1,perc,group) = h&phat>0.5;
            sig_acr_dens_clust_all(:,2,perc,group) = h&phat<0.5;
            nodewise_clust_pval(:,perc,group) = pval_bin;
            nodewise_clust_prop_more(:,perc,group) = phat;
        end
    end

    % Named outputs for young and old together
    sig_acr_dens_degreeHealthy = sig_acr_dens_degree_all(:,:,:,1);
    sig_acr_dens_clustHealthy = sig_acr_dens_clust_all(:,:,:,1);

    % Named outputs for young participants only
    sig_acr_dens_degreeYoung = sig_acr_dens_degree_all(:,:,:,2);
    sig_acr_dens_clustYoung = sig_acr_dens_clust_all(:,:,:,2);

    % Named outputs for old participants only
    sig_acr_dens_degreeOld = sig_acr_dens_degree_all(:,:,:,3);
    sig_acr_dens_clustOld = sig_acr_dens_clust_all(:,:,:,3);

    % Retain the original variable names for the eligible healthy sample
    sig_acr_dens_degree = sig_acr_dens_degreeHealthy;
    sig_acr_dens_clust = sig_acr_dens_clustHealthy;
end

%%  %%% STEP 11 %%% Save combined virtual lesion results

if save_results == 1
    Combined_results = [];
    Combined_results.dens_range = dens_range;
    Combined_results.order = order;
    Combined_results.TempSeq_Selection = TempSeq_Selection;
    Combined_results.AudMemDem_Selection = AudMemDem_Selection;
    Combined_results.Y_Selection = Y_Selection;
    Combined_results.O_Selection = O_Selection;
    Combined_results.D_Selection = D_Selection;
    Combined_results.Healthy_Selection = Healthy_Selection;
    Combined_results.Age = Age;
    Combined_results.Edge_Mask_acr_dens = Edge_Mask_acr_dens;
    Combined_results.Eff_acr_dens = Eff_acr_dens;
    Combined_results.Eff_acr_dens_for_plot = Eff_acr_dens_for_plot;
    Combined_results.Clust_Coef_acr_dens = Clust_Coef_acr_dens;
    Combined_results.Deg_dep_acr_dens = Deg_dep_acr_dens;
    Combined_results.HDI_dens = HDI_dens;
    Combined_results.Efficiency_scores_acr_dens = Efficiency_scores_acr_dens;
    Combined_results.Efficiency_scores_les_acr_dens = Efficiency_scores_les_acr_dens;
    Combined_results.Degree_acr_dens = Degree_acr_dens;
    Combined_results.Degree_les_acr_dens = Degree_les_acr_dens;
    Combined_results.Clustering_coef_acr_dens = Clustering_coef_acr_dens;
    Combined_results.Clustering_coef_les_acr_dens = Clustering_coef_les_acr_dens;
    Combined_results.Or_Les_eff = Or_Les_eff;
    Combined_results.FDR_singif_clustCoef = FDR_singif_clustCoef;
    Combined_results.FDR_singif_degree = FDR_singif_degree;

    if run_nodewise_summary == 1
        Combined_results.Nodewise_Group = Nodewise_Group;
        Combined_results.Nodewise_Selection = Nodewise_Selection;
        Combined_results.sig_acr_dens_degree_all = sig_acr_dens_degree_all;
        Combined_results.sig_acr_dens_clust_all = sig_acr_dens_clust_all;
        Combined_results.nodewise_degree_pval = nodewise_degree_pval;
        Combined_results.nodewise_clust_pval = nodewise_clust_pval;
        Combined_results.nodewise_degree_prop_more = nodewise_degree_prop_more;
        Combined_results.nodewise_clust_prop_more = nodewise_clust_prop_more;
        Combined_results.sig_acr_dens_degreeHealthy = sig_acr_dens_degreeHealthy;
        Combined_results.sig_acr_dens_clustHealthy = sig_acr_dens_clustHealthy;
        Combined_results.sig_acr_dens_degreeYoung = sig_acr_dens_degreeYoung;
        Combined_results.sig_acr_dens_clustYoung = sig_acr_dens_clustYoung;
        Combined_results.sig_acr_dens_degreeOld = sig_acr_dens_degreeOld;
        Combined_results.sig_acr_dens_clustOld = sig_acr_dens_clustOld;
        Combined_results.sig_acr_dens_degree = sig_acr_dens_degree;
        Combined_results.sig_acr_dens_clust = sig_acr_dens_clust;
    end

    save([outdir 'Combined_virtual_lesion_graph_results.mat'],'Combined_results','-v7.3')
end

%%  %%% STEP 12 %%% Visualise combined virtual lesion results

if make_figures == 1
    figure
    imagesc(dens_range,1:npart,Eff_acr_dens)
    xlabel('Connectivity density (%)')
    ylabel('Participant')
    title('Efficiency change permutation p-values')
    colorbar

    figure
    imagesc(dens_range,1:npart,HDI_dens)
    xlabel('Connectivity density (%)')
    ylabel('Participant')
    title('Hub disruption index')
    colorbar

    if run_nodewise_summary == 1
        figure
        for group = 1:nNodewiseGroups
            subplot(nNodewiseGroups,2,(group-1)*2+1)
            imagesc(dens_range,1:90,squeeze(sig_acr_dens_degree_all(:,1,:,group)))
            ylabel('AAL90 node')
            title([Nodewise_Group{group} ': degree loss above expectation'])
            colorbar

            subplot(nNodewiseGroups,2,(group-1)*2+2)
            imagesc(dens_range,1:90,squeeze(sig_acr_dens_degree_all(:,2,:,group)))
            ylabel('AAL90 node')
            title([Nodewise_Group{group} ': degree loss below expectation'])
            colorbar
        end
        xlabel('Connectivity density (%)')

        figure
        for group = 1:nNodewiseGroups
            subplot(nNodewiseGroups,2,(group-1)*2+1)
            imagesc(dens_range,1:90,squeeze(sig_acr_dens_clust_all(:,1,:,group)))
            ylabel('AAL90 node')
            title([Nodewise_Group{group} ': clustering loss above expectation'])
            colorbar

            subplot(nNodewiseGroups,2,(group-1)*2+2)
            imagesc(dens_range,1:90,squeeze(sig_acr_dens_clust_all(:,2,:,group)))
            ylabel('AAL90 node')
            title([Nodewise_Group{group} ': clustering loss below expectation'])
            colorbar
        end
        xlabel('Connectivity density (%)')
    end
end

%%  %%% STEP 13 %%% Rich-club nodes and cluster dependency

if run_rich_club_analysis == 1

    rich_club_norm_acr_dens = nan(npart,90,dens_tested);
    rich_club_pval_acr_dens = ones(npart,90,dens_tested);
    degree_bin_acr_dens = nan(npart,90,dens_tested);

    rich_club_member_acr_dens = false(90,dens_tested);
    rich_club_frequency_acr_dens = nan(90,dens_tested);

    RC_threshold_acr_dens = nan(1,dens_tested);
    RC_range_acr_dens = cell(1,dens_tested);

    clust_dependency_acr_dens = nan(3,npart,dens_tested);
    RC_origpval = nan(3,dens_tested,2);

    RC_Selection = Healthy_Selection;

    for perc = 1:dens_tested

        Density_thresh = dens_range(perc);
        Edge_Mask = Edge_Mask_acr_dens(:,:,perc);

        %% Participant ROI sizes for normalisation

        Roi_Sub_mat = zeros(90,90,npart);

        for participant = 1:npart

            Roi_subj = ROI_sizes(:,participant);
            Roi_subj = Roi_subj + Roi_subj.';

            Roi_Sub_mat(:,:,participant) = Roi_subj;

        end

        Raw_streams_norm = ...
            (Raw_streams .* Edge_Mask) ./ Roi_Sub_mat;

        Cluster_streams_norm = ...
            (Cluster_streams .* Edge_Mask) ./ Roi_Sub_mat;

        Raw_streams_norm(~isfinite(Raw_streams_norm)) = 0;
        Cluster_streams_norm(~isfinite(Cluster_streams_norm)) = 0;


        %% Rich-club permutations

        tmp_RCTempSeq = load( ...
            [RC_perm_dirTempSeq num2str(Density_thresh) '.mat'], ...
            'RC_perm_comb');

        tmp_RCAudMemDem = load( ...
            [RC_perm_dirAudMemDem num2str(Density_thresh) '.mat'], ...
            'RC_perm_comb');

        RC_Perms = cat( ...
            3, ...
            tmp_RCTempSeq.RC_perm_comb, ...
            tmp_RCAudMemDem.RC_perm_comb);

        assert(size(RC_Perms,3)==npart, ...
            'Rich-club permutations do not match the combined participant order.');


        %% Participant-level rich-club coefficients

        rich_club_norm = nan(npart,90);
        rich_club_pval = ones(npart,90);
        degree_bin = zeros(npart,90);

        for ii = 1:npart

            M3 = Raw_streams_norm(:,:,ii);

            rich_club_thispart = rich_club_wu(M3);
            rich_club_thispart = rich_club_thispart(:).';

            nRich = length(rich_club_thispart);

            this_subj_perm = ...
                RC_Perms(:,1:nRich,ii);

            for kk = 1:nRich

                if ~isfinite(rich_club_thispart(kk))
                    continue
                end

                perm_k = this_subj_perm(:,kk);
                perm_k = perm_k(isfinite(perm_k));

                if isempty(perm_k)
                    continue
                end

                perm_mean = mean(perm_k);

                if perm_mean > 0
                    rich_club_norm(ii,kk) = ...
                        rich_club_thispart(kk) / perm_mean;
                end

                rich_club_pval(ii,kk) = ...
                    (sum(perm_k >= rich_club_thispart(kk)) + 1) / ...
                    (length(perm_k) + 1);

            end

            degree_bin(ii,:) = ...
                degrees_und(M3>0);

        end

        rich_club_norm_acr_dens(:,:,perc) = ...
            rich_club_norm;

        rich_club_pval_acr_dens(:,:,perc) = ...
            rich_club_pval;

        degree_bin_acr_dens(:,:,perc) = ...
            degree_bin;


        %% Rich-club threshold based on young participants

        rich_club_sig = ...
            rich_club_pval < 0.05 & ...
            rich_club_norm > 1;

        candidate_thresholds = ...
            find(mean(rich_club_sig(Y_Selection,:),1) > 0.7);

        if ~isempty(candidate_thresholds)

            %% Longest contiguous significant range

            breaks = [ ...
                0 ...
                find(diff(candidate_thresholds)>1) ...
                length(candidate_thresholds)];

            run_length = diff(breaks);

            [~,longest_run] = max(run_length);

            this_run = candidate_thresholds( ...
                breaks(longest_run)+1 : ...
                breaks(longest_run+1));

            RC_threshold = this_run(1);

            RC_threshold_acr_dens(perc) = ...
                RC_threshold;

            RC_range_acr_dens{perc} = ...
                this_run;


            %% Group rich-club node definition

            young_degree = ...
                degree_bin(Y_Selection,:);

            rich_club_frequency = ...
                mean(young_degree > RC_threshold,1);

            rich_club_member = ...
                rich_club_frequency > 0.7;

            rich_club_frequency_acr_dens(:,perc) = ...
                rich_club_frequency.';

            rich_club_member_acr_dens(:,perc) = ...
                rich_club_member.';


            %% Connection-class masks

            UppTri = triu(true(90),1);

            RR_con_mask = ...
                rich_club_member & ...
                rich_club_member.' & ...
                UppTri;

            PP_con_mask = ...
                ~rich_club_member & ...
                ~rich_club_member.' & ...
                UppTri;

            RP_con_mask = ...
                xor( ...
                    rich_club_member, ...
                    rich_club_member.') & ...
                UppTri;


            %% Cluster dependency

            clust_dependency = nan(3,npart);

            for ii = 1:npart

                Mraw = Raw_streams_norm(:,:,ii);
                Mclust = Cluster_streams_norm(:,:,ii);

                denominator = ...
                    sum(Mraw(RR_con_mask));

                if denominator > 0

                    clust_dependency(1,ii) = ...
                        sum(Mclust(RR_con_mask)) / denominator;

                end

                denominator = ...
                    sum(Mraw(PP_con_mask));

                if denominator > 0

                    clust_dependency(2,ii) = ...
                        sum(Mclust(PP_con_mask)) / denominator;

                end

                denominator = ...
                    sum(Mraw(RP_con_mask));

                if denominator > 0

                    clust_dependency(3,ii) = ...
                        sum(Mclust(RP_con_mask)) / denominator;

                end

            end

            clust_dependency_acr_dens(:,:,perc) = ...
                clust_dependency;


            %% Paired comparisons in healthy participants

            [~,p,~,stats] = ttest( ...
                clust_dependency(1,RC_Selection), ...
                clust_dependency(2,RC_Selection));

            RC_origpval(1,perc,1) = p;
            RC_origpval(1,perc,2) = stats.tstat;


            [~,p,~,stats] = ttest( ...
                clust_dependency(1,RC_Selection), ...
                clust_dependency(3,RC_Selection));

            RC_origpval(2,perc,1) = p;
            RC_origpval(2,perc,2) = stats.tstat;


            [~,p,~,stats] = ttest( ...
                clust_dependency(2,RC_Selection), ...
                clust_dependency(3,RC_Selection));

            RC_origpval(3,perc,1) = p;
            RC_origpval(3,perc,2) = stats.tstat;


            %% Summary

            fprintf( ...
                ['Density %d%% | RC range %d-%d | threshold %d | ' ...
                 'RC nodes %d | RR %d | PP %d | RP %d\n'], ...
                Density_thresh, ...
                this_run(1), ...
                this_run(end), ...
                RC_threshold, ...
                sum(rich_club_member), ...
                nnz(RR_con_mask), ...
                nnz(PP_con_mask), ...
                nnz(RP_con_mask));

        else

            warning( ...
                ['No rich-club threshold passed the young-participant ' ...
                 'criterion at density ' ...
                 num2str(Density_thresh) '%.']);

        end

        disp( ...
            ['RICH-CLUB DENSITY ' ...
             num2str(Density_thresh) ...
             '% COMPLETE'])

    end


    %% Save

    if save_results == 1

        RC_results = [];

        RC_results.dens_range = dens_range;

        RC_results.rich_club_norm_acr_dens = ...
            rich_club_norm_acr_dens;

        RC_results.rich_club_pval_acr_dens = ...
            rich_club_pval_acr_dens;

        RC_results.degree_bin_acr_dens = ...
            degree_bin_acr_dens;

        RC_results.rich_club_member_acr_dens = ...
            rich_club_member_acr_dens;

        RC_results.rich_club_frequency_acr_dens = ...
            rich_club_frequency_acr_dens;

        RC_results.RC_threshold_acr_dens = ...
            RC_threshold_acr_dens;

        RC_results.RC_range_acr_dens = ...
            RC_range_acr_dens;

        RC_results.clust_dependency_acr_dens = ...
            clust_dependency_acr_dens;

        RC_results.RC_origpval = ...
            RC_origpval;

        RC_results.RC_Selection = ...
            RC_Selection;

        save( ...
            [outdir 'Combined_rich_club_results.mat'], ...
            'RC_results', ...
            '-v7.3');

    end

end
%%

%%  %%% STEP 14 %%% Save supplementary tables
outdir = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new/Results/Combined/';
supp_dir = fullfile(outdir,'Supplementary');

if ~exist(supp_dir,'dir')
    mkdir(supp_dir)
end

root_new = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';

nperm_total = 1000;

%% Participant labels

tmp = load(fullfile(root_new,'TempSeq','permutation_input.mat'),'cache');
subjectsTempSeq = tmp.cache.subjects;

tmp = load(fullfile(root_new,'AudMemDem','permutation_input.mat'),'cache');
subjectsAudMemDem = tmp.cache.subjects;

subjects = [subjectsTempSeq(:); subjectsAudMemDem(:)];

Dataset = [ ...
    repmat({'TempSeq'},npartTempSeq,1); ...
    repmat({'AudMemDem'},npartAudMemDem,1)];

Group = repmat({''},npart,1);
Group(Y_Selection) = {'Young'};
Group(O_Selection) = {'OlderHealthy'};
Group(D_Selection) = {'Dementia'};

%% Table S1: network masks and settings

mask_data = load( ...
    fullfile(root_new,'combined_edge_masks.mat'), ...
    'edge_masks','dens_range','consistency_thresholds','n_edges');

Density = dens_range(:);
EligibleEdges = zeros(dens_tested,1);

for perc = 1:dens_tested
    EligibleEdges(perc) = ...
        nnz(triu(mask_data.edge_masks(:,:,perc),1));
end

PossibleEdges = repmat(90*89/2,dens_tested,1);
EdgeProportion = EligibleEdges ./ PossibleEdges;
Nodes = repmat(90,dens_tested,1);
PermutationsPerParticipant = repmat(nperm_total,dens_tested,1);

ConsistencyThreshold = mask_data.consistency_thresholds(:);

Table_S1 = table( ...
    Density, ...
    Nodes, ...
    EligibleEdges, ...
    PossibleEdges, ...
    EdgeProportion, ...
    ConsistencyThreshold, ...
    PermutationsPerParticipant);

writetable( ...
    Table_S1, ...
    fullfile(supp_dir,'Table_S1_network_masks_and_settings.csv'));


%% Tables S2 and S3: permutation QC and participant summaries

Dataset_S2 = {};
Density_S2 = [];
NParticipants_S2 = [];
MeanRejection_S2 = [];
MedianRejection_S2 = [];
MaxRejection_S2 = [];
MeanLesionFraction_S2 = [];
MinLesionFraction_S2 = [];
MaxLesionFraction_S2 = [];

Dataset_S3 = {};
Subject_S3 = {};
Group_S3 = {};
Age_S3 = [];
ParticipantIndex_S3 = [];
Density_S3 = [];

TotalRawUnits_S3 = [];
ObservedLesionUnits_S3 = [];
LesionFraction_S3 = [];
ObservedAffectedEdges_S3 = [];
MeanRejectionRate_S3 = [];

OriginalEfficiency_S3 = [];
LesionedEfficiency_S3 = [];
EfficiencyPermutationP_S3 = [];
HDI_S3 = [];

D(1).name = 'TempSeq';
D(1).combined_dir = fullfile(root_new,'TempSeq','Combined');
D(1).participant_ind = find(Healthy_Selection(1:npartTempSeq));
D(1).offset = 0;

D(2).name = 'AudMemDem';
D(2).combined_dir = fullfile(root_new,'AudMemDem','Combined');
D(2).participant_ind = ...
    find(Healthy_Selection(npartTempSeq+1:end));
D(2).offset = npartTempSeq;

for dd = 1:length(D)

    for perc = 1:dens_tested

        density = dens_range(perc);

        this_rejection = [];
        this_lesion_fraction = [];

        for hh = 1:length(D(dd).participant_ind)

            participant = D(dd).participant_ind(hh);
            combined_participant = participant + D(dd).offset;

            filename = fullfile( ...
                D(dd).combined_dir, ...
                ['permut_les_mat' ...
                 num2str(participant) ...
                 '_density_' ...
                 num2str(density) ...
                 '_combined.mat']);

            tmp = load(filename,'combined');
            P = tmp.combined;

            Dataset_S3{end+1,1} = D(dd).name;
            Subject_S3{end+1,1} = subjects{combined_participant};
            Group_S3{end+1,1} = Group{combined_participant};

            Age_S3(end+1,1) = Age(combined_participant);
            ParticipantIndex_S3(end+1,1) = participant;
            Density_S3(end+1,1) = density;

            TotalRawUnits_S3(end+1,1) = ...
                P.total_raw_units;

            ObservedLesionUnits_S3(end+1,1) = ...
                P.observed_lesion_units;

            LesionFraction_S3(end+1,1) = ...
                P.lesion_fraction;

            ObservedAffectedEdges_S3(end+1,1) = ...
                P.observed_affected_edges;

            MeanRejectionRate_S3(end+1,1) = ...
                P.mean_rejection_rate;

            OriginalEfficiency_S3(end+1,1) = ...
                Efficiency_scores_acr_dens(combined_participant,perc);

            LesionedEfficiency_S3(end+1,1) = ...
                Efficiency_scores_les_acr_dens(combined_participant,perc);

            EfficiencyPermutationP_S3(end+1,1) = ...
                Eff_acr_dens(combined_participant,perc);

            HDI_S3(end+1,1) = ...
                HDI_dens(combined_participant,perc);

            this_rejection(end+1,1) = ...
                P.mean_rejection_rate;

            this_lesion_fraction(end+1,1) = ...
                P.lesion_fraction;

        end

        Dataset_S2{end+1,1} = D(dd).name;
        Density_S2(end+1,1) = density;
        NParticipants_S2(end+1,1) = length(D(dd).participant_ind);

        MeanRejection_S2(end+1,1) = mean(this_rejection);
        MedianRejection_S2(end+1,1) = median(this_rejection);
        MaxRejection_S2(end+1,1) = max(this_rejection);

        MeanLesionFraction_S2(end+1,1) = mean(this_lesion_fraction);
        MinLesionFraction_S2(end+1,1) = min(this_lesion_fraction);
        MaxLesionFraction_S2(end+1,1) = max(this_lesion_fraction);

    end
end

Table_S2 = table( ...
    Dataset_S2, ...
    Density_S2, ...
    NParticipants_S2, ...
    MeanRejection_S2, ...
    MedianRejection_S2, ...
    MaxRejection_S2, ...
    MeanLesionFraction_S2, ...
    MinLesionFraction_S2, ...
    MaxLesionFraction_S2);

writetable( ...
    Table_S2, ...
    fullfile(supp_dir,'Table_S2_permutation_QC.csv'));


Table_S3 = table( ...
    Dataset_S3, ...
    Subject_S3, ...
    Group_S3, ...
    Age_S3, ...
    ParticipantIndex_S3, ...
    Density_S3, ...
    TotalRawUnits_S3, ...
    ObservedLesionUnits_S3, ...
    LesionFraction_S3, ...
    ObservedAffectedEdges_S3, ...
    MeanRejectionRate_S3, ...
    OriginalEfficiency_S3, ...
    LesionedEfficiency_S3, ...
    EfficiencyPermutationP_S3, ...
    HDI_S3);

writetable( ...
    Table_S3, ...
    fullfile(supp_dir,'Table_S3_participant_lesion_summary.csv'));


%% Tables S4-S5: rich-club results
if run_rich_club_analysis == 1

    %% Table S4: rich-club summary

    Density_S4 = dens_range(:);
    RCThreshold_S4 = RC_threshold_acr_dens(:);

    RCRange_S4 = cell(dens_tested,1);
    RCNodeIndices_S4 = cell(dens_tested,1);

    NRCNodes_S4 = zeros(dens_tested,1);
    NRREdges_S4 = zeros(dens_tested,1);
    NPPEdges_S4 = zeros(dens_tested,1);
    NRPEdges_S4 = zeros(dens_tested,1);

    RRvsPP_p = nan(dens_tested,1);
    RRvsPP_t = nan(dens_tested,1);

    RRvsRP_p = nan(dens_tested,1);
    RRvsRP_t = nan(dens_tested,1);

    PPvsRP_p = nan(dens_tested,1);
    PPvsRP_t = nan(dens_tested,1);

    for perc = 1:dens_tested

        this_range = RC_range_acr_dens{perc};

        if isempty(this_range)

            RCRange_S4{perc} = '';
            RCNodeIndices_S4{perc} = '';

        else

            RCRange_S4{perc} = ...
                strtrim(sprintf('%d ',this_range));

            rich_nodes = ...
                find(rich_club_member_acr_dens(:,perc));

            RCNodeIndices_S4{perc} = ...
                strtrim(sprintf('%d ',rich_nodes));

            NRCNodes_S4(perc) = length(rich_nodes);

            rich_member = ...
                rich_club_member_acr_dens(:,perc);

            UppTri = triu(true(90),1);

            RR = ...
                rich_member & rich_member.' & UppTri;

            PP = ...
                ~rich_member & ~rich_member.' & UppTri;

            RP = ...
                xor(rich_member,rich_member.') & UppTri;

            NRREdges_S4(perc) = nnz(RR);
            NPPEdges_S4(perc) = nnz(PP);
            NRPEdges_S4(perc) = nnz(RP);

        end

        RRvsPP_p(perc) = RC_origpval(1,perc,1);
        RRvsPP_t(perc) = RC_origpval(1,perc,2);

        RRvsRP_p(perc) = RC_origpval(2,perc,1);
        RRvsRP_t(perc) = RC_origpval(2,perc,2);

        PPvsRP_p(perc) = RC_origpval(3,perc,1);
        PPvsRP_t(perc) = RC_origpval(3,perc,2);

    end

    Table_S4 = table( ...
        Density_S4, ...
        RCThreshold_S4, ...
        RCRange_S4, ...
        NRCNodes_S4, ...
        RCNodeIndices_S4, ...
        NRREdges_S4, ...
        NPPEdges_S4, ...
        NRPEdges_S4, ...
        RRvsPP_t, ...
        RRvsPP_p, ...
        RRvsRP_t, ...
        RRvsRP_p, ...
        PPvsRP_t, ...
        PPvsRP_p);

    writetable( ...
        Table_S4, ...
        fullfile(supp_dir,'Table_S4_rich_club_summary.csv'));


    %% Table S5: participant-level rich-club dependency

    Dataset_S5 = {};
    Subject_S5 = {};
    Group_S5 = {};
    Age_S5 = [];
    Density_S5 = [];
    RCThreshold_S5 = [];

    RRDependency_S5 = [];
    PPDependency_S5 = [];
    RPDependency_S5 = [];

    for perc = 1:dens_tested

        for participant = 1:npart

            Dataset_S5{end+1,1} = Dataset{participant};
            Subject_S5{end+1,1} = subjects{participant};
            Group_S5{end+1,1} = Group{participant};

            Age_S5(end+1,1) = Age(participant);
            Density_S5(end+1,1) = dens_range(perc);

            RCThreshold_S5(end+1,1) = ...
                RC_threshold_acr_dens(perc);

            RRDependency_S5(end+1,1) = ...
                clust_dependency_acr_dens(1,participant,perc);

            PPDependency_S5(end+1,1) = ...
                clust_dependency_acr_dens(2,participant,perc);

            RPDependency_S5(end+1,1) = ...
                clust_dependency_acr_dens(3,participant,perc);

        end
    end

    Table_S5 = table( ...
        Dataset_S5, ...
        Subject_S5, ...
        Group_S5, ...
        Age_S5, ...
        Density_S5, ...
        RCThreshold_S5, ...
        RRDependency_S5, ...
        PPDependency_S5, ...
        RPDependency_S5);

    writetable( ...
        Table_S5, ...
        fullfile(supp_dir,'Table_S5_participant_rich_club_dependency.csv'));

    disp('Supplementary tables saved')

else
    disp('Rich-club supplementary tables skipped: current rich-club nulls not enabled.')
end


%%
%% ------------------------------------------------------------------------
% DISTANCE-MATCHED SENSITIVITY ANALYSIS
% -------------------------------------------------------------------------
% Uses the distance-matched permutation outputs but the same observed lesion
% and predefined combined healthy-sample edge masks.

%%  %%% STEP 8 %%% Distance-matched analysis settings

dens_range = [10 20 30 40 50];
dens_tested = length(dens_range);

root_new = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';

%% Distance-matched virtual lesion permutations

perm_dirTempSeq = ...
    [root_new '/TempSeq_Dist/Permutation_Combined/perm_str_'];

perm_dirAudMemDem = ...
    [root_new '/AudMemDem_Dist/Permutation_Combined/perm_str_'];

%% Output

outdir = ...
    [root_new '/Results/Combined_Dist/'];

if ~exist(outdir,'dir')
    mkdir(outdir)
end

run_virtual_lesion_analysis = 1;
run_nodewise_summary = 1;
make_figures = 1;
save_results = 1;


%%  %%% STEP 9 %%% Distance-matched virtual lesion analysis

if run_virtual_lesion_analysis == 1

    Eff_acr_dens = nan(npart,dens_tested);
    Eff_acr_dens_for_plot = nan(npart,dens_tested,3);

    Clust_Coef_acr_dens = ...
        nan(npart,90,dens_tested,2);

    Deg_dep_acr_dens = ...
        nan(npart,90,dens_tested,2);

    HDI_dens = nan(npart,dens_tested);

    Efficiency_scores_acr_dens = ...
        nan(npart,dens_tested);

    Efficiency_scores_les_acr_dens = ...
        nan(npart,dens_tested);

    Degree_acr_dens = ...
        nan(90,npart,dens_tested);

    Degree_les_acr_dens = ...
        nan(90,npart,dens_tested);

    Clustering_coef_acr_dens = ...
        nan(90,npart,dens_tested);

    Clustering_coef_les_acr_dens = ...
        nan(90,npart,dens_tested);

    Or_Les_eff = nan(4,dens_tested);

    FDR_singif_clustCoef = ...
        false(90,dens_tested);

    FDR_singif_degree = ...
        false(90,dens_tested);

    Edge_Mask_acr_dens = ...
        false(90,90,dens_tested);


    %% Predefined edge masks used in permutations

    mask_data = load( ...
        fullfile(root_new,'combined_edge_masks.mat'), ...
        'edge_masks','dens_range');


    %% Primary analysis sample

    Analysis_Selection = Healthy_Selection;


    for perc = 1:dens_tested

        Density_thresh = dens_range(perc);

        density_ind = ...
            find(mask_data.dens_range == Density_thresh,1);

        Edge_Mask = ...
            logical(mask_data.edge_masks(:,:,density_ind));

        Edge_Mask_acr_dens(:,:,perc) = ...
            Edge_Mask;

        Edge_Mask_tri = ...
            triu(Edge_Mask,1);

        My_edge_ind = find(Edge_Mask_tri);


        %% ROI size normalisation

        Roi_Sub_mat = zeros(90,90,npart);

        for participant = 1:npart

            Roi_subj = ...
                ROI_sizes(:,participant);

            Roi_subj = ...
                Roi_subj + Roi_subj.';

            Roi_Sub_mat(:,:,participant) = ...
                Roi_subj;

        end


        %% Apply mask

        Raw_streams_norm = ...
            Raw_streams .* Edge_Mask;

        Cluster_streams_norm = ...
            Cluster_streams .* Edge_Mask;

        Lesion_mats_norm = ...
            Lesion_mats .* Edge_Mask;


        %% ROI normalisation

        Raw_streams_norm = ...
            Raw_streams_norm ./ Roi_Sub_mat;

        Cluster_streams_norm = ...
            Cluster_streams_norm ./ Roi_Sub_mat;

        Lesion_mats_norm = ...
            Lesion_mats_norm ./ Roi_Sub_mat;

        Raw_streams_norm(~isfinite(Raw_streams_norm)) = 0;
        Cluster_streams_norm(~isfinite(Cluster_streams_norm)) = 0;
        Lesion_mats_norm(~isfinite(Lesion_mats_norm)) = 0;


        %% Original graph measures

        degree = zeros(90,npart);
        Efficiency_scores = zeros(npart,1);
        clustering_coef = zeros(90,npart);

        for ii = 1:npart

            M3 = Raw_streams_norm(:,:,ii);

            degree(:,ii) = ...
                sum(M3,2);

            Efficiency_scores(ii,1) = ...
                efficiency_wei(M3);

            clustering_coef(:,ii) = ...
                clustering_coef_wu(M3);

        end

        disp('Measures original computed')


        %% Observed lesion graph measures

        degree_les = zeros(90,npart);
        Efficiency_scores_les = zeros(npart,1);
        clustering_coef_les = zeros(90,npart);

        for ii = 1:npart

            M3 = Lesion_mats_norm(:,:,ii);

            degree_les(:,ii) = ...
                sum(M3,2);

            Efficiency_scores_les(ii,1) = ...
                efficiency_wei(M3);

            clustering_coef_les(:,ii) = ...
                clustering_coef_wu(M3);

        end

        disp('Measures lesioned computed')


        %% Observed lesion effects

        Eff_change = ...
            Efficiency_scores - Efficiency_scores_les;

        Perc_Eff_Change = ...
            Eff_change ./ Efficiency_scores;

        Perc_Eff_Change(~isfinite(Perc_Eff_Change)) = NaN;


        clust_coef_change = ...
            clustering_coef - clustering_coef_les;

        degree_change = ...
            degree - degree_les;

        degree_change_prop = ...
            degree_change ./ degree;

        degree_change_prop(~isfinite(degree_change_prop)) = NaN;


        %% Hub disruption index

        HDI = nan(npart,1);

        for participant = 1:npart

            k_raw = ...
                degree(:,participant);

            delta = ...
                degree_change_prop(:,participant);

            ValidNodes = ...
                isfinite(k_raw) & isfinite(delta);

            if sum(ValidNodes)>2 && ...
                    std(k_raw(ValidNodes))>0

                b = polyfit( ...
                    zscore(k_raw(ValidNodes)), ...
                    delta(ValidNodes), ...
                    1);

                HDI(participant,1) = b(1);

            end

        end


        %% Save graph measures across densities

        HDI_dens(:,perc) = HDI;

        Efficiency_scores_acr_dens(:,perc) = ...
            Efficiency_scores;

        Efficiency_scores_les_acr_dens(:,perc) = ...
            Efficiency_scores_les;

        Degree_acr_dens(:,:,perc) = ...
            degree;

        Degree_les_acr_dens(:,:,perc) = ...
            degree_les;

        Clustering_coef_acr_dens(:,:,perc) = ...
            clustering_coef;

        Clustering_coef_les_acr_dens(:,:,perc) = ...
            clustering_coef_les;


        %% Original versus lesioned tests
        % Healthy participants only

        [~,p,~,stats] = ttest( ...
            Efficiency_scores(Analysis_Selection), ...
            Efficiency_scores_les(Analysis_Selection));

        Or_Les_eff(1,perc) = p;
        Or_Les_eff(2,perc) = stats.tstat;


        [~,p,~,stats] = ttest( ...
            mean(clustering_coef(:,Analysis_Selection),1), ...
            mean(clustering_coef_les(:,Analysis_Selection),1));

        Or_Les_eff(3,perc) = p;
        Or_Les_eff(4,perc) = stats.tstat;


        clust_pvals_curr = nan(90,1);
        degree_pvals_curr = nan(90,1);

        for ii = 1:90

            [~,clust_pvals_curr(ii)] = ttest( ...
                clustering_coef(ii,Analysis_Selection), ...
                clustering_coef_les(ii,Analysis_Selection));

            [~,degree_pvals_curr(ii)] = ttest( ...
                degree(ii,Analysis_Selection), ...
                degree_les(ii,Analysis_Selection));

        end

        clust_pvals_curr(~isfinite(clust_pvals_curr)) = 1;
        degree_pvals_curr(~isfinite(degree_pvals_curr)) = 1;

        FDR_singif_clustCoef(:,perc) = ...
            fdr_bh(clust_pvals_curr,0.05,'pdep','yes');

        FDR_singif_degree(:,perc) = ...
            fdr_bh(degree_pvals_curr,0.05,'pdep','yes');


        %%%%%%
        % Distance-matched permutation comparison
        %%%%%%

        perm_fileTempSeq = ...
            [perm_dirTempSeq num2str(Density_thresh) '.mat'];

        perm_fileAudMemDem = ...
            [perm_dirAudMemDem num2str(Density_thresh) '.mat'];

        tmp_PermTempSeq = load(perm_fileTempSeq);
        tmp_PermAudMemDem = load(perm_fileAudMemDem);

        perm_strTempSeq = ...
            tmp_PermTempSeq.perc_perm_str;

        perm_strAudMemDem = ...
            tmp_PermAudMemDem.perc_perm_str;


        %% Combine permutations in matrix participant order

        this_perc_perm_eff = cat( ...
            2, ...
            perm_strTempSeq.eff, ...
            perm_strAudMemDem.eff);

        this_perc_perm_clust = cat( ...
            3, ...
            perm_strTempSeq.clust, ...
            perm_strAudMemDem.clust);

        this_perc_perm_deg = cat( ...
            3, ...
            perm_strTempSeq.deg, ...
            perm_strAudMemDem.deg);


        nperm = ...
            size(this_perc_perm_eff,1);

        assert( ...
            size(this_perc_perm_eff,2)==npart, ...
            'Efficiency permutations do not match participant order.')

        assert( ...
            size(this_perc_perm_clust,3)==npart && ...
            size(this_perc_perm_deg,3)==npart, ...
            'Nodewise permutations do not match participant order.')


        %% Efficiency permutation significance

        Eff_change_perm = ...
            Efficiency_scores - this_perc_perm_eff.';

        pvals = nan(npart,1);

        pvals(Healthy_Selection) = ...
            (sum( ...
                Eff_change_perm(Healthy_Selection,:) >= ...
                Eff_change(Healthy_Selection), ...
                2) + 1) ./ ...
            (nperm + 1);


        Max_perm_eff = ...
            max(this_perc_perm_eff,[],1);

        Eff_acr_dens_for_plot(:,perc,1) = ...
            Efficiency_scores;

        Eff_acr_dens_for_plot(:,perc,2) = ...
            Efficiency_scores_les;

        Eff_acr_dens_for_plot(:,perc,3) = ...
            Max_perm_eff.';


        %% Clustering coefficient permutation significance

        this_perc_perm_clust = ...
            permute(this_perc_perm_clust,[1 3 2]);

        clust_coef_change_perms = ...
            clustering_coef - this_perc_perm_clust;

        pvals_clust_coef_change_more = ...
            nan(90,npart);

        pvals_clust_coef_change_less = ...
            nan(90,npart);


        pvals_clust_coef_change_more(:,Healthy_Selection) = ...
            (sum( ...
                clust_coef_change_perms(:,Healthy_Selection,:) >= ...
                clust_coef_change(:,Healthy_Selection), ...
                3) + 1) ./ ...
            (nperm + 1);


        pvals_clust_coef_change_less(:,Healthy_Selection) = ...
            (sum( ...
                clust_coef_change_perms(:,Healthy_Selection,:) <= ...
                clust_coef_change(:,Healthy_Selection), ...
                3) + 1) ./ ...
            (nperm + 1);


        %% Weighted degree permutation significance

        this_perc_perm_deg = ...
            permute(this_perc_perm_deg,[1 3 2]);

        degree_dep_perm = ...
            degree - this_perc_perm_deg;

        pvals_deg_dep_change_more = ...
            nan(90,npart);

        pvals_deg_dep_change_less = ...
            nan(90,npart);


        pvals_deg_dep_change_more(:,Healthy_Selection) = ...
            (sum( ...
                degree_dep_perm(:,Healthy_Selection,:) >= ...
                degree_change(:,Healthy_Selection), ...
                3) + 1) ./ ...
            (nperm + 1);


        pvals_deg_dep_change_less(:,Healthy_Selection) = ...
            (sum( ...
                degree_dep_perm(:,Healthy_Selection,:) <= ...
                degree_change(:,Healthy_Selection), ...
                3) + 1) ./ ...
            (nperm + 1);


        %% Save permutation significance

        Eff_acr_dens(:,perc) = ...
            pvals;

        Clust_Coef_acr_dens(:,:,perc,1) = ...
            pvals_clust_coef_change_more.';

        Clust_Coef_acr_dens(:,:,perc,2) = ...
            pvals_clust_coef_change_less.';

        Deg_dep_acr_dens(:,:,perc,1) = ...
            pvals_deg_dep_change_more.';

        Deg_dep_acr_dens(:,:,perc,2) = ...
            pvals_deg_dep_change_less.';


        disp( ...
            ['DISTANCE-MATCHED DENSITY ' ...
             num2str(Density_thresh) ...
             '% COMPLETE'])

    end
end


%%  %%% STEP 10 %%% Nodewise effects shared across participants

if run_nodewise_summary == 1

    Nodewise_Selection = { ...
        Healthy_Selection, ...
        Y_Selection, ...
        O_Selection};

    Nodewise_Group = { ...
        'YoungOld', ...
        'Young', ...
        'Old'};

    nNodewiseGroups = ...
        length(Nodewise_Selection);

    assert( ...
        any(Healthy_Selection) && ...
        any(Y_Selection) && ...
        any(O_Selection), ...
        'Young and/or old participants are missing.')

    assert( ...
        ~any(Healthy_Selection & D_Selection), ...
        'Dementia participants overlap with healthy sample.')


    sig_acr_dens_degree_all = ...
        false(90,2,dens_tested,nNodewiseGroups);

    sig_acr_dens_clust_all = ...
        false(90,2,dens_tested,nNodewiseGroups);

    nodewise_degree_pval = ...
        ones(90,dens_tested,nNodewiseGroups);

    nodewise_clust_pval = ...
        ones(90,dens_tested,nNodewiseGroups);

    nodewise_degree_prop_more = ...
        nan(90,dens_tested,nNodewiseGroups);

    nodewise_clust_prop_more = ...
        nan(90,dens_tested,nNodewiseGroups);


    for group = 1:nNodewiseGroups

        CurrentSelection = ...
            Nodewise_Selection{group};

        disp( ...
            [Nodewise_Group{group} ...
             ' participants: ' ...
             num2str(sum(CurrentSelection))])


        for perc = 1:dens_tested

            %% Weighted degree

            P_more = ...
                squeeze(Deg_dep_acr_dens(:,:,perc,1)).' < 0.05;

            P_less = ...
                squeeze(Deg_dep_acr_dens(:,:,perc,2)).' < 0.05;

            P_more = ...
                P_more(:,CurrentSelection);

            P_less = ...
                P_less(:,CurrentSelection);


            Eff_present = ...
                sum(P_more | P_less,2);

            more_count = ...
                sum(P_more,2);

            pval_bin = ones(90,1);

            idx = Eff_present > 0;

            k = more_count(idx);
            n = Eff_present(idx);

            pval_bin(idx) = ...
                2 .* min( ...
                    binocdf(k,n,0.5), ...
                    1-binocdf(k-1,n,0.5));

            pval_bin = ...
                min(pval_bin,1);


            h = ...
                fdr_bh(pval_bin,0.01,'pdep','yes');

            phat = ...
                more_count ./ max(Eff_present,1);


            sig_acr_dens_degree_all(:,1,perc,group) = ...
                h & phat > 0.5;

            sig_acr_dens_degree_all(:,2,perc,group) = ...
                h & phat < 0.5;

            nodewise_degree_pval(:,perc,group) = ...
                pval_bin;

            nodewise_degree_prop_more(:,perc,group) = ...
                phat;


            %% Clustering coefficient

            P_more = ...
                squeeze(Clust_Coef_acr_dens(:,:,perc,1)).' < 0.05;

            P_less = ...
                squeeze(Clust_Coef_acr_dens(:,:,perc,2)).' < 0.05;

            P_more = ...
                P_more(:,CurrentSelection);

            P_less = ...
                P_less(:,CurrentSelection);


            Eff_present = ...
                sum(P_more | P_less,2);

            more_count = ...
                sum(P_more,2);

            pval_bin = ones(90,1);

            idx = ...
                Eff_present > 0;

            k = more_count(idx);
            n = Eff_present(idx);

            pval_bin(idx) = ...
                2 .* min( ...
                    binocdf(k,n,0.5), ...
                    1-binocdf(k-1,n,0.5));

            pval_bin = ...
                min(pval_bin,1);


            h = ...
                fdr_bh(pval_bin,0.01,'pdep','yes');

            phat = ...
                more_count ./ max(Eff_present,1);


            sig_acr_dens_clust_all(:,1,perc,group) = ...
                h & phat > 0.5;

            sig_acr_dens_clust_all(:,2,perc,group) = ...
                h & phat < 0.5;

            nodewise_clust_pval(:,perc,group) = ...
                pval_bin;

            nodewise_clust_prop_more(:,perc,group) = ...
                phat;

        end
    end


    sig_acr_dens_degreeHealthy = ...
        sig_acr_dens_degree_all(:,:,:,1);

    sig_acr_dens_clustHealthy = ...
        sig_acr_dens_clust_all(:,:,:,1);


    sig_acr_dens_degreeYoung = ...
        sig_acr_dens_degree_all(:,:,:,2);

    sig_acr_dens_clustYoung = ...
        sig_acr_dens_clust_all(:,:,:,2);


    sig_acr_dens_degreeOld = ...
        sig_acr_dens_degree_all(:,:,:,3);

    sig_acr_dens_clustOld = ...
        sig_acr_dens_clust_all(:,:,:,3);


    sig_acr_dens_degree = ...
        sig_acr_dens_degreeHealthy;

    sig_acr_dens_clust = ...
        sig_acr_dens_clustHealthy;

end


%%  %%% STEP 11 %%% Save distance-matched results

if save_results == 1

    Combined_results = [];

    Combined_results.null_model = ...
        'distance_matched';

    Combined_results.dens_range = ...
        dens_range;

    Combined_results.order = ...
        order;

    Combined_results.TempSeq_Selection = ...
        TempSeq_Selection;

    Combined_results.AudMemDem_Selection = ...
        AudMemDem_Selection;

    Combined_results.Y_Selection = ...
        Y_Selection;

    Combined_results.O_Selection = ...
        O_Selection;

    Combined_results.D_Selection = ...
        D_Selection;

    Combined_results.Healthy_Selection = ...
        Healthy_Selection;

    Combined_results.Age = ...
        Age;

    Combined_results.Edge_Mask_acr_dens = ...
        Edge_Mask_acr_dens;

    Combined_results.Eff_acr_dens = ...
        Eff_acr_dens;

    Combined_results.Eff_acr_dens_for_plot = ...
        Eff_acr_dens_for_plot;

    Combined_results.Clust_Coef_acr_dens = ...
        Clust_Coef_acr_dens;

    Combined_results.Deg_dep_acr_dens = ...
        Deg_dep_acr_dens;

    Combined_results.HDI_dens = ...
        HDI_dens;

    Combined_results.Efficiency_scores_acr_dens = ...
        Efficiency_scores_acr_dens;

    Combined_results.Efficiency_scores_les_acr_dens = ...
        Efficiency_scores_les_acr_dens;

    Combined_results.Degree_acr_dens = ...
        Degree_acr_dens;

    Combined_results.Degree_les_acr_dens = ...
        Degree_les_acr_dens;

    Combined_results.Clustering_coef_acr_dens = ...
        Clustering_coef_acr_dens;

    Combined_results.Clustering_coef_les_acr_dens = ...
        Clustering_coef_les_acr_dens;

    Combined_results.Or_Les_eff = ...
        Or_Les_eff;

    Combined_results.FDR_singif_clustCoef = ...
        FDR_singif_clustCoef;

    Combined_results.FDR_singif_degree = ...
        FDR_singif_degree;


    if run_nodewise_summary == 1

        Combined_results.Nodewise_Group = ...
            Nodewise_Group;

        Combined_results.Nodewise_Selection = ...
            Nodewise_Selection;

        Combined_results.sig_acr_dens_degree_all = ...
            sig_acr_dens_degree_all;

        Combined_results.sig_acr_dens_clust_all = ...
            sig_acr_dens_clust_all;

        Combined_results.nodewise_degree_pval = ...
            nodewise_degree_pval;

        Combined_results.nodewise_clust_pval = ...
            nodewise_clust_pval;

        Combined_results.nodewise_degree_prop_more = ...
            nodewise_degree_prop_more;

        Combined_results.nodewise_clust_prop_more = ...
            nodewise_clust_prop_more;

        Combined_results.sig_acr_dens_degreeHealthy = ...
            sig_acr_dens_degreeHealthy;

        Combined_results.sig_acr_dens_clustHealthy = ...
            sig_acr_dens_clustHealthy;

        Combined_results.sig_acr_dens_degreeYoung = ...
            sig_acr_dens_degreeYoung;

        Combined_results.sig_acr_dens_clustYoung = ...
            sig_acr_dens_clustYoung;

        Combined_results.sig_acr_dens_degreeOld = ...
            sig_acr_dens_degreeOld;

        Combined_results.sig_acr_dens_clustOld = ...
            sig_acr_dens_clustOld;

    end


    save( ...
        [outdir ...
         'Combined_virtual_lesion_graph_results_distance.mat'], ...
        'Combined_results', ...
        '-v7.3');

end


%%  %%% STEP 12 %%% Visualise distance-matched results

if make_figures == 1

    figure

    imagesc( ...
        dens_range, ...
        1:npart, ...
        Eff_acr_dens)

    xlabel('Connectivity density (%)')
    ylabel('Participant')
    title('Distance-matched efficiency permutation p-values')
    colorbar


    figure

    imagesc( ...
        dens_range, ...
        1:npart, ...
        HDI_dens)

    xlabel('Connectivity density (%)')
    ylabel('Participant')
    title('Hub disruption index')
    colorbar


    if run_nodewise_summary == 1

        figure

        for group = 1:nNodewiseGroups

            subplot( ...
                nNodewiseGroups, ...
                2, ...
                (group-1)*2+1)

            imagesc( ...
                dens_range, ...
                1:90, ...
                squeeze( ...
                    sig_acr_dens_degree_all(:,1,:,group)))

            ylabel('AAL90 node')

            title( ...
                [Nodewise_Group{group} ...
                 ': degree loss above distance-matched expectation'])

            colorbar


            subplot( ...
                nNodewiseGroups, ...
                2, ...
                (group-1)*2+2)

            imagesc( ...
                dens_range, ...
                1:90, ...
                squeeze( ...
                    sig_acr_dens_degree_all(:,2,:,group)))

            ylabel('AAL90 node')

            title( ...
                [Nodewise_Group{group} ...
                 ': degree loss below distance-matched expectation'])

            colorbar

        end

        xlabel('Connectivity density (%)')


        figure

        for group = 1:nNodewiseGroups

            subplot( ...
                nNodewiseGroups, ...
                2, ...
                (group-1)*2+1)

            imagesc( ...
                dens_range, ...
                1:90, ...
                squeeze( ...
                    sig_acr_dens_clust_all(:,1,:,group)))

            ylabel('AAL90 node')

            title( ...
                [Nodewise_Group{group} ...
                 ': clustering loss above distance-matched expectation'])

            colorbar


            subplot( ...
                nNodewiseGroups, ...
                2, ...
                (group-1)*2+2)

            imagesc( ...
                dens_range, ...
                1:90, ...
                squeeze( ...
                    sig_acr_dens_clust_all(:,2,:,group)))

            ylabel('AAL90 node')

            title( ...
                [Nodewise_Group{group} ...
                 ': clustering loss below distance-matched expectation'])

            colorbar

        end

        xlabel('Connectivity density (%)')

    end

end
%%


%% Distance-matched permutation supplementary table
% Participant-level QC and distance-bin summaries for the sensitivity null.

clear
clc
%%
root_new = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';

outdir = ...
    fullfile(root_new,'Results','Combined_Dist','Supplementary');

if ~exist(outdir,'dir')
    mkdir(outdir)
end

dens_range = [10 20 30 40 50];

%% Load distance-matched analysis results

tmp = load( ...
    fullfile(root_new,'Results','Combined_Dist', ...
    'Combined_virtual_lesion_graph_results_distance.mat'), ...
    'Combined_results');

R = tmp.Combined_results;

Age = R.Age;
Y_Selection = R.Y_Selection;
O_Selection = R.O_Selection;
D_Selection = R.D_Selection;
Healthy_Selection = R.Healthy_Selection;

%% Load participant information from caches

tmp = load(fullfile(root_new,'TempSeq','permutation_input.mat'),'cache');
cacheTemp = tmp.cache;

tmp = load(fullfile(root_new,'AudMemDem','permutation_input.mat'),'cache');
cacheAud = tmp.cache;

npartTempSeq = length(cacheTemp.subjects);
npartAudMemDem = length(cacheAud.subjects);

subjects = [cacheTemp.subjects(:); cacheAud.subjects(:)];

Dataset = [ ...
    repmat({'TempSeq'},npartTempSeq,1); ...
    repmat({'AudMemDem'},npartAudMemDem,1)];

Group = repmat({''},length(subjects),1);
Group(Y_Selection) = {'Young'};
Group(O_Selection) = {'OlderHealthy'};
Group(D_Selection) = {'Dementia'};

%% Dataset setup

D(1).name = 'TempSeq';
D(1).combined_dir = ...
    fullfile(root_new,'TempSeq_Dist','Combined');
D(1).participant_ind = ...
    find(cacheTemp.healthy_selection);
D(1).offset = 0;

D(2).name = 'AudMemDem';
D(2).combined_dir = ...
    fullfile(root_new,'AudMemDem_Dist','Combined');
D(2).participant_ind = ...
    find(cacheAud.healthy_selection);
D(2).offset = npartTempSeq;

%% Initialise table variables

Dataset_T = {};
Subject_T = {};
Group_T = {};

Age_T = [];
ParticipantIndex_T = [];
Density_T = [];
NPerm_T = [];

TotalRawUnits_T = [];
ObservedLesionUnits_T = [];
LesionFraction_T = [];

ObservedAffectedEdges_T = [];
MeanPermAffectedEdges_T = [];
MeanRejectionRate_T = [];

EfficiencyPermutationP_T = [];

DistanceP20_T = [];
DistanceP40_T = [];
DistanceP60_T = [];
DistanceP80_T = [];

NEdgesBin1_T = [];
NEdgesBin2_T = [];
NEdgesBin3_T = [];
NEdgesBin4_T = [];
NEdgesBin5_T = [];

LesionUnitsBin1_T = [];
LesionUnitsBin2_T = [];
LesionUnitsBin3_T = [];
LesionUnitsBin4_T = [];
LesionUnitsBin5_T = [];

ObservedAffectedBin1_T = [];
ObservedAffectedBin2_T = [];
ObservedAffectedBin3_T = [];
ObservedAffectedBin4_T = [];
ObservedAffectedBin5_T = [];

MeanPermAffectedBin1_T = [];
MeanPermAffectedBin2_T = [];
MeanPermAffectedBin3_T = [];
MeanPermAffectedBin4_T = [];
MeanPermAffectedBin5_T = [];

%% Build table

for dd = 1:length(D)

    for pp = 1:length(dens_range)

        density = dens_range(pp);

        for hh = 1:length(D(dd).participant_ind)

            participant = D(dd).participant_ind(hh);
            combined_participant = participant + D(dd).offset;

            filename = fullfile( ...
                D(dd).combined_dir, ...
                ['permut_les_dist_mat' ...
                 num2str(participant) ...
                 '_density_' ...
                 num2str(density) ...
                 '_combined.mat']);

            if ~exist(filename,'file')
                error('Missing file: %s',filename);
            end

            tmp = load(filename,'combined');
            P = tmp.combined;

            Dataset_T{end+1,1} = D(dd).name;
            Subject_T{end+1,1} = subjects{combined_participant};
            Group_T{end+1,1} = Group{combined_participant};

            Age_T(end+1,1) = Age(combined_participant);
            ParticipantIndex_T(end+1,1) = participant;
            Density_T(end+1,1) = density;
            NPerm_T(end+1,1) = P.nperm;

            TotalRawUnits_T(end+1,1) = P.total_raw_units;
            ObservedLesionUnits_T(end+1,1) = P.observed_lesion_units;
            LesionFraction_T(end+1,1) = P.lesion_fraction;

            ObservedAffectedEdges_T(end+1,1) = ...
                P.observed_affected_edges;

            MeanPermAffectedEdges_T(end+1,1) = ...
                mean(P.affected_edges);

            MeanRejectionRate_T(end+1,1) = ...
                P.mean_rejection_rate;

            EfficiencyPermutationP_T(end+1,1) = ...
                R.Eff_acr_dens(combined_participant,pp);

            %% Distance boundaries

            DistanceP20_T(end+1,1) = P.distance_percentiles(1);
            DistanceP40_T(end+1,1) = P.distance_percentiles(2);
            DistanceP60_T(end+1,1) = P.distance_percentiles(3);
            DistanceP80_T(end+1,1) = P.distance_percentiles(4);

            %% Eligible edges per distance bin

            NEdgesBin1_T(end+1,1) = P.n_edges_per_bin(1);
            NEdgesBin2_T(end+1,1) = P.n_edges_per_bin(2);
            NEdgesBin3_T(end+1,1) = P.n_edges_per_bin(3);
            NEdgesBin4_T(end+1,1) = P.n_edges_per_bin(4);
            NEdgesBin5_T(end+1,1) = P.n_edges_per_bin(5);

            %% Observed lesion mass per distance bin

            LesionUnitsBin1_T(end+1,1) = ...
                P.observed_lesion_units_per_bin(1);

            LesionUnitsBin2_T(end+1,1) = ...
                P.observed_lesion_units_per_bin(2);

            LesionUnitsBin3_T(end+1,1) = ...
                P.observed_lesion_units_per_bin(3);

            LesionUnitsBin4_T(end+1,1) = ...
                P.observed_lesion_units_per_bin(4);

            LesionUnitsBin5_T(end+1,1) = ...
                P.observed_lesion_units_per_bin(5);

            %% Observed affected edges per bin

            ObservedAffectedBin1_T(end+1,1) = ...
                P.observed_affected_edges_per_bin(1);

            ObservedAffectedBin2_T(end+1,1) = ...
                P.observed_affected_edges_per_bin(2);

            ObservedAffectedBin3_T(end+1,1) = ...
                P.observed_affected_edges_per_bin(3);

            ObservedAffectedBin4_T(end+1,1) = ...
                P.observed_affected_edges_per_bin(4);

            ObservedAffectedBin5_T(end+1,1) = ...
                P.observed_affected_edges_per_bin(5);

            %% Mean randomized affected edges per bin

            MeanPermAffectedBin1_T(end+1,1) = ...
                mean(P.affected_edges_bin(1,:));

            MeanPermAffectedBin2_T(end+1,1) = ...
                mean(P.affected_edges_bin(2,:));

            MeanPermAffectedBin3_T(end+1,1) = ...
                mean(P.affected_edges_bin(3,:));

            MeanPermAffectedBin4_T(end+1,1) = ...
                mean(P.affected_edges_bin(4,:));

            MeanPermAffectedBin5_T(end+1,1) = ...
                mean(P.affected_edges_bin(5,:));

        end
    end
end

%% Create table

Table_Dist = table( ...
    Dataset_T, ...
    Subject_T, ...
    Group_T, ...
    Age_T, ...
    ParticipantIndex_T, ...
    Density_T, ...
    NPerm_T, ...
    TotalRawUnits_T, ...
    ObservedLesionUnits_T, ...
    LesionFraction_T, ...
    ObservedAffectedEdges_T, ...
    MeanPermAffectedEdges_T, ...
    MeanRejectionRate_T, ...
    EfficiencyPermutationP_T, ...
    DistanceP20_T, ...
    DistanceP40_T, ...
    DistanceP60_T, ...
    DistanceP80_T, ...
    NEdgesBin1_T, ...
    NEdgesBin2_T, ...
    NEdgesBin3_T, ...
    NEdgesBin4_T, ...
    NEdgesBin5_T, ...
    LesionUnitsBin1_T, ...
    LesionUnitsBin2_T, ...
    LesionUnitsBin3_T, ...
    LesionUnitsBin4_T, ...
    LesionUnitsBin5_T, ...
    ObservedAffectedBin1_T, ...
    ObservedAffectedBin2_T, ...
    ObservedAffectedBin3_T, ...
    ObservedAffectedBin4_T, ...
    ObservedAffectedBin5_T, ...
    MeanPermAffectedBin1_T, ...
    MeanPermAffectedBin2_T, ...
    MeanPermAffectedBin3_T, ...
    MeanPermAffectedBin4_T, ...
    MeanPermAffectedBin5_T);

%% Save

writetable( ...
    Table_Dist, ...
    fullfile(outdir, ...
    'Table_distance_matched_permutation_summary.csv'));

save( ...
    fullfile(outdir, ...
    'Table_distance_matched_permutation_summary.mat'), ...
    'Table_Dist');

disp('Distance-matched supplementary table saved')


%%  %%% STEP 14 %%% Final combined analysis summary

disp('Combined graph theory analysis complete')
disp(['Densities analysed: ' num2str(dens_range)])
disp(['Results directory: ' outdir])
