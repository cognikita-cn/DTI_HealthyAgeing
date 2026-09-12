function output = ttest_perm_func(S)
% TTEST_PERM_FUNC Generate TempSeq voxelwise age-permutation cluster sizes.
%   output = ttest_perm_func(S)
%
%   S.nperm     Number of permutations in this batch.
%   S.iteration Batch identifier used in the output filename.
%
%   Data and output locations are set inside the function. Subject 1 is
%   excluded; subject identifiers index the fourth dimension of skeletonised
%   FA. Young and older groups use age cutoffs of <27 and >55 years.
%   One subject permutation is shared across voxels within each iteration.
%   Clusters use abs(t) above a two-sided 0.01 threshold and 26-connectivity.
%
%   Saves perm_clust_matrix_Age_ttest, a S.nperm-by-3000 matrix of descending
%   cluster sizes, zero-padded when fewer clusters are present. Column 1 is
%   the largest cluster size. Output is 1 on success or -1 on a caught error.
%   Requires NIfTI helpers, statistical/image-processing functions and the
%   Leonardo helper environment. Replace /path/to roots before execution.

try
    % Initialise the random-number generator from the current clock.
    rng(sum(100*clock));
    participant_info = readtable('/path/to/tempseq_metadata/participant_info.xlsx');

    % Initialise the Leonardo helper environment.
    pathl = '/path/to/Leonardo_FunctionsPhD'; % Leonardo helper-function directory.
    addpath(pathl);
    LBPD_startup_D(pathl);

    % Load the skeletonised FA volumes.
    % 1mm FA mask in TBSS
    fname = '/path/to/tempseq_dti/tbss/stats/all_FA_skeletonised.nii';
    skel = load_nii(fname); % Load all_FA_skeletonised file (4D file)
    skel=skel.img;
    SS = size(skel(:,:,:,1));
    %

    ind = find(skel(:,:,:,1) ~= 0); % Use nonzero voxels from the first loaded FA volume.
    [i1,i2,i3] = ind2sub(size(skel(:,:,:,1)),ind); % Convert linear voxel indices to three coordinate vectors.

    % Convert zero-valued behavioural scores to NaN.
    participant_info.WMCombined(participant_info.WMCombined==0)=NaN;
    participant_info.STMMelody(participant_info.STMMelody==0)=NaN;
    participant_info.STMRhythm(participant_info.STMRhythm==0)=NaN;

    % ValidRows = ~isnan(participant_info.WMCombined);
    % Exclude subject 1; behavioural completeness is not the active row filter.
    ValidRows = participant_info.Subject~=1;
    Valid_dat = participant_info(ValidRows,:); % Retain the rows selected by ValidRows.

    % subsetting only relevant images
    % Use subject identifiers as FA-volume indices.
    skel = skel(:,:,:,Valid_dat.Subject);
    % Normalising the memory scores
    Valid_dat.WMCombined = Valid_dat.WMCombined./max(Valid_dat.WMCombined);
    Valid_dat.STMMelody = Valid_dat.STMMelody./max(Valid_dat.STMMelody);
    Valid_dat.STMRhythm = Valid_dat.STMRhythm./max(Valid_dat.STMRhythm);
    Valid_dat.SumMemory_total = Valid_dat.WMCombined+Valid_dat.STMMelody+Valid_dat.STMRhythm;
    Valid_dat.SumMemory_Musical = Valid_dat.STMMelody+Valid_dat.STMRhythm;

    % setting degrees of freedom for future cluster thresholding
    % Set the cluster-forming threshold using the selected table-row count.
    df_age = length(Valid_dat.SumMemory_Musical) - 2;
    t_Thresh_age = tinv(1-0.01/2,df_age);
    % Define young and older groups using the age cutoffs.
    Young = Valid_dat.Subject((Valid_dat.Age < 27));
    Old = Valid_dat.Subject((Valid_dat.Age > 55));
    Y_Selection = ismember(Valid_dat.Subject,Young);
    O_Selection = ismember(Valid_dat.Subject,Old);

    Y_WM = Valid_dat.WMCombined(Y_Selection);
    Y_tot_Mem = Valid_dat.SumMemory_total(Y_Selection);
    Y_Mus_Mem = Valid_dat.SumMemory_Musical(Y_Selection);
    O_WM = Valid_dat.WMCombined(O_Selection);
    O_tot_Mem = Valid_dat.SumMemory_total(O_Selection);
    O_Mus_Mem = Valid_dat.SumMemory_Musical(O_Selection);

    All_WM = [Y_WM;O_WM];
    All_tot_Mem = [Y_tot_Mem;O_tot_Mem];
    All_Mus_Mem = [Y_Mus_Mem;O_Mus_Mem];
    All_Age = [repmat({'Young'},length(Y_WM),1);repmat({'Old'},length(O_WM),1)];

    nperm = S.nperm;
    perm_clust_matrix_Age_ttest = zeros((nperm),3000);

    % Shuffle subjects once per iteration and use that ordering at every voxel.
    for i = 1:nperm

        % Create the subject permutation shared across all voxels.
        FA = squeeze(skel(i1(1),i2(1),i3(1),:));
        perm_indexing = randperm(length(FA));

        P = zeros(length(i1),1);
        TVAL = zeros(length(i1),1);

        for nn = 1:length(ind) % Loop across voxels
            FA = squeeze(skel(i1(nn),i2(nn),i3(nn),:));
            FA_shuffled = FA(perm_indexing);
            a = FA_shuffled(Y_Selection); % take 1st group of subjects
            b = FA_shuffled(O_Selection); % take 2nd group of subjects
            % Alternative rank-sum call:
            % [p,~,stats] = ranksum(a,b);
            [~,p,~,stats] = ttest2(a,b,'Vartype','equal'); % compute t-test between the two groups
            P(nn,1) = p; % store p-values in P
            TVAL(nn,1) = stats.tstat; % store t-values in TVAL
        end
        % reshaping them into 3D matrices
        TVAL2 = zeros(SS(1),SS(2),SS(3));
        P2 = zeros(SS(1),SS(2),SS(3));
        for ii =1:length(i1)
            TVAL2(i1(ii),i2(ii),i3(ii)) = TVAL(ii);
            P2(i1(ii),i2(ii),i3(ii)) = P(ii);
        end

        conn = 26;
        Bin_TVAL_AgeONLY3D = zeros(size((TVAL2)));
        % Pool positive and negative suprathreshold voxels into one binary image.
        Bin_TVAL_AgeONLY3D((abs(TVAL2)>t_Thresh_age)) = 1;
        CC_Age = bwconncomp(Bin_TVAL_AgeONLY3D,conn); % Identify connected suprathreshold voxel clusters.

        a_Age = CC_Age.PixelIdxList;
        lengths_Age = cellfun(@length,a_Age);
        b_Age = sort(lengths_Age,'descend');

        % Store up to 3000 cluster sizes, in descending order.
        if length(b_Age) > 3000
            perm_clust_matrix_Age_ttest(i,:) = b_Age(1:3000);
        else
            perm_clust_matrix_Age_ttest(i,1:length(b_Age)) = b_Age;
        end

    end

    this_iteration = S.iteration;

    % Age only effect
    outdir = '/path/to/tempseq_dti/tbss/stats/perm_res_final/AgeOnlyttest_full';
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end

    filename_Age = fullfile(outdir,['Permutation_res_ttest_Age_Clusters_iteration_' num2str(this_iteration) '.mat']);
    % Save the complete batch; column 1 contains the maximum cluster size.
    save(filename_Age,'perm_clust_matrix_Age_ttest');

    output = 1;
% Report a caught error and return the failure status.
catch ME
    fprintf('ERROR in Perm_function_ANOVA_continious: %s/m',ME.message);
    output = -1;
end

end
