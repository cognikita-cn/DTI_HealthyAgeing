function output = perm_les_func_streams(S)
% PERM_LES_FUNC_STREAMS Generate multinomial streamline-lesion permutations.
%   output = perm_les_func_streams(S)
%
%   S.perc      Percentile of the edge-consistency distribution (not an index).
%   S.nperm     Accepted permutations per participant.
%   S.iteration Batch identifier used in output filenames.
%
%   Uses the TempSeq paths specified below and loops over 78 participants.
%   Other S fields, including path, Roi_sizes, outdir, subdir_raw,
%   subdir_clust and partStart, are not read by this implementation.
%   Matrices are symmetrised and reordered to 90 AAL regions. Edge selection
%   uses the coefficient of variation of waytotal-normalised connectivity.
%   Streamline counts are max(round(weight)-1,0). The multinomial draw uses
%   floor(cluster weight) removals; draws with negative residuals are rejected.
%
%   Saves one this_part_struct per participant and batch: efficiency
%   (S.nperm-by-1), clustering_coef and degree_les (90-by-S.nperm).
%   degree_les contains node strength, the sum of remaining edge weights.
%   Output is 1 on success or -1 on a caught error. Requires Brain Connectivity
%   Toolbox and mnrnd. Replace /path/to roots before execution.

try
    % Initialise the random-number generator from the current clock.
    rng(sum(100*clock));

    my_prctile = S.perc;
    order = [1:2:90 90:-2:2];
    Cluster_dir = '/path/to/tempseq_dti/Graph_Theory_Age_Diff/ThroughClust';
    Cluster_dir = dir(Cluster_dir);

    No_Cluster_dir = '/path/to/tempseq_dti/Graph_Theory_Age_Diff/NoCluster';
    No_Cluster_dir = dir(No_Cluster_dir);

    path = dir('/path/to/tempseq_dti/00*'); % path to subjects
    % Allocate ROI-by-ROI-by-participant connectivity arrays.
    Thresholding_Mat = zeros(90,90,78);
    Raw_streams = zeros(90,90,78);
    Prop_Mats = zeros(90,90,78);
    Cluster_streams = zeros(90,90,78);
    order = [1:2:90 90:-2:2];
    for ii = 1:length(Cluster_dir)-2
        subj = ii+2;
        subj_num = Cluster_dir(subj).name;
        % Parse the subject number used to index the subject directory listing.
        subj_num = str2double(subj_num(6:9));
        M = readtable([path(subj_num).folder '/' path(subj_num).name '/Tractography_AgeDiff/NoCluster/Age_ttest/AAL90_5000stream/fdt_network_matrix']);
        M = table2array(M);
        M = M(1:90,1:90);
        % Average reciprocal raw connections before reordering the AAL regions.
        M_raw = (M+M.')/2;
        M_raw = M_raw(order,order);

        Raw_streams(:,:,ii) = M_raw;

        M_Clust = readtable([path(subj_num).folder '/' path(subj_num).name '/Tractography_AgeDiff_005/AAL_Clust/Age_ttest/AAL90_5000stream/fdt_network_matrix']);
        M_Clust = table2array(M_Clust);
        M_Clust = M_Clust(1:90,1:90);
        M_Clust = (M_Clust+M_Clust.')/2;
        M_Clust = M_Clust (order,order);
        Cluster_streams(:,:,ii) = M_Clust;
        a = readtable([path(subj_num).folder '/' path(subj_num).name '/Tractography_AgeDiff/NoCluster/Age_ttest/AAL90_5000stream/waytotal']);
        a = table2array(a);
        a = a(1:90,1);
        % Normalise each seed row by its waytotal before symmetrising.
        M = M./a;
        M = (M+M.')/2;
        Thresholding_Mat(:,:,ii) = M;
    end

    Density_thresh = my_prctile;

    % Select consistent edges using between-participant SD divided by mean.
    sd_mat = std(Thresholding_Mat,0,3);
    mean_mat = mean(Thresholding_Mat,3);
    Consistecy_mat = sd_mat./mean_mat;

    Edge_Mask = Consistecy_mat < prctile(Consistecy_mat(:),Density_thresh);
    Edge_Mask = Edge_Mask(order,order);
    % Restrict sampling to the upper triangle, including any selected diagonal entries.
    Edge_Mask_tri = Edge_Mask&~tril(true(size(Edge_Mask)),-1);
    My_edge_ind = find(Edge_Mask_tri);

    % Sum cluster-constrained weights across selected upper-triangle edges.
    clust_sums = zeros(78,1);
    for ii = 1:length(Cluster_dir)-2
        a = Cluster_streams(:,:,ii);
        clust_sums(ii,1) = sum(a(Edge_Mask_tri));
    end

    outdir = ['/path/to/tempseq_dti/Graph_Theory_Age_Diff/permuts_Lesion_streamlines_005/density_' num2str(my_prctile)]
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end

    tmp = load('/path/to/tempseq_dti/Graph_Theory_Age_Diff/ROI_sizes.mat')
    ROI_sizes = tmp.ROI_sizes;
    this_iteration = S.iteration;

    nperm = S.nperm;
    % Generate and save a separate permutation batch for each participant.
    for participant = 1:78
        efficiency = zeros(nperm,1);
        clustering_coef = zeros(90,nperm);
        degree_les = zeros(90,nperm);
        Roi_subj = ROI_sizes(:,participant);
        % Form pairwise sums of ROI sizes for edge-weight normalisation.
        Roi_subj = Roi_subj+Roi_subj.';
        Roi_subj(~Edge_Mask_tri) = 0;
        this_raw_mat = Raw_streams(:,:,participant);
        this_raw_mat(~Edge_Mask_tri) = 0;
        % Convert symmetrised weights to the integer counts used by the sampler.
        raw_mat_counter = max(round(this_raw_mat)-1,0);
        total_streams = sum(raw_mat_counter(Edge_Mask_tri));
        % Use the floored cluster-weight sum as the multinomial removal count.
        this_lesion = floor(clust_sums(participant,1));

        raw_vals = raw_mat_counter(My_edge_ind);
        Multinonm_prob = raw_vals/total_streams;
        Multinonm_prob = Multinonm_prob/sum(Multinonm_prob);
        Multinonm_prob = Multinonm_prob';
        permuting = 1;
        perm = 1;
        % Accept draws until S.nperm nonnegative residual networks are obtained.
        while permuting == 1
            % Draw removal counts with probabilities proportional to integer edge counts.
            this_les = mnrnd(this_lesion,Multinonm_prob);
            this_les = raw_vals - this_les';
            if any(this_les<0), continue; end
            Lesion_perm = this_raw_mat;
            Lesion_perm(My_edge_ind) = this_les;

            % Normalise residual counts by pairwise ROI sizes.
            Lesion_perm = Lesion_perm./Roi_subj;
            % Mirror the upper triangle to restore a symmetric connectivity matrix.
            Lesion_perm = triu(Lesion_perm)+triu(Lesion_perm,1)';
            Lesion_perm(find(isnan(Lesion_perm))) = 0;
            efficiency(perm,1) = efficiency_wei(Lesion_perm);
            clustering_coef(:,perm) = clustering_coef_wu(Lesion_perm);
            degree_les(:,perm) = sum(Lesion_perm,2);
            perm = perm+1;
            if perm == nperm+1
                permuting = 0;
            end
        end
        % Collect the permutation summaries for this participant.
        this_part_struct = [];
        this_part_struct.efficiency = efficiency;
        this_part_struct.clustering_coef = clustering_coef;
        this_part_struct.degree_les = degree_les;
        filename = fullfile(outdir,['permut_les_mat' num2str(participant) '_iter_' num2str(this_iteration) '.mat']);
        save(filename,'this_part_struct');
    end

    output = 1;
% Report a caught error and return the failure status.
catch ME
    fprintf('ERROR in Perm_function_ANOVA_continious: %s/m',ME.message);
    output = -1;
end

end
