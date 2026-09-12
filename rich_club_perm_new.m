function output = rich_club_perm_new(S)
% RICH_CLUB_PERM_NEW Generate random-network rich-club coefficient curves.
%   output = rich_club_perm_new(S)
%
%   S.perc      Percentile of the edge-consistency distribution (not an index).
%   S.nperm     Random networks per participant.
%   S.iteration Batch identifier used in the output filename.
%
%   Uses the TempSeq paths specified below. Connectivity and ROI-size arrays
%   are prepared for 78 participants; saved output follows the NoCluster
%   subject list. Additional dataset/path fields in S are not read.
%   Each random network is generated with randmio_und_connected(M3,25),
%   followed by rich_club_wu. Short coefficient curves are zero-padded to 90.
%
%   Saves perm_mat with dimensions S.nperm-by-90-by-number of listed subjects.
%   Output is 1 on success or -1 on a caught error. Requires Brain Connectivity
%   Toolbox. Replace /path/to roots before execution.

try
    % Initialise the random-number generator from the current clock.
    rng(sum(100*clock));
    my_prctile = S.perc;
    order = [1:2:90 90:-2:2];
    Cluster_dir = '/path/to/tempseq_dti/Graph_Theory_Age_Diff/ThroughClust';
    Cluster_dir = dir(Cluster_dir);

    No_Cluster_dir = '/path/to/tempseq_dti/Graph_Theory_Age_Diff/NoCluster';
    addpath('/path/to/BCT/2019_03_03_BCT'); % add BCT functions to path (FFM)
    No_Cluster_dir = dir(No_Cluster_dir);

    path = dir('/path/to/tempseq_dti/00*'); % path to subjects
    % Allocate ROI-by-ROI-by-participant connectivity arrays.
    Thresholding_Mat = zeros(90,90,78);
    Raw_streams = zeros(90,90,78);
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

        a = readtable([path(subj_num).folder '/' path(subj_num).name '/Tractography_AgeDiff/NoCluster/Age_ttest/AAL90_5000stream/waytotal']);
        a = table2array(a);
        a = a(1:90,1);
        % Normalise each seed row by its waytotal before symmetrising.
        M = M./a;
        M = (M+M.')/2;
        Thresholding_Mat(:,:,ii) = M;
        disp(ii)
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

    tmp = load('/path/to/tempseq_dti/Graph_Theory_Age_Diff/ROI_sizes.mat')
    ROI_sizes = tmp.ROI_sizes;

    % Construct the stack of participant-specific ROI-size matrices.
    Roi_Sub_mat = zeros(90,90,78);
    for participant = 1:78
        Roi_subj = ROI_sizes(:,participant);
        % Form pairwise sums of ROI sizes for edge-weight normalisation.
        Roi_subj = Roi_subj+Roi_subj.';
        Roi_subj(~Edge_Mask) = 0; % Zero ROI-size denominators outside the selected edge mask.
        Roi_Sub_mat (:,:,participant) = Roi_subj;
    end

    Raw_streams = Raw_streams.* Edge_Mask;

    % This expression uses Roi_subj from the final loop iteration for every slice.
    Raw_streams_norm = Raw_streams./ Roi_Sub_mat;
    Raw_streams_norm(find(isnan(Raw_streams_norm))) = 0;

    No_Cluster_dir = '/path/to/tempseq_dti/Graph_Theory_Age_Diff/NoCluster';
    No_Cluster_dir = dir(No_Cluster_dir);

    outdir = ['/path/to/tempseq_dti/Graph_Theory_Age_Diff/permuts_Rich_club/density_' num2str(my_prctile)]
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end

    list = dir('/path/to/tempseq_dti/Graph_Theory_Age_Diff/NoCluster/Subj*'); % list of subjects
    nperm = S.nperm;
    perm_mat = zeros(nperm,90,length(list));

    % Compute a rich-club coefficient curve for each random network.
    for ii = 1:length(list) % over subjects
        % Load the saved subject variables before selecting the normalised connectivity slice.
        load([list(ii).folder '/' list(ii).name]); % Load variables from the subject MAT file.
        % M3(~Edge_Mask) = 0;
        % M3=Lesion_Mat;
        M3 = Raw_streams_norm(:,:,ii);
        nRand=nperm;
        phi_rand = zeros(nRand,90);
        for i = 1:nRand
            % Generate the randomised network using the BCT routine.
            rand_net = randmio_und_connected(M3,25);
            this_dist = rich_club_wu(rand_net);
            phi_rand(i,1:length(this_dist)) = this_dist;

        end
        perm_mat(:,:,ii) = phi_rand;
    end

    this_iteration = S.iteration;

    filename = fullfile(outdir,['permut_richclub' num2str(this_iteration) '.mat']);
    % Save this iteration of rich-club permutation curves.
    save(filename,'perm_mat');

    output = 1;
% Report a caught error and return the failure status.
catch ME
    fprintf('ERROR in Perm_function_ANOVA_continious: %s/m',ME.message);
    output = -1;
end

end
