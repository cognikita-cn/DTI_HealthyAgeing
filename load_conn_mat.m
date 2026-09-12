function output = load_conn_mat(~)
% LOAD_CONN_MAT Load TempSeq raw and cluster-constrained connectivity.
%   output = load_conn_mat(~)
%
%   The input argument is ignored. Subject, tractography and ROI-size paths
%   are specified below. Connectivity arrays are preallocated for 78 subjects.
%   Subject numbers parsed from Subj* filenames index the subject directory
%   listing; output slices follow the Cluster_dir listing order.
%
%   Returned fields:
%     Raw_streams       Symmetrised raw weights in the specified AAL order.
%     Cluster_streams   Symmetrised cluster-constrained weights in that order.
%     Lesion_mats       Raw_streams minus Cluster_streams.
%     Prop_Mats         Cluster_streams divided elementwise by Raw_streams.
%     Thresholding_Mat  Row-wise waytotal-normalised weights, then symmetrised;
%                       this field remains in the initial AAL order.
%     ROI_sizes         ROI-size array loaded from disk without reordering.
%
%   Connectivity fields have dimensions 90-by-90-by-participant. Division
%   results are returned directly; this function does not replace NaN/Inf.
%   Replace /path/to roots before execution.

output = [];
Cluster_dir = '/path/to/tempseq_dti/Graph_Theory_Age_Diff/ThroughClust/Subj*';
Cluster_dir = dir(Cluster_dir);

path = dir('/path/to/tempseq_dti/00*'); % path to subjects
% Allocate ROI-by-ROI-by-participant connectivity arrays.
Thresholding_Mat = zeros(90,90,78);
Raw_streams = zeros(90,90,78);
Lesion_mats = zeros(90,90,78);
Cluster_streams = zeros(90,90,78);
Prop_Mats = zeros(90,90,78);
order = [1:2:90 90:-2:2];
for ii = 1:length(Cluster_dir)
    subj = ii;
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
    % Subtract cluster-constrained weights to obtain the virtual-lesion matrix.
    Lesion_mats(:,:,ii) = M_raw - M_Clust;
    Prop_Mats(:,:,ii) = M_Clust./M_raw;
    a = readtable([path(subj_num).folder '/' path(subj_num).name '/Tractography_AgeDiff/NoCluster/Age_ttest/AAL90_5000stream/waytotal']);
    a = table2array(a);
    a = a(1:90,1);
    % Normalise each seed row by its waytotal before symmetrising.
    M = M./a;
    M = (M+M.')/2;
    Thresholding_Mat(:,:,ii) = M;
end

tmp = load('/path/to/tempseq_dti/Graph_Theory_Age_Diff/ROI_sizes.mat');
ROI_sizes = tmp.ROI_sizes;
% Return connectivity fields and the ROI-size array.
output.Thresholding_Mat = Thresholding_Mat;
output.Raw_streams = Raw_streams;
output.Cluster_streams = Cluster_streams;
output.Lesion_mats = Lesion_mats;
output.Prop_Mats = Prop_Mats;
output.ROI_sizes = ROI_sizes;
disp('Matrices Loaded Successfully')
end
