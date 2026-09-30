function output = load_conn_mat_dem(S)
% LOAD_CONN_MAT_DEM Load AAL90 connectivity matrices in valid-DTI order.
%
% The function preserves the original analysis definitions:
% - valid DTI participants are subjects with a 20* acquisition entry;
% - raw and cluster matrices are symmetrised and reordered to the AAL90
%   left-right ordering used throughout the graph-theory pipeline;
% - lesion matrices are raw minus cluster-through streamlines;
% - thresholding matrices are raw fdt_network_matrix values divided by
%   waytotal before symmetrisation. They remain in native AAL file order
%   because the combined edge-mask code applies the AAL reordering once.

output = [];

path = dir(S.path);
path = path([path.isdir]);

has_DTI = false(length(path),1);

for gg = 1:length(path)

    path6 = dir(fullfile( ...
        path(gg).folder, ...
        path(gg).name, ...
        '20*'));

    has_DTI(gg) = ~isempty(path6);

end

path = path(has_DTI);
DTI_Count = length(path);

Thresholding_Mat = zeros(90,90,DTI_Count);
Raw_streams = zeros(90,90,DTI_Count);
Lesion_mats = zeros(90,90,DTI_Count);
Cluster_streams = zeros(90,90,DTI_Count);
Prop_Mats = zeros(90,90,DTI_Count);

subjects = cell(DTI_Count,1);

order = [1:2:90 90:-2:2];

subdir_raw = regexprep(S.subdir_raw,'^[\\/]+','');
subdir_clust = regexprep(S.subdir_clust,'^[\\/]+','');

for ii = 1:DTI_Count

    subject_dir = fullfile(path(ii).folder,path(ii).name);

    raw_file = fullfile(subject_dir,subdir_raw,'fdt_network_matrix');
    clust_file = fullfile(subject_dir,subdir_clust,'fdt_network_matrix');
    waytotal_file = fullfile(subject_dir,subdir_raw,'waytotal');

    if ~exist(raw_file,'file')
        error('Missing raw connectivity matrix for %s.',path(ii).name);
    end

    if ~exist(clust_file,'file')
        error('Missing cluster connectivity matrix for %s.',path(ii).name);
    end

    if ~exist(waytotal_file,'file')
        error('Missing waytotal for %s.',path(ii).name);
    end

    M = readtable(raw_file);
    M = table2array(M);
    M = M(1:90,1:90);

    M_raw = (M + M.')/2;
    M_raw = M_raw(order,order);

    Raw_streams(:,:,ii) = M_raw;

    M_Clust = readtable(clust_file);
    M_Clust = table2array(M_Clust);
    M_Clust = M_Clust(1:90,1:90);

    M_Clust = (M_Clust + M_Clust.')/2;
    M_Clust = M_Clust(order,order);

    Cluster_streams(:,:,ii) = M_Clust;
    Lesion_mats(:,:,ii) = M_raw - M_Clust;
    Prop_Mats(:,:,ii) = M_Clust ./ M_raw;

    waytotal = readtable(waytotal_file);
    waytotal = table2array(waytotal);
    waytotal = waytotal(1:90,1);

    M_thresh = M ./ waytotal;
    M_thresh = (M_thresh + M_thresh.')/2;

    Thresholding_Mat(:,:,ii) = M_thresh;
    subjects{ii} = path(ii).name;

end

tmp = load(S.Roi_sizes,'ROI_sizes');
ROI_sizes = tmp.ROI_sizes;

if size(ROI_sizes,2) ~= DTI_Count
    error(['ROI_sizes contains %d participants, but %d valid-DTI ' ...
           'participants were found.'], ...
           size(ROI_sizes,2),DTI_Count);
end

output.Thresholding_Mat = Thresholding_Mat;
output.Raw_streams = Raw_streams;
output.Cluster_streams = Cluster_streams;
output.Lesion_mats = Lesion_mats;
output.Prop_Mats = Prop_Mats;
output.ROI_sizes = ROI_sizes;
output.subjects = subjects;

disp('Matrices Loaded Successfully')

end
