function cache = make_permutation_cache(S,healthy_selection)

order = [1:2:90 90:-2:2];

%% Participants with DTI data

path_all = dir(S.path);
path_all = path_all([path_all.isdir]);

has_DTI = false(length(path_all),1);

for ii = 1:length(path_all)

    dti_data = dir(fullfile( ...
        path_all(ii).folder, ...
        path_all(ii).name, ...
        '20*'));

    has_DTI(ii) = ~isempty(dti_data);

end

path = path_all(has_DTI);
npart = length(path);

fprintf('Subject folders: %d\n',length(path_all));
fprintf('Subjects with DTI: %d\n',npart);

%% ROI sizes

tmp = load(S.Roi_sizes,'ROI_sizes');
ROI_sizes = tmp.ROI_sizes;

if size(ROI_sizes,2) ~= npart
    error(['ROI_sizes contains %d participants, but %d participants ' ...
           'with DTI data were found.'], ...
           size(ROI_sizes,2),npart);
end

healthy_selection = logical(healthy_selection(:));

if length(healthy_selection) ~= npart
    error(['Healthy selection contains %d participants, but %d ' ...
           'participants with DTI data were found.'], ...
           length(healthy_selection),npart);
end

%% Allocate

Raw_streams = zeros(90,90,npart);
Cluster_streams = zeros(90,90,npart);

Centroid_distance = zeros(90,90,npart);
Centroids = zeros(90,3,npart);

subjects = cell(npart,1);

%% Load participant data

for ii = 1:npart

    subject_dir = fullfile( ...
        path(ii).folder, ...
        path(ii).name);

    raw_file = fullfile( ...
        subject_dir, ...
        S.subdir_raw, ...
        'fdt_network_matrix');

    clust_file = fullfile( ...
        subject_dir, ...
        S.subdir_clust, ...
        'fdt_network_matrix');

    masks_file = fullfile( ...
        subject_dir, ...
        'masks.txt');

    if ~exist(raw_file,'file')
        error('Missing raw connectivity matrix for %s.',path(ii).name);
    end

    if ~exist(clust_file,'file')
        error('Missing cluster connectivity matrix for %s.',path(ii).name);
    end

    if ~exist(masks_file,'file')
        error('Missing masks.txt for %s.',path(ii).name);
    end

    %% Raw connectivity

    M = readtable(raw_file);
    M = table2array(M);
    M = M(1:90,1:90);

    M = (M + M.') / 2;
    M = M(order,order);

    Raw_streams(:,:,ii) = M;

    %% Cluster connectivity

    M = readtable(clust_file);
    M = table2array(M);
    M = M(1:90,1:90);

    M = (M + M.') / 2;
    M = M(order,order);

    Cluster_streams(:,:,ii) = M;

    %% AAL centroid distances

    [distance_mat,centroids] = ...
        get_aal_centroid_distances(masks_file);

    distance_mat = distance_mat(order,order);
    centroids = centroids(order,:);

    Centroid_distance(:,:,ii) = distance_mat;
    Centroids(:,:,ii) = centroids;

    subjects{ii} = path(ii).name;

    fprintf('%d/%d | %s\n',ii,npart,path(ii).name);

end

%% Cache

cache = struct();

cache.Raw_streams = Raw_streams;
cache.Cluster_streams = Cluster_streams;

cache.ROI_sizes = ROI_sizes;

cache.Centroid_distance = Centroid_distance;
cache.Centroids = Centroids;

cache.subjects = subjects;
cache.healthy_selection = healthy_selection;

cache.has_DTI = has_DTI;
cache.all_subjects = {path_all.name}';

end