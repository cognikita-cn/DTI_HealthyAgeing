function [distance_mat, centroids] = get_aal_centroid_distances(masks_file)
% GET_AAL_CENTROID_DISTANCES Calculate Euclidean distances between AAL ROIs.
%
% Inputs:
%   masks_file    Path to masks.txt containing one NIfTI mask per line.
%
% Outputs:
%   distance_mat  90 x 90 Euclidean distance matrix in mm.
%   centroids     90 x 3 ROI centroid coordinates in mm.

fid = fopen(masks_file,'r');

if fid == -1
    error('Could not open %s.',masks_file);
end

C = textscan(fid,'%s','Delimiter','\n');
fclose(fid);

mask_files = strtrim(C{1});
mask_files = mask_files(~cellfun('isempty',mask_files));

nroi = length(mask_files);

if nroi ~= 90
    error('Expected 90 ROI masks, found %d in %s.',nroi,masks_file);
end

centroids = zeros(nroi,3);

tmpdir = tempname;
mkdir(tmpdir);

cleanupObj = onCleanup(@() rmdir(tmpdir,'s'));

for roi = 1:nroi

    mask_file = mask_files{roi};

    if ~exist(mask_file,'file')
        error('Mask not found: %s',mask_file);
    end

    %% Uncompress NIfTI if required

    if length(mask_file) > 7 && strcmp(mask_file(end-6:end),'.nii.gz')

        [status,~] = system(sprintf( ...
            'gunzip -c "%s" > "%s"', ...
            mask_file, ...
            fullfile(tmpdir,sprintf('roi_%03d.nii',roi))));

        if status ~= 0
            error('Could not decompress %s.',mask_file);
        end

        nii_file = fullfile(tmpdir,sprintf('roi_%03d.nii',roi));

    else

        nii_file = mask_file;

    end

    %% Load mask

    V = spm_vol(nii_file);
    Y = spm_read_vols(V);

    ind = find(Y > 0);

    if isempty(ind)
        error('ROI %d is empty: %s',roi,mask_file);
    end

    %% Voxel coordinates

    [i,j,k] = ind2sub(size(Y),ind);

    xyz_vox = [ ...
        double(i)'; ...
        double(j)'; ...
        double(k)'; ...
        ones(1,length(i))];

    %% Native-space coordinates in mm

    xyz_mm = V.mat * xyz_vox;

    centroids(roi,:) = mean(xyz_mm(1:3,:),2)';

end

%% Pairwise Euclidean distances

distance_mat = zeros(nroi,nroi);

for roi = 1:nroi

    delta = centroids - centroids(roi,:);

    distance_mat(roi,:) = ...
        sqrt(sum(delta.^2,2));

end

distance_mat(1:nroi+1:end) = 0;

end