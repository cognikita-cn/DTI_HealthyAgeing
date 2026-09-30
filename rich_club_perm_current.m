function output = rich_club_perm_current(S)
% RICH_CLUB_PERM_CURRENT Generate rich-club null networks for the current pipeline.
%
% Required fields:
%   S.cache_file
%   S.edge_mask_file
%   S.outdir
%   S.bct_dir
%   S.perc
%   S.nperm
%   S.iteration
%   S.seed
%
% Optional:
%   S.partStart
%
% The null model itself is unchanged from the historical analysis:
% null_model_und_sign(M3,5,0.1). The changes are reproducible seeding,
% current cached matrices, and the current saved edge masks.

try

    %% BCT

    bct_paths = strsplit(genpath(S.bct_dir),pathsep);

    for jj = 1:length(bct_paths)

        if isempty(bct_paths{jj})
            continue
        end

        if exist(fullfile(bct_paths{jj},'rich_club_wu.m'),'file') || ...
           exist(fullfile(bct_paths{jj},'null_model_und_sign.m'),'file')

            addpath(bct_paths{jj});

        end
    end

    if isempty(which('rich_club_wu'))
        error('BCT function rich_club_wu not found.');
    end

    if isempty(which('null_model_und_sign'))
        error('BCT function null_model_und_sign not found.');
    end


    %% Settings

    density = S.perc;
    nperm = S.nperm;
    iteration = S.iteration;

    if isfield(S,'partStart')
        partStart = S.partStart;
    else
        partStart = 1;
    end


    %% Load current cache

    tmp = load(S.cache_file,'cache');
    cache = tmp.cache;

    Raw_streams = cache.Raw_streams;
    ROI_sizes = cache.ROI_sizes;
    subjects = cache.subjects;

    npart = size(Raw_streams,3);

    if size(ROI_sizes,2) ~= npart
        error('ROI sizes do not match participant count.');
    end

    if length(subjects) ~= npart
        error('Subject list does not match participant count.');
    end


    %% Current saved edge mask

    mask_data = load(S.edge_mask_file,'edge_masks','dens_range');

    density_ind = find(mask_data.dens_range == density,1);

    if isempty(density_ind)
        error('Density %g not found in edge-mask file.',density);
    end

    Edge_Mask = logical(mask_data.edge_masks(:,:,density_ind));

    if ~isequal(size(Edge_Mask),[90 90]) || ...
            ~isequal(Edge_Mask,Edge_Mask.')

        error('Edge mask must be a symmetric 90 x 90 matrix.');
    end


    %% Output

    outdir = fullfile( ...
        S.outdir, ...
        'RichClub_Permuts', ...
        ['density_' num2str(density)]);

    if ~exist(outdir,'dir')
        mkdir(outdir);
    end

    perm_mat = nan(nperm,90,npart);


    %% Participant-level null networks

    for participant = partStart:npart

        participant_seed = ...
            S.seed + ...
            1000000*iteration + ...
            1000*density + ...
            participant;

        stream = RandStream('mt19937ar','Seed',participant_seed);
        RandStream.setGlobalStream(stream);

        Roi_subj = ROI_sizes(:,participant);
        Roi_pair = Roi_subj + Roi_subj.';

        M3 = ...
            (Raw_streams(:,:,participant) .* Edge_Mask) ./ ...
            Roi_pair;

        M3(~isfinite(M3)) = 0;
        M3(1:91:end) = 0;
        M3 = (M3 + M3.')/2;

        for perm = 1:nperm

            rand_net = null_model_und_sign(M3,5,0.1);
            this_rc = rich_club_wu(rand_net);
            this_rc = this_rc(:).';

            perm_mat(perm,1:length(this_rc),participant) = ...
                this_rc;

        end

        fprintf( ...
            '%d/%d | %s | density %g | iteration %d\n', ...
            participant,npart,subjects{participant},density,iteration);

    end


    %% Save

    filename = fullfile( ...
        outdir, ...
        ['permut_richclub_iter_' num2str(iteration) '.mat']);

    save(filename,'perm_mat','-v7.3');

    output = 1;

catch ME

    fprintf(2,'ERROR in rich_club_perm_current: %s\n',ME.message);

    for ii = 1:length(ME.stack)
        fprintf(2,'  %s, line %d\n', ...
            ME.stack(ii).name, ...
            ME.stack(ii).line);
    end

    output = -1;

end

end
