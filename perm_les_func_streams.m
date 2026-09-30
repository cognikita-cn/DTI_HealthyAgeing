function output = perm_les_func_streams(S)
% PERM_LES_FUNC_STREAMS Generate mass-matched streamline permutations.
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

try

    %% BCT

    bct_paths = strsplit(genpath(S.bct_dir),pathsep);

    for jj = 1:length(bct_paths)

        if isempty(bct_paths{jj})
            continue
        end

        if exist(fullfile(bct_paths{jj},'efficiency_wei.m'),'file') || ...
           exist(fullfile(bct_paths{jj},'clustering_coef_wu.m'),'file')
            addpath(bct_paths{jj});
        end
    end

    if isempty(which('efficiency_wei'))
        error('BCT function efficiency_wei not found.');
    end

    if isempty(which('clustering_coef_wu'))
        error('BCT function clustering_coef_wu not found.');
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

    %% Load input data

    tmp = load(S.cache_file,'cache');
    cache = tmp.cache;

    Raw_streams = cache.Raw_streams;
    Cluster_streams = cache.Cluster_streams;
    ROI_sizes = cache.ROI_sizes;
    subjects = cache.subjects;
    healthy_selection = logical(cache.healthy_selection(:));

    npart = size(Raw_streams,3);

    if size(Cluster_streams,3) ~= npart
        error('Raw and cluster matrices have different participant counts.');
    end

    if size(ROI_sizes,2) ~= npart
        error('ROI sizes do not match connectivity participant count.');
    end

    if length(subjects) ~= npart
        error('Subject list does not match connectivity participant count.');
    end

    if length(healthy_selection) ~= npart
        error('Healthy selection does not match participant count.');
    end

    participant_ind = find(healthy_selection);
    participant_ind = participant_ind(participant_ind >= partStart);

    if isempty(participant_ind)
        error('No healthy participants selected.');
    end

    %% Edge mask

    mask_data = load( ...
        S.edge_mask_file, ...
        'edge_masks', ...
        'dens_range');

    density_ind = find(mask_data.dens_range == density,1);

    if isempty(density_ind)
        error('Density %g not found in edge-mask file.',density);
    end

    Edge_Mask = logical(mask_data.edge_masks(:,:,density_ind));

    if ~isequal(size(Edge_Mask),[90 90])
        error('Edge mask must be 90 x 90.');
    end

    if ~isequal(Edge_Mask,Edge_Mask.')
        error('Edge mask is not symmetric.');
    end

    Edge_Mask_tri = triu(Edge_Mask,1);
    edge_ind = find(Edge_Mask_tri);

    %% Lesion mass

    clust_sums = zeros(npart,1);

    for ii = participant_ind'

        cluster_units = round(2 * Cluster_streams(:,:,ii));
        clust_sums(ii) = sum(cluster_units(edge_ind));

    end

    %% Output directory

    outdir = fullfile( ...
        S.outdir, ...
        'Streamline_Permuts', ...
        ['density_' num2str(density)]);

    if ~exist(outdir,'dir')
        mkdir(outdir);
    end

    %% Permutations

    for hh = 1:length(participant_ind)

        participant = participant_ind(hh);

        participant_seed = ...
            S.seed + ...
            1000000 * iteration + ...
            1000 * density + ...
            participant;

        stream = RandStream('mt19937ar','Seed',participant_seed);
        RandStream.setGlobalStream(stream);

        efficiency = zeros(nperm,1);
        clustering_coef = zeros(90,nperm);
        degree_les = zeros(90,nperm);

        Roi_subj = ROI_sizes(:,participant);
        Roi_pair = Roi_subj + Roi_subj.';

        this_raw_mat = Raw_streams(:,:,participant);

        raw_units = round(2 * this_raw_mat);
        raw_units(~Edge_Mask_tri) = 0;

        raw_vals = raw_units(edge_ind);

        total_streams = sum(raw_vals);
        this_lesion = clust_sums(participant);

        if total_streams == 0
            error( ...
                'Participant %d has no eligible connectivity.', ...
                participant);
        end

        if this_lesion > total_streams
            error( ...
                ['Participant %d: lesion mass exceeds available ' ...
                 'connectivity.'], ...
                participant);
        end

        Multinom_prob = raw_vals / total_streams;
        Multinom_prob = Multinom_prob / sum(Multinom_prob);
        Multinom_prob = Multinom_prob';

        %% Bookkeeping

        bookkeeping = struct();

        bookkeeping.participant_index = participant;
        bookkeeping.subject = subjects{participant};
        bookkeeping.density = density;
        bookkeeping.iteration = iteration;
        bookkeeping.seed = participant_seed;

        bookkeeping.nperm = nperm;
        bookkeeping.n_eligible_edges = length(edge_ind);

        bookkeeping.total_raw_units = total_streams;
        bookkeeping.observed_lesion_units = this_lesion;
        bookkeeping.lesion_fraction = ...
            this_lesion / total_streams;

        this_cluster_units = ...
            round(2 * Cluster_streams(:,:,participant));

        bookkeeping.observed_affected_edges = ...
            nnz(this_cluster_units(edge_ind) > 0);

        bookkeeping.attempted_draws = 0;
        bookkeeping.rejected_draws = 0;

        %% Generate permutations

        perm = 1;

        while perm <= nperm

            bookkeeping.attempted_draws = ...
                bookkeeping.attempted_draws + 1;

            removed_units = ...
                mnrnd(this_lesion,Multinom_prob);

            residual_units = ...
                raw_vals - removed_units';

            if any(residual_units < 0)

                bookkeeping.rejected_draws = ...
                    bookkeeping.rejected_draws + 1;

                continue

            end

            residual_weights = residual_units / 2;

            Lesion_perm = zeros(90,90);
            Lesion_perm(edge_ind) = residual_weights;

            Lesion_perm(edge_ind) = ...
                Lesion_perm(edge_ind) ./ ...
                Roi_pair(edge_ind);

            Lesion_perm = ...
                Lesion_perm + Lesion_perm.';

            efficiency(perm) = ...
                efficiency_wei(Lesion_perm);

            clustering_coef(:,perm) = ...
                clustering_coef_wu(Lesion_perm);

            degree_les(:,perm) = ...
                sum(Lesion_perm,2);

            perm = perm + 1;

        end

        %% Finalise bookkeeping

        bookkeeping.accepted_draws = nperm;

        bookkeeping.rejection_rate = ...
            bookkeeping.rejected_draws / ...
            bookkeeping.attempted_draws;

        %% Save

        this_part_struct = struct();

        this_part_struct.efficiency = efficiency;
        this_part_struct.clustering_coef = clustering_coef;
        this_part_struct.degree_les = degree_les;
        this_part_struct.bookkeeping = bookkeeping;

        filename = fullfile( ...
            outdir, ...
            ['permut_les_mat' num2str(participant) ...
             '_iter_' num2str(iteration) '.mat']);

        save(filename,'this_part_struct');

        fprintf( ...
            ['%d/%d | participant %d | %s | density %g | ' ...
             'iteration %d | rejected %d/%d\n'], ...
            hh, ...
            length(participant_ind), ...
            participant, ...
            subjects{participant}, ...
            density, ...
            iteration, ...
            bookkeeping.rejected_draws, ...
            bookkeeping.attempted_draws);

    end

    output = 1;

catch ME

    fprintf(2,'ERROR in perm_les_func_streams: %s\n',ME.message);

    for ii = 1:length(ME.stack)
        fprintf(2,'  %s, line %d\n', ...
            ME.stack(ii).name, ...
            ME.stack(ii).line);
    end

    output = -1;

end

end