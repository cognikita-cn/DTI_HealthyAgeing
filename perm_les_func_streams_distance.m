function output = perm_les_func_streams_distance(S)
% PERM_LES_FUNC_STREAMS_DISTANCE
% Generate distance-matched streamline lesion permutations.
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
    nbins = 5;

    if isfield(S,'partStart')
        partStart = S.partStart;
    else
        partStart = 1;
    end

    %% Load cache

    tmp = load(S.cache_file,'cache');
    cache = tmp.cache;

    if ~isfield(cache,'Centroid_distance')
        error('Centroid_distance not found in permutation cache.');
    end

    Raw_streams = cache.Raw_streams;
    Cluster_streams = cache.Cluster_streams;
    ROI_sizes = cache.ROI_sizes;
    Centroid_distance = cache.Centroid_distance;
    subjects = cache.subjects;
    healthy_selection = logical(cache.healthy_selection(:));

    npart = size(Raw_streams,3);

    if size(Cluster_streams,3) ~= npart
        error('Raw and cluster matrices have different participant counts.');
    end

    if size(Centroid_distance,3) ~= npart
        error('Distance matrices do not match participant count.');
    end

    if size(ROI_sizes,2) ~= npart
        error('ROI sizes do not match participant count.');
    end

    if length(subjects) ~= npart
        error('Subject list does not match participant count.');
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

    mask_data = load(S.edge_mask_file,'edge_masks','dens_range');

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

    %% Output directory

    outdir = fullfile( ...
        S.outdir, ...
        'Streamline_Permuts', ...
        ['density_' num2str(density)]);

    if ~exist(outdir,'dir')
        mkdir(outdir);
    end

    %% Participants

    for hh = 1:length(participant_ind)

        participant = participant_ind(hh);

        participant_seed = ...
            S.seed + ...
            1000000 * iteration + ...
            1000 * density + ...
            participant;

        stream = RandStream( ...
            'mt19937ar', ...
            'Seed',participant_seed);

        RandStream.setGlobalStream(stream);

        efficiency = zeros(nperm,1);
        clustering_coef = zeros(90,nperm);
        degree_les = zeros(90,nperm);

        affected_edges = zeros(nperm,1);
        affected_edges_bin = zeros(nbins,nperm);

        %% Participant data

        Roi_subj = ROI_sizes(:,participant);
        Roi_pair = Roi_subj + Roi_subj.';

        raw_units = round(2 * Raw_streams(:,:,participant));
        raw_units(~Edge_Mask_tri) = 0;
        raw_vals = raw_units(edge_ind);

        cluster_units = ...
            round(2 * Cluster_streams(:,:,participant));

        cluster_vals = cluster_units(edge_ind);

        distance_mat = Centroid_distance(:,:,participant);
        distance_vals = distance_mat(edge_ind);

        if any(~isfinite(distance_vals))
            error('Non-finite distances for participant %d.',participant);
        end

        %% Distance quintiles

        distance_percentiles = ...
            prctile(distance_vals,[20 40 60 80]);

        distance_bin = ones(length(edge_ind),1);

        distance_bin(distance_vals > distance_percentiles(1)) = 2;
        distance_bin(distance_vals > distance_percentiles(2)) = 3;
        distance_bin(distance_vals > distance_percentiles(3)) = 4;
        distance_bin(distance_vals > distance_percentiles(4)) = 5;

        n_edges_per_bin = zeros(nbins,1);
        lesion_units_per_bin = zeros(nbins,1);
        raw_units_per_bin = zeros(nbins,1);
        observed_affected_edges_per_bin = zeros(nbins,1);

        Multinom_prob = cell(nbins,1);

        for bb = 1:nbins

            bin_ind = distance_bin == bb;

            n_edges_per_bin(bb) = sum(bin_ind);

            if n_edges_per_bin(bb) == 0
                error('Empty distance bin %d for participant %d.', ...
                    bb,participant);
            end

            raw_bin = raw_vals(bin_ind);
            cluster_bin = cluster_vals(bin_ind);

            raw_units_per_bin(bb) = sum(raw_bin);
            lesion_units_per_bin(bb) = sum(cluster_bin);

            observed_affected_edges_per_bin(bb) = ...
                nnz(cluster_bin > 0);

            if lesion_units_per_bin(bb) > raw_units_per_bin(bb)

                error(['Participant %d, bin %d: lesion mass exceeds ' ...
                       'available connectivity.'], ...
                       participant,bb);

            end

            if lesion_units_per_bin(bb) > 0

                if raw_units_per_bin(bb) == 0
                    error('No eligible connectivity in bin %d.',bb);
                end

                p = raw_bin / raw_units_per_bin(bb);
                Multinom_prob{bb} = p';

            end

        end

        %% Bookkeeping

        bookkeeping = struct();

        bookkeeping.participant_index = participant;
        bookkeeping.subject = subjects{participant};
        bookkeeping.density = density;
        bookkeeping.iteration = iteration;
        bookkeeping.seed = participant_seed;
        bookkeeping.nperm = nperm;

        bookkeeping.n_distance_bins = nbins;
        bookkeeping.distance_percentiles = distance_percentiles;
        bookkeeping.n_edges_per_bin = n_edges_per_bin;

        bookkeeping.total_raw_units = sum(raw_vals);
        bookkeeping.observed_lesion_units = sum(cluster_vals);

        bookkeeping.raw_units_per_bin = raw_units_per_bin;
        bookkeeping.observed_lesion_units_per_bin = ...
            lesion_units_per_bin;

        bookkeeping.observed_affected_edges = ...
            nnz(cluster_vals > 0);

        bookkeeping.observed_affected_edges_per_bin = ...
            observed_affected_edges_per_bin;

        bookkeeping.attempted_draws = 0;
        bookkeeping.rejected_draws = 0;

        %% Permutations

        perm = 1;

        while perm <= nperm

            bookkeeping.attempted_draws = ...
                bookkeeping.attempted_draws + 1;

            removed_units = zeros(length(edge_ind),1);
            valid_draw = true;

            for bb = 1:nbins

                bin_ind = distance_bin == bb;

                if lesion_units_per_bin(bb) == 0
                    continue
                end

                removed_bin = mnrnd( ...
                    lesion_units_per_bin(bb), ...
                    Multinom_prob{bb});

                removed_units(bin_ind) = removed_bin';

                if any(raw_vals(bin_ind) - removed_bin' < 0)

                    valid_draw = false;
                    break

                end

            end

            if ~valid_draw

                bookkeeping.rejected_draws = ...
                    bookkeeping.rejected_draws + 1;

                continue

            end

            residual_units = raw_vals - removed_units;

            %% Reconstruct network

            residual_weights = residual_units / 2;

            Lesion_perm = zeros(90,90);

            Lesion_perm(edge_ind) = residual_weights;

            Lesion_perm(edge_ind) = ...
                Lesion_perm(edge_ind) ./ Roi_pair(edge_ind);

            Lesion_perm = Lesion_perm + Lesion_perm.';

            %% Graph measures

            efficiency(perm) = ...
                efficiency_wei(Lesion_perm);

            clustering_coef(:,perm) = ...
                clustering_coef_wu(Lesion_perm);

            degree_les(:,perm) = ...
                sum(Lesion_perm,2);

            %% Affected edges

            affected_edges(perm) = ...
                nnz(removed_units > 0);

            for bb = 1:nbins

                bin_ind = distance_bin == bb;

                affected_edges_bin(bb,perm) = ...
                    nnz(removed_units(bin_ind) > 0);

            end

            perm = perm + 1;

        end

        %% Final bookkeeping

        bookkeeping.accepted_draws = nperm;

        bookkeeping.rejection_rate = ...
            bookkeeping.rejected_draws / ...
            bookkeeping.attempted_draws;

        %% Save

        this_part_struct = struct();

        this_part_struct.efficiency = efficiency;
        this_part_struct.clustering_coef = clustering_coef;
        this_part_struct.degree_les = degree_les;

        this_part_struct.affected_edges = affected_edges;
        this_part_struct.affected_edges_bin = affected_edges_bin;

        this_part_struct.bookkeeping = bookkeeping;

        filename = fullfile( ...
            outdir, ...
            ['permut_les_dist_mat' ...
             num2str(participant) ...
             '_iter_' ...
             num2str(iteration) ...
             '.mat']);

        save(filename,'this_part_struct');

        fprintf([ ...
            '%d/%d | participant %d | %s | density %g | ' ...
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

    fprintf(2, ...
        'ERROR in perm_les_func_streams_distance: %s\n', ...
        ME.message);

    for ii = 1:length(ME.stack)

        fprintf(2, ...
            '  %s, line %d\n', ...
            ME.stack(ii).name, ...
            ME.stack(ii).line);

    end

    output = -1;

end

end