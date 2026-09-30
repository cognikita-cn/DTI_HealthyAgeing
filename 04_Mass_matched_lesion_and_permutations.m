%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% POST-REVIEW PERMUTATION NULL MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run sections individually. The active manuscript analyses use the
% healthy-only caches, the saved combined edge masks, reproducible seeds,
% and the current mass- and distance-matched permutation functions.

%% ------------------------------------------------------------------------
% POST-REVIEW MASS-MATCHED PERMUTATION QC / TEST UTILITIES
% -------------------------------------------------------------------------
% These blocks were previously stored in the combined-analysis setup script.
% Run the relevant section independently when testing the permutation environment.

%%
addpath('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/Functions')
root_dir = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';
S = struct();

S.cache_file = ...
    fullfile(root_dir,'TempSeq','permutation_input.mat');

S.edge_mask_file = ...
    fullfile(root_dir,'combined_edge_masks.mat');

S.outdir = ...
    fullfile(root_dir,'TempSeq');

S.perc = 10;
S.nperm = 5;
S.iteration = 99;
S.seed = 100000000;
S.partStart = 1;

clusterconfig('scheduler','cluster');
clusterconfig('long_running',1);
clusterconfig('slot',2);

test_job = job2cluster(@perm_les_func_streams,S);

%%
which perm_les_func_streams -all
dbtype perm_les_func_streams 118:130
which mnrnd -all
which efficiency_wei -all
which clustering_coef_wu -all
which rng -all
%%
seed = 12345;

stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);
a = rand(1,5);

stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);
b = rand(1,5);

isequal(a,b)

%%
bct_dir = ...
    ['/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/' ...
     'Graph_Theory/BCT/2019_03_03_BCT'];
addpath(genpath(bct_dir))

which efficiency_wei -all
which clustering_coef_wu -all

%%
%% Test streamline lesion permutations

clear
clc

%% Paths

root_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';

function_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/Functions';

bct_dir = ...
    ['/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/' ...
     'Graph_Theory/BCT/2019_03_03_BCT'];

addpath(function_dir)

edge_mask_file = fullfile(root_dir,'combined_edge_masks.mat');

%% Test settings

S = struct();

S.cache_file = ...
    fullfile(root_dir,'TempSeq','permutation_input.mat');

S.edge_mask_file = edge_mask_file;

S.outdir = ...
    fullfile(root_dir,'TestRun','TempSeq');

S.bct_dir = bct_dir;

S.perc = 10;
S.nperm = 5;
S.iteration = 1;
S.seed = 100000000;
S.partStart = 1;

%% Validate input

if ~exist(S.cache_file,'file')
    error('Permutation input not found.');
end

if ~exist(S.edge_mask_file,'file')
    error('Combined edge-mask file not found.');
end

if ~exist(S.bct_dir,'dir')
    error('BCT directory not found.');
end

if ~exist(S.outdir,'dir')
    mkdir(S.outdir);
end

%% Cluster settings

clusterconfig('scheduler','cluster');
clusterconfig('long_running',1);
clusterconfig('slot',2);

%% Submit test

test_job = job2cluster(@perm_les_func_streams,S);

fprintf('Test job submitted: %s\n',num2str(test_job));

%% ------------------------------------------------------------------------
% POST-REVIEW MASS-MATCHED PERMUTATION PRODUCTION LAUNCHER
% -------------------------------------------------------------------------

%% Streamline lesion permutations

clear
clc

%% Paths

root_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';

function_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/Functions';

bct_dir = ...
    ['/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/' ...
     'Graph_Theory/BCT/2019_03_03_BCT'];

addpath(function_dir)

edge_mask_file = ...
    fullfile(root_dir,'combined_edge_masks.mat');

%% Settings

dens_range = [10 20 30 40 50];

nperm = 50;
niter = 20;

%% Datasets

D(1).name = 'TempSeq';

D(1).cache_file = ...
    fullfile(root_dir,'TempSeq','permutation_input.mat');

D(1).outdir = ...
    fullfile(root_dir,'TempSeq');

D(1).seed = 100000000;


D(2).name = 'AudMemDem';

D(2).cache_file = ...
    fullfile(root_dir,'AudMemDem','permutation_input.mat');

D(2).outdir = ...
    fullfile(root_dir,'AudMemDem');

D(2).seed = 200000000;

%% Validate input

if ~exist(edge_mask_file,'file')
    error('Combined edge-mask file not found.');
end

if ~exist(bct_dir,'dir')
    error('BCT directory not found.');
end

for dd = 1:length(D)

    if ~exist(D(dd).cache_file,'file')
        error('Permutation input not found for %s.',D(dd).name);
    end

    if ~exist(D(dd).outdir,'dir')
        mkdir(D(dd).outdir);
    end

    tmp = load(D(dd).cache_file,'cache');

    fprintf('%s: %d total DTI participants, %d healthy participants\n', ...
        D(dd).name, ...
        size(tmp.cache.Raw_streams,3), ...
        sum(tmp.cache.healthy_selection));

    clear tmp

end

%% Cluster settings

clusterconfig('scheduler','cluster');
clusterconfig('long_running',1);
clusterconfig('slot',2);

%% Submit jobs

jobids = cell(length(D),length(dens_range),niter);

for dd = 1:length(D)

    for pp = 1:length(dens_range)

        for iter = 1:niter

            S = struct();

            S.cache_file = D(dd).cache_file;
            S.edge_mask_file = edge_mask_file;
            S.outdir = D(dd).outdir;
            S.bct_dir = bct_dir;

            S.perc = dens_range(pp);
            S.nperm = nperm;
            S.iteration = iter;

            S.seed = D(dd).seed;
            S.partStart = 1;

            jobids{dd,pp,iter} = ...
                job2cluster(@perm_les_func_streams,S);

            fprintf('%s | density %d | iteration %d/%d | job %s\n', ...
                D(dd).name, ...
                S.perc, ...
                iter, ...
                niter, ...
                num2str(jobids{dd,pp,iter}));

        end

    end

end

%% Save job IDs

save( ...
    fullfile(root_dir,'permutation_jobids.mat'), ...
    'jobids', ...
    'D', ...
    'dens_range', ...
    'nperm', ...
    'niter', ...
    'edge_mask_file', ...
    'bct_dir');

fprintf('\nSubmitted %d jobs.\n',numel(jobids));

%% ------------------------------------------------------------------------
% POST-REVIEW DISTANCE-MATCHED PERMUTATION TEST
% -------------------------------------------------------------------------

%% Test distance-matched lesion permutations

clear
clc
root_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';

function_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/Functions';

bct_dir = ...
    ['/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/' ...
     'Graph_Theory/BCT/2019_03_03_BCT'];

addpath(function_dir)

S = struct();

S.cache_file = ...
    fullfile(root_dir,'TempSeq','permutation_input.mat');

S.edge_mask_file = ...
    fullfile(root_dir,'combined_edge_masks.mat');

S.outdir = ...
    fullfile(root_dir,'TempSeq_Dist');

S.bct_dir = bct_dir;

S.perc = 10;
S.nperm = 2;
S.iteration = 1;

S.seed = 300000000;
S.partStart = 1;

if ~exist(S.outdir,'dir')
    mkdir(S.outdir);
end

tmp = load(S.cache_file,'cache');

if ~isfield(tmp.cache,'Centroid_distance')
    error('Centroid_distance not found in cache.');
end

clear tmp

clusterconfig('scheduler','cluster');
clusterconfig('long_running',0);
clusterconfig('slot',2);

test_job = ...
    job2cluster(@perm_les_func_streams_distance,S);

fprintf('Distance test job: %s\n',num2str(test_job));

%% ------------------------------------------------------------------------
% POST-REVIEW DISTANCE-MATCHED PERMUTATION PRODUCTION LAUNCHER
% -------------------------------------------------------------------------

%% Distance-matched streamline lesion permutations

clear
clc

%% Paths

root_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';

function_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/Functions';

bct_dir = ...
    ['/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/' ...
     'Graph_Theory/BCT/2019_03_03_BCT'];

addpath(function_dir)

edge_mask_file = ...
    fullfile(root_dir,'combined_edge_masks.mat');

%% Settings

dens_range = [10 20 30 40 50];

nperm = 50;
niter = 20;

%% Datasets

D(1).name = 'TempSeq';

D(1).cache_file = ...
    fullfile(root_dir,'TempSeq','permutation_input.mat');

D(1).outdir = ...
    fullfile(root_dir,'TempSeq_Dist');

D(1).seed = 300000000;


D(2).name = 'AudMemDem';

D(2).cache_file = ...
    fullfile(root_dir,'AudMemDem','permutation_input.mat');

D(2).outdir = ...
    fullfile(root_dir,'AudMemDem_Dist');

D(2).seed = 400000000;

%% Validate

if ~exist(edge_mask_file,'file')
    error('Combined edge-mask file not found.');
end

if ~exist(bct_dir,'dir')
    error('BCT directory not found.');
end

for dd = 1:length(D)

    if ~exist(D(dd).cache_file,'file')
        error('Permutation input not found for %s.',D(dd).name);
    end

    if ~exist(D(dd).outdir,'dir')
        mkdir(D(dd).outdir);
    end

    tmp = load(D(dd).cache_file,'cache');

    if ~isfield(tmp.cache,'Centroid_distance')
        error('Centroid_distance missing for %s.',D(dd).name);
    end

    fprintf('%s: %d total DTI participants, %d healthy participants\n', ...
        D(dd).name, ...
        size(tmp.cache.Raw_streams,3), ...
        sum(tmp.cache.healthy_selection));

    clear tmp

end

%% Cluster

clusterconfig('scheduler','cluster');
clusterconfig('long_running',1);
clusterconfig('slot',2);

%% Submit

jobids = cell(length(D),length(dens_range),niter);

for dd = 1:length(D)

    for pp = 1:length(dens_range)

        for iter = 1:niter

            S = struct();

            S.cache_file = D(dd).cache_file;
            S.edge_mask_file = edge_mask_file;
            S.outdir = D(dd).outdir;
            S.bct_dir = bct_dir;

            S.perc = dens_range(pp);
            S.nperm = nperm;
            S.iteration = iter;

            S.seed = D(dd).seed;
            S.partStart = 1;

            jobids{dd,pp,iter} = ...
                job2cluster( ...
                    @perm_les_func_streams_distance, ...
                    S);

            fprintf( ...
                '%s | density %d | iteration %d/%d | job %s\n', ...
                D(dd).name, ...
                S.perc, ...
                iter, ...
                niter, ...
                num2str(jobids{dd,pp,iter}));

        end

    end

end

%% Save job IDs

save( ...
    fullfile(root_dir,'distance_permutation_jobids.mat'), ...
    'jobids', ...
    'D', ...
    'dens_range', ...
    'nperm', ...
    'niter', ...
    'edge_mask_file', ...
    'bct_dir');

fprintf('\nSubmitted %d distance-matched jobs.\n',numel(jobids));

%% ------------------------------------------------------------------------
% COMBINE AND CONVERT DISTANCE-MATCHED PERMUTATIONS
% -------------------------------------------------------------------------

%% Combine and convert distance-matched permutations

clear
clc

root_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';

dens_range = [10 20 30 40 50];

niter = 20;
nperm_iter = 50;
nperm_total = niter * nperm_iter;

%% Datasets

D(1).name = 'TempSeq';
D(1).cache_file = ...
    fullfile(root_dir,'TempSeq','permutation_input.mat');
D(1).perm_dir = ...
    fullfile(root_dir,'TempSeq_Dist','Streamline_Permuts');
D(1).combined_dir = ...
    fullfile(root_dir,'TempSeq_Dist','Combined');
D(1).legacy_dir = ...
    fullfile(root_dir,'TempSeq_Dist','Permutation_Combined');

D(2).name = 'AudMemDem';
D(2).cache_file = ...
    fullfile(root_dir,'AudMemDem','permutation_input.mat');
D(2).perm_dir = ...
    fullfile(root_dir,'AudMemDem_Dist','Streamline_Permuts');
D(2).combined_dir = ...
    fullfile(root_dir,'AudMemDem_Dist','Combined');
D(2).legacy_dir = ...
    fullfile(root_dir,'AudMemDem_Dist','Permutation_Combined');

%% Combine iterations

for dd = 1:length(D)

    tmp = load(D(dd).cache_file,'cache');
    cache = tmp.cache;

    participant_ind = find(cache.healthy_selection);

    if ~exist(D(dd).combined_dir,'dir')
        mkdir(D(dd).combined_dir);
    end

    if ~exist(D(dd).legacy_dir,'dir')
        mkdir(D(dd).legacy_dir);
    end

    fprintf('\n%s\n',D(dd).name);

    for pp = 1:length(dens_range)

        density = dens_range(pp);

        rejection_rate = zeros(length(participant_ind),niter);

        fprintf('Density %d\n',density);

        for hh = 1:length(participant_ind)

            participant = participant_ind(hh);

            efficiency = zeros(nperm_total,1);
            clustering_coef = zeros(90,nperm_total);
            degree_les = zeros(90,nperm_total);

            affected_edges = zeros(nperm_total,1);
            affected_edges_bin = zeros(5,nperm_total);

            first_perm = 1;

            for iter = 1:niter

                filename = fullfile( ...
                    D(dd).perm_dir, ...
                    ['density_' num2str(density)], ...
                    ['permut_les_dist_mat' ...
                     num2str(participant) ...
                     '_iter_' ...
                     num2str(iter) ...
                     '.mat']);

                if ~exist(filename,'file')
                    error('Missing file: %s',filename);
                end

                tmp_perm = load(filename,'this_part_struct');
                P = tmp_perm.this_part_struct;

                if length(P.efficiency) ~= nperm_iter
                    error('Unexpected permutation count in %s.',filename);
                end

                last_perm = first_perm + nperm_iter - 1;

                efficiency(first_perm:last_perm) = ...
                    P.efficiency;

                clustering_coef(:,first_perm:last_perm) = ...
                    P.clustering_coef;

                degree_les(:,first_perm:last_perm) = ...
                    P.degree_les;

                affected_edges(first_perm:last_perm) = ...
                    P.affected_edges;

                affected_edges_bin(:,first_perm:last_perm) = ...
                    P.affected_edges_bin;

                rejection_rate(hh,iter) = ...
                    P.bookkeeping.rejection_rate;

                bookkeeping = P.bookkeeping;

                first_perm = last_perm + 1;

            end

            %% Save combined participant

            combined = struct();

            combined.participant_index = participant;
            combined.subject = cache.subjects{participant};
            combined.density = density;
            combined.nperm = nperm_total;

            combined.efficiency = efficiency;
            combined.clustering_coef = clustering_coef;
            combined.degree_les = degree_les;

            combined.affected_edges = affected_edges;
            combined.affected_edges_bin = affected_edges_bin;

            combined.total_raw_units = ...
                bookkeeping.total_raw_units;

            combined.observed_lesion_units = ...
                bookkeeping.observed_lesion_units;

            combined.lesion_fraction = ...
                bookkeeping.observed_lesion_units / ...
                bookkeeping.total_raw_units;

            combined.observed_affected_edges = ...
                bookkeeping.observed_affected_edges;

            combined.observed_affected_edges_per_bin = ...
                bookkeeping.observed_affected_edges_per_bin;

            combined.observed_lesion_units_per_bin = ...
                bookkeeping.observed_lesion_units_per_bin;

            combined.raw_units_per_bin = ...
                bookkeeping.raw_units_per_bin;

            combined.distance_percentiles = ...
                bookkeeping.distance_percentiles;

            combined.n_edges_per_bin = ...
                bookkeeping.n_edges_per_bin;

            combined.mean_rejection_rate = ...
                mean(rejection_rate(hh,:));

            outfile = fullfile( ...
                D(dd).combined_dir, ...
                ['permut_les_dist_mat' ...
                 num2str(participant) ...
                 '_density_' ...
                 num2str(density) ...
                 '_combined.mat']);

            save(outfile,'combined','-v7.3');

        end

        fprintf('  mean rejection rate = %.4f\n', ...
            mean(rejection_rate(:)));

        fprintf('  maximum rejection rate = %.4f\n', ...
            max(rejection_rate(:)));

    end

    %% Convert to format used by next analysis

    npart = size(cache.Raw_streams,3);

    for pp = 1:length(dens_range)

        density = dens_range(pp);

        perc_perm_str = struct();

        perc_perm_str.eff = ...
            nan(nperm_total,npart);

        perc_perm_str.clust = ...
            nan(90,nperm_total,npart);

        perc_perm_str.deg = ...
            nan(90,nperm_total,npart);

        for hh = 1:length(participant_ind)

            participant = participant_ind(hh);

            filename = fullfile( ...
                D(dd).combined_dir, ...
                ['permut_les_dist_mat' ...
                 num2str(participant) ...
                 '_density_' ...
                 num2str(density) ...
                 '_combined.mat']);

            if ~exist(filename,'file')
                error('Missing combined file: %s',filename);
            end

            tmp_perm = load(filename,'combined');
            P = tmp_perm.combined;

            if length(P.efficiency) ~= nperm_total
                error('Unexpected combined permutation count in %s.',filename);
            end

            perc_perm_str.eff(:,participant) = ...
                P.efficiency;

            perc_perm_str.clust(:,:,participant) = ...
                P.clustering_coef;

            perc_perm_str.deg(:,:,participant) = ...
                P.degree_les;

        end

        perc_perm_str.healthy_selection = ...
            cache.healthy_selection;

        perc_perm_str.subjects = ...
            cache.subjects;

        perc_perm_str.density = density;
        perc_perm_str.nperm = nperm_total;

        save( ...
            fullfile( ...
                D(dd).legacy_dir, ...
                ['perm_str_' num2str(density) '.mat']), ...
            'perc_perm_str', ...
            '-v7.3');

    end

end

fprintf('\nDistance-matched permutation combination complete.\n');


%% ------------------------------------------------------------------------
% CURRENT RICH-CLUB NULL PERMUTATIONS
% -------------------------------------------------------------------------
% IMPORTANT:
% The old rich-club permutation files were generated with the pre-review
% thresholding workflow. Recompute them with this section before enabling
% rich-club inference in 06_Combined_Graph_theory_analysis.m.
%
% This keeps the original rich-club null model (null_model_und_sign) and
% the original total of 10,000 null networks per participant/density, but
% uses the current cached raw matrices and the saved healthy-derived masks.

clear
clc

root_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';

function_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/Functions';

bct_dir = ...
    ['/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/' ...
     'Graph_Theory/BCT/2019_03_03_BCT'];

addpath(function_dir)

dens_range = [10 20 30 40 50];
nperm = 500;
niter = 20;

edge_mask_file = fullfile(root_dir,'combined_edge_masks.mat');

D = struct([]);

D(1).name = 'TempSeq';
D(1).cache_file = fullfile(root_dir,'TempSeq','permutation_input.mat');
D(1).outdir = fullfile(root_dir,'TempSeq_RichClub');
D(1).seed = 500000000;

D(2).name = 'AudMemDem';
D(2).cache_file = fullfile(root_dir,'AudMemDem','permutation_input.mat');
D(2).outdir = fullfile(root_dir,'AudMemDem_RichClub');
D(2).seed = 600000000;

clusterconfig('scheduler','cluster');
clusterconfig('long_running',2);
clusterconfig('slot',2);

for dd = 1:length(D)

    for pp = 1:length(dens_range)

        for iter = 1:niter

            S = struct();

            S.cache_file = D(dd).cache_file;
            S.edge_mask_file = edge_mask_file;
            S.outdir = D(dd).outdir;
            S.bct_dir = bct_dir;

            S.perc = dens_range(pp);
            S.nperm = nperm;
            S.iteration = iter;
            S.seed = D(dd).seed;

            job2cluster(@rich_club_perm_current,S);

        end
    end
end


%% Combine rich-club permutation iterations

clear
clc

root_dir = ...
    '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Cross_new';

dens_range = [10 20 30 40 50];

niter = 20;
nperm_iter = 500;
nperm_total = niter*nperm_iter;

D = struct([]);

D(1).name = 'TempSeq';
D(1).cache_file = fullfile(root_dir,'TempSeq','permutation_input.mat');
D(1).perm_dir = fullfile(root_dir,'TempSeq_RichClub','RichClub_Permuts');
D(1).outdir = fullfile(root_dir,'TempSeq_RichClub','Permutation_Combined');

D(2).name = 'AudMemDem';
D(2).cache_file = fullfile(root_dir,'AudMemDem','permutation_input.mat');
D(2).perm_dir = fullfile(root_dir,'AudMemDem_RichClub','RichClub_Permuts');
D(2).outdir = fullfile(root_dir,'AudMemDem_RichClub','Permutation_Combined');

for dd = 1:length(D)

    tmp = load(D(dd).cache_file,'cache');
    npart = size(tmp.cache.Raw_streams,3);

    if ~exist(D(dd).outdir,'dir')
        mkdir(D(dd).outdir);
    end

    for pp = 1:length(dens_range)

        density = dens_range(pp);
        RC_perm_comb = nan(nperm_total,90,npart);

        first_perm = 1;

        for iter = 1:niter

            filename = fullfile( ...
                D(dd).perm_dir, ...
                ['density_' num2str(density)], ...
                ['permut_richclub_iter_' num2str(iter) '.mat']);

            if ~exist(filename,'file')
                error('Missing rich-club permutation file: %s',filename);
            end

            tmp_perm = load(filename,'perm_mat');
            P = tmp_perm.perm_mat;

            if size(P,1) ~= nperm_iter || ...
                    size(P,2) ~= 90 || ...
                    size(P,3) ~= npart

                error('Unexpected dimensions in %s.',filename);
            end

            last_perm = first_perm + nperm_iter - 1;
            RC_perm_comb(first_perm:last_perm,:,:) = P;
            first_perm = last_perm + 1;

        end

        save( ...
            fullfile(D(dd).outdir, ...
            ['RC_comb_' num2str(density) '.mat']), ...
            'RC_perm_comb', ...
            '-v7.3');

        fprintf('%s | density %d%% complete\n',D(dd).name,density);

    end
end

disp('Rich-club permutation combination complete')
