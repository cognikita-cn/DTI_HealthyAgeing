%% TBSS and cluster-based permutation analysis
% TBSS preprocessing, voxel-wise age-group statistics and cluster-size permutation thresholding.
% Existing paths, parameters and analysis code are retained unchanged.

%% TBSS - Tract-Based Spatial Statistics

%To compare FA values between two groups of subjects
%ALL the FA files (one per subject) must be in the SAME directory

%Make the directory for TBSS analysis
path_tbss='/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/tbss'; %set the path to the folder where you will save your converted data
mkdir(path_tbss)
%%
path5=dir([path_out '0*']);
%0) Copy and rename *_FA.nii.gz files from each subject into a new folder, where tbss will be run - fast
for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        path_FA=dir([path5(gg).folder '/' path5(gg).name '/*FA.nii.gz']); %loop across each subject and get the content of each subject's FA file
        copyfile([path_FA(1).folder '/' path_FA(1).name],[path_tbss '/' path5(gg).name  '_FA.nii.gz']); %copy acqparams.txt from DTI001 in every subject's folder and rename the file
        disp(gg)
    end
end

%% 1) Removing likely outliers by removing brain-edge artefacts and the zero end slices (tbss_1_preproc) - FAST process

% tbss_1_preproc *nii.gz
%It creates a new folder called FA, containing the sub-directory "origdata", with the original images
cd(path_tbss); %set tbss folder as current path
cmd = ['tbss_1_preproc *nii.gz']; %run tbss preprocessing
system(cmd)
% index.html shows the slices of every single subject

%% 2) Non-linear registration: aligning all the FA data across subjects (tbss_2_reg ) - moderately to very SLOW process (see options)
% tbss_2_reg
% It estimates warping parameters to standardize space (later applied in step 3)
% Two ways to do it:
% a) align to FSL's "FMRIB58_FA" Template (about 10 min/participant) - RECOMMENDED by FSL guys. Use -T flag
% b) automatic search for the most representative image to be used as a template across subjects (days/weeks). Use -t flag to choose your own target image

% 2.1) Creating new FMRIB58_FA template with the wanted resolution (flirt)
% - not necessary
flag_2mm=0; %Select template resolution: 1=2mm or 0=1mm (WE CHOSE 1mm, as suggested by FSL)
if flag_2mm~=1
    varT='T';
else
    template_path='/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/tbss/template';% path to the templates
    cmd=['flirt -in ' template_path '/FMRIB58_FA_1mm.nii.gz -ref ' template_path '/FMRIB58_FA_1mm.nii.gz -out ' template_path '/FMRIB58_FA_2mm.nii.gz -applyisoxfm 2']; %convert the 1mm template into a 2mm one (matching the DTI acquisition parameters)
    system(cmd);
    varT=['t ' template_path '/FMRIB58_FA_2mm.nii.gz']; %select the desired template by using -t and path to template
end
%%
%2.2)Running non-linear registration (tbss_2_reg)
cd(path_tbss); %set tbss folder as current path
cmd=['submit_to_cluster -q short.q -n 1 -p MINDLAB2023_MEG-AuditMemDement "/usr/local/fsl/bin/tbss_2_reg -' varT '"']; %run non-linear registration
system(cmd);

%% 3) Post-registration: applying the previous registration to take all subjects into 1x1x1mm standard space (tbss_3_postreg) - relatively fast
%Transforming the original FA image into MNI152 (1x1x1mm) space and merges all of the subjects' images into a single 4D image called "all_FA" (saved in "stats" folder) and calculates the mean of all images (mean_FA), which is used to create the "mean_FA_skeleton"

%OBS! It's always better to VISUALLY INSPECT that the skeleton is well aligned with MNI152 image!! --> FSLview or FSLeyes
% Options:
% a) -S : use the mean across all subjects as a skeleton (recommended).
% b) -T : use the FMRIB58_FA mean FA image and its derived skeleton.
cd(path_tbss); %set tbss folder as current path
cmd=['tbss_3_postreg -S'];
system(cmd);
%  Check MATLAB command window for info about the status
%when visualizing the output on fsl, remember to add also the MNI152
%template and change the colourscale for the skeleton

%% 4) Project the pre-aligned data into the skeleton (tbss_4_prestats)

%Thresholding the mean FA skeleton image at the chosen threshold (0.3).
% the script takes the "all_FA" image (containing all subjects' aligned FA data) and, for each subject, projects the FA data onto the mean FA skeleton.
%The output is a binary skeleton mask (4D image) with the projected skeletonised FA data = set of voxels that will be used for the actual statistics
cd(path_tbss); %set tbss folder as current path
cmd=['tbss_4_prestats 0.2']; %0.2 is a value that generally works well
system(cmd);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% SECTION 2 -- VOXEL-WISE statistics + clust-based permutations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Logical steps:
% - Sending permutations to the cluster over several iterations for
% parallelisation.
% - Defining permutations based cluster threshold
% - Running ttest on the original data contrasting old and young FA values
% per voxel
% - Defining significantly big clusters
% - Exporting them to a NII file for future probabilistic tractography and
% visual inspection in FSLeyes

%% %%% STEP 1 %%% Cluster based permutations and defining cluster size threshold


% !!Preparation: do once
% Parallelizing and sending jobs with permutations function to the cluster

S = []; % creating a structure to feed to the permutation function

S.nperm = 100; % amount of permutations per iteration
S.iteration = 0; % preallocating iteration value
S.path_out='/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/'; %set the path to the folder where you will save your converted data
S.fname = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/tbss/stats/all_FA_skeletonised.nii.gz';
% Sending jobs to the cluster
addpath('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/Functions') % path to the scripts folder
for iter = 1:100
    clusterconfig('scheduler', 'cluster'); %set the cluster
    % clusterconfig('scheduler', 'none'); %set locally
    clusterconfig('long_running', 0); %set the most typical cue (see our cluster guide in labook for details)
    clusterconfig('slot', 2); %set amount of memory; between 1 and 12 (each slot is 8gb of ram)
    S.iteration = iter;
    addpath('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/Functions') % path to the scripts folder
    jobid = job2cluster(@ttest_perm_func,S)
end

%% %%% STEP 2 %%%  ttest on the original data

%uploading the allFa
fname = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/tbss/stats/all_FA_skeletonised.nii.gz';
skel1 = load_nii(fname); %Load all_FA_skeletonised file (4D file)
refNii = skel1;
skel1=skel1.img;
SS = size(skel1(:,:,:,1));
ind = find(skel1(:,:,:,1) ~= 0); %getting indices of non-0 voxels (indeces that only belong to brain areas, not empty space)
[i1,i2,i3] = ind2sub(size(skel1(:,:,:,1)),ind); %reshaping the vector into a 3D matrix shaped as the original data (skel) = get one vector for each dimension

%%
participant_info = readtable('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/tbss/Descript_Dat/Participant_categories_AuditMemDement_ages.xlsx'); % loading behavioural data
%
path_out='/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/'; %set the path to the folder where you will save your converted data

path5=dir([path_out '0*']);
DTI = zeros(length(participant_info.Group),1);
%0) Copy and rename *_FA.nii.gz files from each subject into a new folder, where tbss will be run - fast
a = 0;
for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        a = a+1;
        DTI(gg,1) = 1;
    end
end
participant_info.DTI = DTI;
%
%
Valid_participant_info  = participant_info(participant_info.DTI==1,:);
Y_Selection1 = ismember(Valid_participant_info.Group,{'A'});
O_Selection1 = ismember(Valid_participant_info.Group,{'B'});
%O_Selection = ismember(Valid_participant_info.Group,{'B','B_MCI'});
D_Selection = ismember(Valid_participant_info.Group,{'C'});
%%
SubSet = Y_Selection|O_Selection;
Y_Selection = Y_Selection(SubSet);
O_Selection = O_Selection(SubSet);
%%
figure
hist(Valid_participant_info.Age(Y_Selection))
hold on
hist(Valid_participant_info.Age(O_Selection))
hold off
%%
participant_info = readtable('/aux/MINDLAB2021_MEG-TempSeqAges/Nikita/TSA2021_Nikita.xlsx'); % loading behavioural data

%starting up some of the functions by LB for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl);


%Converting 0 no NaN
participant_info.WMCombined(participant_info.WMCombined==0)=NaN;
participant_info.STMMelody(participant_info.STMMelody==0)=NaN;
participant_info.STMRhythm(participant_info.STMRhythm==0)=NaN;

ValidRows = ~isnan(participant_info.WMCombined); %IF you want to test for WM relations afterwards
%ValidRows = participant_info.Subject~=1;
Valid_dat = participant_info(ValidRows,:); %defining a dataset with participants that have all the behavioural tests
Young = Valid_dat.Subject(Valid_dat.Age < 27);
Old = Valid_dat.Subject(Valid_dat.Age > 55);
Y_Selection = ismember(Valid_dat.Subject,Young);
O_Selection = ismember(Valid_dat.Subject,Old);
%%
figure
histogram(Valid_dat.Age(O_Selection))
hold on
histogram(Valid_participant_info.Age(O_Selection1), 'FaceColor', [0.85 0.33 0.10])
hold off
%%
%starting up some of the functions by LB for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl);

%subsetting only relevant images
%skel = skel1(:,:,:,:); % need this step to dublicate the skel file, otherwise running out of memory when trying to load it at later stages
skel = skel1(:,:,:,SubSet);
%%

% setting tvalue threshold for future binarisation
num_part = size(skel,4);
df_age = num_part - 2;
t_Thresh_age = tinv(1-0.005/2,df_age);

%%

% preallocating vectors to store pvalue and tvalue per voxel
P = zeros(length(i1),1);
TVAL = zeros(length(i1),1);
for nn = 1:length(ind) %Loop across voxels
    a = squeeze(skel(i1(nn),i2(nn),i3(nn),Y_Selection)); %take 1st group of subjects
    b = squeeze(skel(i1(nn),i2(nn),i3(nn),O_Selection)); %take 2nd group of subjects
    [~,p,~,stats] = ttest2(a,b,'Vartype','equal'); %compute t-test between the two groups
    P(nn,1) = p; %store p-values in P
    TVAL(nn,1) = stats.tstat; %store t-values in TVAL
    %disp([num2str(nn) '/' num2str(length(ind))])
end
%reshaping them into 3D matrices
TVAL2 = zeros(SS(1),SS(2),SS(3));
P2 = zeros(SS(1),SS(2),SS(3));
for ii =1:length(i1)
    TVAL2(i1(ii),i2(ii),i3(ii)) = TVAL(ii);
    P2(i1(ii),i2(ii),i3(ii)) = P(ii);
end

% Using bwconncomp -- efficient function for detecting clusters in 2d or 3d
% data
conn = 26; % for detecting neighbouring binarised neighbours in 3d space

Bin_PVAL_Age_WM3D = zeros(size((TVAL2))); % preallocating zeros of same size
Bin_PVAL_Age_WM3D((TVAL2>t_Thresh_age)) = 1; % binarising based on the t-value threshold defined earlier
CC_Age = bwconncomp(Bin_PVAL_Age_WM3D,conn); % detecting clusters

a_Age = CC_Age.PixelIdxList; % extracting cell list with all of my clusters
lengths_Age = cellfun(@length,a_Age); % Creating a vector with corresponding size for each of the detected vectors
disp('DONE')
%%  %%% STEP 3 %%%  defining permutation based cluster size threshold
my_dir = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/tbss_combined/stats/Age_WM_ClustPerm_Combined'; % directory with permutation results

% setting up permutation matrices (combining several iterations into one matrix)

AgeEffClust = load([my_dir, '/Permutation_res_ttest_Age_Clusters_iteration_',num2str(1),'.mat']);
names = fieldnames(AgeEffClust);
AgePerms = AgeEffClust.(names{1});

for i = 2:100
    AgeEffClust = load([my_dir, '/Permutation_res_ttest_Age_Clusters_iteration_',num2str(i),'.mat']);
    names = fieldnames(AgeEffClust);
    AgePerms_iter = AgeEffClust.(names{1});
    AgePerms = [AgePerms;AgePerms_iter];
end

% Thresholding
AgePerms_T = AgePerms(:,1);% looking at the biggest cluster per permutation
AgeThreshold = prctile(AgePerms_T,95);% finding 95 percentile cluster size
%%
Orig_stat = load('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/tbss_combined/stats/Age_ttest_clust_Origin.mat');
CC_Age = Orig_stat.CC_Age;
a_Age = CC_Age.PixelIdxList; % extracting cell list with all of my clusters
lengths_Age = cellfun(@length,a_Age); % Creating a vector with corresponding size for each of the detected vectors


%%

Sign_Age_Clust = a_Age(find(lengths_Age>AgeThreshold)); %finding sigificantly big clusters in the original data

% sanity check: looking if mean FA value among my significant clusters
% differs between age groups. Take into account !circularity!: these clusters
% were already identified using a ttest between age groups, so that is
% purely a sanity check to see that nothing went wrong
nSub = size(skel,4);
nClust_Age = numel(Sign_Age_Clust);
meanFA_Age = nan(nSub,nClust_Age);
for i = 1:nClust_Age
    idx = Sign_Age_Clust{i};
    faVals = reshape(skel(idx+(0:nSub-1)*prod(SS)),[],nSub);
    meanFA_Age(:,i) = sum(faVals,1)';
end
Grand_Mean_FA = mean(meanFA_Age,2);
[~,p,~,stats] = ttest2(Grand_Mean_FA(Y_Selection),Grand_Mean_FA(O_Selection))

% Saving significant clusters as nifti images

% combining all voxel indices from my significant clusters into one list,
% just to save my significant clusters as one file. Inspect in FSLeyes for
% visual check, they remain as independent clusters spacially.
all_idx_Age = Sign_Age_Clust{1};

for i = 2:length(Sign_Age_Clust)
    idx = Sign_Age_Clust{i};
    all_idx_Age = [all_idx_Age;idx];
end
%%
% Saving my significant white matter clusters to a nifti image
maskVol_age = zeros (SS,'uint8');
maskVol_age(all_idx_Age) = 1;

%refNii = load_nii(fname);
refNii.img = maskVol_age;
refNii.hdr.dime.dim(1) = 3;
refNii.hdr.dime.dim(2:4) = SS;
refNii.hdr.dime.dim(5) = 1;
refNii.hdr.dime.datatype = 2;
refNii.hdr.dime.bitpix = 8;
%
save_nii(refNii,'Age_ClustOY2.nii.gz');
