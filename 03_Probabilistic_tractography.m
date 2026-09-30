%% FINAL MANUSCRIPT TRACTOGRAPHY NOTE
% The final graph-theory analysis uses cross-dataset age-sensitive clusters:
% an AudMemDem-derived cluster is applied to TempSeqAges and a TempSeqAges-
% derived cluster is applied to AudMemDem. Historical within-dataset and
% combined-cluster tractography blocks may also remain below for provenance.
% Participant 0001 was a pilot and was processed separately where loops begin at 2.

%% Probabilistic tractography
% Builds the baseline AAL90 connectivity matrices and the ageing-cluster-constrained matrices.
% Existing tractography commands, paths and participant loops are retained unchanged.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% SECTION 3 -- PROBABILISTIC TRACTOGRAPHY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running 2 forms: First creates a usual 90x90 connectivity matrix (AAL
% parcellation). Second only keeps seeds that touched at least one voxel in
% the ageing white matter cluster
%% PROBABILISTIC TRACTOGRAPHY - FSL
%BedpostX (using markov chain Monte Carlo) (slow, about 15 hours)

path1=dir([path_out '/0*']); %query the content of the folder "raw": it contains one folder for each subject
% A) Create and input directory for BedpostX (DTI001_BedpostX in this case)
% for each subject
for ii=2:length(path1) %loop across all subjects
    mkdir([path1(1).folder '/' path1(ii).name], [path1(ii).name '_BedpostX']);  %create a new folder for each subject (ii) and name it 'subjectID_Bedpostx'
end

%%  B) Add the following files to the folder
% 1) "data.nii.gz" which is the renamed version of "eddy_unwarped_images.nii.gz"
% 2) "nodif_brain_mask.nii.gz" mask (you already have it)
% 3) "bvals" obtained from initial nifti convertion (here we renamed "20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bval" as "bvals"
% 4) "bvecs" obtained from initial nifti convertion (here we renamed "20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bvec" as "bvecs"

path1=dir([path_out '/0*']); %query the content of the folder "raw": it contains one folder for each subject
path5=dir([path_out '0*']);
for gg=1:length(path1)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        path_bed=dir([path1(1).folder '/' path1(gg).name '/' path1(gg).name '_BedpostX' ]); %path to BedpostX folder
        path_datanii=[path1(1).folder '/' path1(gg).name '/eddy_unwarped_images.nii.gz']; %path(gg) to eddy unwarped images
        path_nodif=[path1(1).folder '/' path1(gg).name '/hifi_nodif_brain_mask.nii.gz']; %path to nodif_brain_mask
        path_sub=dir([path1(1).folder '/' path1(gg).name '/20*']); %path to each subject's main folder and relative content
        path_bvals=[path_sub(1).folder '/' path_sub(1).name];%path to each subject's .bval file
        path_bvecs= [path_sub(1).folder '/' path_sub(2).name];%path to each subject's .bvec file
        copyfile(path_datanii,[path_bed(1).folder '/data.nii.gz']); %copy eddy_unwarped_images and rename it as "data.nii.gz"
        delete(path_datanii); %remove eddy_unwarped_images from subject's main folder
        %remove oiginal eddy
        copyfile(path_nodif,[path_bed(1).folder '/nodif_brain_mask.nii.gz']); %copy hifi_nodif_brain_mask and rename it as nodif_brain_mask
        copyfile(path_bvals,[path_bed(1).folder '/bvals']); %copy .bval file and rename it as "bvals"
        copyfile(path_bvecs,[path_bed(1).folder '/bvecs']); %copy .bvec file and rename it as "bvecs"
        disp(gg)
    end
end

%%  C) calling terminal command with "bedpostx" operating in the desired directory (bedpostx) - slow (More then 70)

path1=dir([path_out '0*']); %query the content of the folder "raw": it contains one folder for each subject
for gg = 1:length(path1)
    path_bed_run = [path1(1).folder '/' path1(gg).name '/' path1(gg).name '_BedpostX' ]; %path to BedpostX folder
    %cd(path_bed); %change directory for each subject
    cmd = ['submit_to_cluster -q long.q -n 1 -p MINDLAB2023_MEG-AuditMemDement "/usr/local/fsl/bin/bedpostx ' path_bed_run '"']; %submit BedpostX to cluster
    system(cmd)
    disp(gg)
end
% cmd = 'bedpostx /projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/DTI001/DTI001_BedpostX';
%bedpostx command
%directory that you created (/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/DTI001/DTI001_BedpostX)
% system(cmd)
%%
path=dir([path_out '0*']); %path to subjects
for ii = 2:length(path) %over subjects
    path6=dir([path(ii).folder '/' path(ii).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        disp(ii)
        %going to subject's directory
        cd([path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX'])
        copyfile('nodif_brain_mask.nii.gz','nodif_brain.nii.gz') %copying mask file (just to rename it..)
        %actual command line for registration with flirt (FSL)
        cmd = 'flirt -in nodif_brain.nii.gz -ref /scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/AAL_Templates/MNI152_T1_2mm_brain.nii.gz -omat diff2stand.mat';
        system(cmd)
        %moving the output file into the bedpostX xfms subfolder
        copyfile('diff2stand.mat',[path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/diff2stand.mat'])
        load([path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/diff2stand.mat'],'-ascii') %loading the matrix computed by flirt
        %then, since we need to get the inverse matrix, it seems that we can simply compute it, as follows:
        stand2diff = inv(diff2stand);
        save([path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/stand2diff.mat'],'stand2diff','-ascii') %saving the inverse matrix
        disp(ii)
    end
end

%%  %%% STEP 1 %%%  Preparation
% Rescaling AAL parcels and white matter cluster to individual diffusion
% space per participant. Then creating a file that stores ROI sizes in
% individual diffusion space (for future normalisation)

% *Step 1* converting my waipoint (white matter cluster) mask to individual diffusion space
path = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/0*'); %path to subjects
ClustMNI = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/tbss/stats/Age_ClustOY2.nii.gz'; % Ageing cluster dir (IN MNI SPACE)
for ii = 2:length(path)
    path6=dir([path(ii).folder '/' path(ii).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        subjPath = fullfile(path(ii).folder, path(ii).name);
        bedDir = fullfile(subjPath,[path(ii).name,'_BedpostX.bedpostX']);
        xfm_mat = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/stand2diff.mat'];
        nodifmaskbrain = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/nodif_brain.nii.gz'];
        Mask_roi_dir = [path(ii).folder '/' path(ii).name '/ttestcluster_OY_ThisData'];
        if ~exist(Mask_roi_dir,'dir')
            mkdir(Mask_roi_dir);
        end
        parcel_path = ClustMNI;
        clust_path = [path(ii).folder '/' path(ii).name '/ttestcluster_OY_ThisData/ttestClust_waypoint_Comb.nii.gz'];
        cmd = sprintf(['/usr/local/fsl/bin/applywarp ' '--in=%s --ref=%s --out=%s --premat=%s --interp=nn'], parcel_path, nodifmaskbrain, clust_path, xfm_mat);
        system(cmd)
    end
    disp([num2str(ii) '/' num2str(length(path))])
end
%%
% *Step 2* converting AAL mask to individual diffusion space
path = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/0*'); %path to subjects
AAL_Masks = '/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/AAL_Templates/AAL_2mm_90ROIs';
AAL_Masks_dir = dir(AAL_Masks);
for ii = 2:length(path)
    path6=dir([path(ii).folder '/' path(ii).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        subjPath = fullfile(path(ii).folder, path(ii).name);
        bedDir = fullfile(subjPath,[path(ii).name,'_BedpostX.bedpostX']);
        xfm_mat = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/stand2diff.mat'];
        nodifmaskbrain = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/nodif_brain.nii.gz'];
        Mask_roi_dir = [path(ii).folder '/' path(ii).name '/AAL_ind_diff_space'];
        if ~exist(Mask_roi_dir,'dir')
            mkdir(Mask_roi_dir);
        end
        for jj = 1:length(AAL_Masks_dir)-2
            parcel_path = [AAL_Masks_dir(jj+2).folder '/' AAL_Masks_dir(jj+2).name];
            clust_path = [path(ii).folder '/' path(ii).name '/AAL_ind_diff_space/'  AAL_Masks_dir(jj+2).name];
            cmd = sprintf(['/usr/local/fsl/bin/applywarp ' '--in=%s --ref=%s --out=%s --premat=%s --interp=nn'], parcel_path, nodifmaskbrain, clust_path, xfm_mat);
            system(cmd)
        end
    end
    disp([num2str(ii) '/' num2str(length(path))])
end
%%
result = [];
for i = 1:9
    result = [result,i];
    for j = 0:9
        num = i*10+j
        result = [result,num]
    end
end
file_list_numbers = result(1:90)


reshap_order = [];
for i = 1:90
    posit = find(file_list_numbers == i);
    reshap_order = [reshap_order,posit];
end
%%
% *Step 3* creating directories for individual masks.txt files (directories to all aal parcels for probabilistic tractography)
path = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/0*');
for ii = 2:length(path) %over subjects
    path6=dir([path(ii).folder '/' path(ii).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        file_name = 'masks.txt';
        masks_path = fullfile([path(ii).folder '/' path(ii).name], file_name);
        fid = fopen(masks_path, 'w'); % just creating the empty txt to adjust later
        fclose(fid);

        masks_folder = [path(ii).folder '/' path(ii).name '/AAL_ind_diff_space'];
        curr_dir = [path(ii).folder '/' path(ii).name];
        list_AAL_files = dir(masks_folder); % get the struckt with all the
        output_file = fullfile(curr_dir, 'masks.txt');
        fid = fopen (output_file, 'w');

        for k = reshap_order %this loop ads directories of every AAL parcellation file to the txt list
            curr_ind = k+2;
            fullpath = fullfile(masks_folder,list_AAL_files(curr_ind).name);
            fprintf(fid, '%s\n', fullpath);
        end
    end
end
%%
NumOfValid_Subj = 0;
for ii = 1:length(path) %over subjects
    path6=dir([path(ii).folder '/' path(ii).name '/20*']); %Check if folder has data
    if ~isempty(path6)
        NumOfValid_Subj=NumOfValid_Subj+1;
    end
end
%%
% *Step 4* Creating a file with all ROI sizes for future normalisation
ROI_sizes = zeros(90,NumOfValid_Subj);
path = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/0*'); %path to subjects
order = [1:2:90 90:-2:2]; %to reshape the matrix symmetricaly
subjCounter = 0;
for ii = 1:length(path) %over subjects
    path6=dir([path(ii).folder '/' path(ii).name '/20*']); %Check if folder has data
    if ~isempty(path6)
        subjCounter = subjCounter+1;
        ROI_sizes_temp = zeros(90,1);
        for pp = 1:90 %over ROI1
            %loading image ROI1
            ROI1 = load_untouch_nii([path(ii).folder '/' path(ii).name '/AAL_ind_diff_space/' num2str(pp) '.nii.gz']);
            sizeROI1 = length(find(ROI1.img==1)); %calculating size of Rx = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - bluecolormap(bluewhitered_PD(0,x))OI1 (number of voxels)
            ROI_sizes_temp(pp,1) = sizeROI1;
            disp(['Subject = ' path(ii).name ' - ROI1 = ' num2str(pp)])
        end
        ROI_sizes (:,subjCounter) = ROI_sizes_temp(order,1);
    end
end

%%
Mask_roi_dir = ['/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/GraphTheory_Ageing_Combined'];
if ~exist(Mask_roi_dir,'dir')
    mkdir(Mask_roi_dir);
end
%%
tmp = load('/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/Graph_Theory_Age_Diff/ROI_sizes.mat');
ROI_sizes = tmp.ROI_sizes;

%%  %%% STEP 2 %%%  Probabilistic tractography itself
% Tractography 1: creating 90 by 90 connectivity matirx with AAL90 parcellation

mask_label=2; %flag for mask type: 1=3559 parcels; 2=2mm AAL
streamlines=5000; %set streamlines number (5000 streamlines for AAL 2mm)


path = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/0*'); %path to subjects
wpNames = {'Age_ttest'};
wp = '/scratch7/MINDLAB2021_MEG-TempSeqAges/Nikita/DTI/tbss/stats/perm_res_final/Age_MainEff_Clust.nii.gz';
for ii = 1:length(path) %over subjects
    path6=dir([path(ii).folder '/' path(ii).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        if mask_label ==1
            mask = ['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/masks.txt']; %3559 parcels
            outdir = [path(ii).folder '/' path(ii).name '/Tractography/parce3559_' num2str(streamlines) 'stream']; %output directory
        elseif mask_label==2
            mask = [path(ii).folder '/' path(ii).name '/masks.txt'];%2-mm AAL
            %UPDATE THE DIRECTORY OF THE MASK
        end
        xfm_mat = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/stand2diff.mat']; %registration matrix
        invxfm_mat = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/diff2stand.mat']; %registration matrix (inverse)
        mergedk = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/merged']; %merged files from bedpostX
        nodifmask = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/nodif_brain_mask']; %brain mask from bedpostX

        wp = [path(ii).folder '/' path(ii).name '/ttestcluster_Combined/ttestClust_waypoint_Comb.nii.gz'];
        foldnam = wpNames{1};
        outdir = [path(ii).folder '/' path(ii).name '/Tractography_AgeDiff/NoCluster/' foldnam '/AAL90_' num2str(streamlines) 'stream'];
        if ~exist(outdir,'dir')
            mkdir(outdir);
        end
        %actual line for probtrackx - FSL
        cmd = ['submit_to_cluster -q long.q -n 2 -p MINDLAB2023_MEG-AuditMemDement "/usr/local/fsl/bin/probtrackx2 --network -x ' mask ' -l --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P ' num2str(streamlines) ' --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0' ' --forcedir --opd -s ' mergedk ' -m ' nodifmask ' --dir=' outdir '"'];
        system(cmd)
        disp([num2str(ii) '/' num2str(length(path))])
    end
end
%%
% Tractography 2: Counting only streamlines passing through the white matter cluster
% (condition: seed passes through at least one voxel of the cluster along its way from one parcel to the next)

mask_label=2; %flag for mask type: 1=3559 parcels; 2=2mm AAL
streamlines=5000; %set streamlines number (5000 streamlines for AAL 2mm)

path = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/0*'); %path to subjects
wpNames = {'Age_ttest'};
for ii = 2:length(path) %over subjects
    path6=dir([path(ii).folder '/' path(ii).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        if mask_label ==1
            mask = ['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/masks.txt']; %3559 parcels
            outdir = [path(ii).folder '/' path(ii).name '/Tractography/parce3559_' num2str(streamlines) 'stream']; %output directory
        elseif mask_label==2
            mask = [path(ii).folder '/' path(ii).name '/masks.txt'];%2-mm AAL
            %UPDATE THE DIRECTORY OF THE MASK
        end
        xfm_mat = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/stand2diff.mat']; %registration matrix
        invxfm_mat = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/diff2stand.mat']; %registration matrix (inverse)
        mergedk = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/merged']; %merged files from bedpostX
        nodifmask = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/nodif_brain_mask']; %brain mask from bedpostX
        %looping over cluster types
        wp = [path(ii).folder '/' path(ii).name '/ttestcluster_OY_ThisData/ttestClust_waypoint_Comb.nii.gz'];
        foldnam = wpNames{1};
        outdir = [path(ii).folder '/' path(ii).name '/Tractography_AgeDiff/WM_Cluster_ThisDat/' foldnam '/AAL90_' num2str(streamlines) 'stream'];
        if ~exist(outdir,'dir')
            mkdir(outdir);
        end
        %actual line for probtrackx - FSL
        cmd = ['submit_to_cluster -q long.q -n 2 -p MINDLAB2023_MEG-AuditMemDement "/usr/local/fsl/bin/probtrackx2 --network -x ' mask ' --waypoints=' wp ' --waycond=OR' ' -l --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P ' num2str(streamlines) ' --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0' ' --forcedir --opd -s ' mergedk ' -m ' nodifmask ' --dir=' outdir '"'];
        system(cmd)
        disp([num2str(ii) '/' num2str(length(path))])
    end
end

%% clear the workspace before the next section
clear all
