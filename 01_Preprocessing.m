%% DTI preprocessing
% Original preprocessing workflow for the DTI Healthy Ageing project.
% Cluster-specific paths, manual setup notes and run ranges are retained intentionally.


%OBS! To set-up the code to running.. only on the first time, follow these steps (source:  labbook)
% 1. open terminal window and type:
%       gedit .bashrc
% 2. copy the following line to the beginning of the file, save and close the window:
%       source /usr/local/common/meeg-cfin/configurations/setup_environment.sh
% 3. in same terminal window, type:
%       gedit .bash_profile
% 4. copy following 3 lines into the file, save and close window:
%       if [ -f $HOME/.bashrc ]; then
%           source $HOME/.bashrc
%       fi

%OBS! before running the code below you need to close matlab, open the terminal and write: 'use anaconda', then open matlab and run the code

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% SECTION 1 -- DTI data PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIPELINE OF ALL SUBJECTS - DTI - FSL %%


%Relavant SCRIPTS at:
'/projects/MINDLAB2021_MEG-TempSeqAges/scripts/Nikita/Ageing_conn_scripts';
path_out='/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/'; %set the path to the folder where you will save your converted data


%% INITIAL CONVERSION FROM DICOM (RAW) FILE TO NIFTI (dcm2nii) - fast

%mkdir('/scratch7/MINDLAB2022_MEG-EncodingMusicSeq/Nikita/DTI/') %creates the directory for storign analysis results
path1=dir([path '/0*']); %query the content of the folder "raw": it contains one folder for each subject
for ii=1:length(path1) %loop across all subjects
    path2=dir([path1(ii).folder '/' path1(ii).name '/20*']); %for each subject, search folders named '2018*'. There should be 1 or 2 folders per participant: one contains the folder 'MR' and 'SR' the other one contains the folder 'MEG'. We need the folder 'MR'.
    outdir = [path_out path1(ii).name];
    if ~exist(outdir,'dir')
        mkdir(path_out, path1(ii).name);  %for each subject (ii), create a new folder in (DTI_Portis) and call it with the same name as the subject's (path1(ii).name)
    end
    for dd=1:length(path2)
        path3=dir([path2(dd).folder '/' path2(dd).name '/MR/*AP']); %choose 'MR' rather than 'SR' - AP
        if length(path3) > 1
            warning(['Subj ' path1(ii).name ' has more than 1 DTI AP'])
        end
        path4=dir([path2(dd).folder '/' path2(dd).name '/MR/*PA']); %choose 'MR' rather than 'SR' - PA
        if length(path4) > 1
            warning(['Subj ' path1(ii).name ' has more than 1 DTI PA'])
        end
        if ~isempty(path3)
            pathAP=[path3.folder '/' path3.name '/files']; %set path to diff data (AP)
            %Convert AP to .nii.gz
            cmd = ['dcm2nii -o ' path_out path1(ii).name ' ' pathAP]; %local
            system(cmd) %send command to the cluster
        end
        if ~isempty(path4)
            pathPA=[path4.folder '/' path4.name '/files'];%set path to PA
            %Convert PA
            cmd = ['dcm2nii -o ' path_out path1(ii).name ' ' pathPA];
            system(cmd)
        end
    end
    disp(ii)
end


%% Getting ready for TOPUP

cd(path_out);

% 1) NODIF (fslroi) - fast
%creating reference volume (nodif) based on the first image of the DTI data
path5=dir([path_out '0*']);
for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder (specify 20* to avoid hidden files)
    if ~isempty(path6)
        cd([path5(gg).folder '/' path5(gg).name]); %change directory for each subject, so that nodif is saved in that of the corresponding subject
        cmd = ['fslroi ' path6(1).folder '/' path6(3).name ' nodif 0 1']; %apply fslroi to AP, which is the 3rd file in the folder
        system(cmd)
        disp(gg)
    end
end

%fslroi command
%input DTI data (AP)
%output name of the reference volume (nodif)
%0 (minimum index of the volumes within the nifti file (indexing as in python=starting from 0))
%1 (number of images to be taken (in this case only first))


%%
% 2) NODIF_PA (copyfile) - fast
%renaming PA as follows:
%1)copy the DTI PA file (20180123_102505cmrrmbep2ddiffmultidirPASubjectNo0001s023a1001.nii.gz)
%2)rename it as "nodif_PA.nii.gz"
for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        cd([path5(gg).folder '/' path5(gg).name]); %change directory for each subject, so that a nodif_PA.nii.gz is created for each subject
        copyfile([path6(4).folder '/' path6(4).name],[path6(4).folder '/nodif_PA.nii.gz']); %copy PA (4th file in the folder) and rename it
        disp(gg)
    end
end
%%
% 3) AP_PA_b0 (fslmerge) - fast
%merging AP with PA data with the reference volume (nodif)
for mm=1:length(path5)
    path6=dir([path5(mm).folder '/' path5(mm).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        cd([path5(mm).folder '/' path5(mm).name]); %change directory for each subject, so that a nodif_PA.nii.gz is created for each subject
        cmd = ['fslmerge -t ' path6(4).folder '/AP_PA_b0 nodif nodif_PA'];
        system(cmd)
        disp(mm)
    end
end
%fslmerge command
%-t (parameter of fslmerge)
%output name (AP_PA_b0)
%reference AP (nodif)
%reference PA (nodif_PA)

%% Creating acquisition parameters file (acqparams.txt) - To be run ONLY for the FIRST SUBJECT!!! - fast
%txt file with information about acquisition parameters (the same for all subjects)
%In this case: 2x4 matrix
%row 1: 0 -1 (meaning AP) 0 (from the beginning of time) 0.104 (echo time in seconds (104ms))
%!!!!!! NEED TO FIX WHERE THE FILE IS BEING CREATED AND CLEAR THE FIRST
%LINE OF THE CREATED MATRIX, IT SHOULD ONLY BE NUMBERS WITH SPACES

cd(path_out);
matrixacqpar = [0 -1 0 0.065; 0 1 0 0.065]; %creating double matrix with values
dlmwrite('acqparams.txt', matrixacqpar,'delimiter',' ')
%t = table(matrixacqpar); %converting matrix to table
%writetable(t,'acqparams.txt') %saving table as .txt file (acqparams.txt MUST be the name)

%% Copying acqparams.txt to each subject's folder - To be run for ALL THE OTHER SUBJECTS

path_acqparams='/scratch7/MINDLAB2023_MEG-AuditMemDement/nikita/DTI/acqparams.txt';
%set the path to the folder where you will save your converted data
path5=dir([path_out '0*']);
for gg=99:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %here it is just to check for empty folders and not add anything to them
    if ~isempty(path6)
        copyfile(path_acqparams,[path5(gg).folder '/' path5(gg).name '/acqparams.txt']); %copy acqparams.txt from DTI001 in every subject's folder
        disp(gg)
    end
end

%% Call to TOPUP command (topup) - quite fast

%set the path to the folder where you will save your converted data
path5=dir([path_out '0*']);
for gg=99:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    %cmd = 'topup --imain=AP_PA_b0 --datain=acqparams.txt --config=b02b0.cnf --out=topup_AP_PA_b0 --iout=topup_AP_PA_bo_iout --fout=topup_AP_PA_bo_fout'; %apply topup to every subject
    if ~isempty(path6)
        cmd = ['submit_to_cluster -q all.q -n 1 -p MINDLAB2023_MEG-AuditMemDement "/usr/local/fsl/bin/topup --imain=' path6(1).folder '/AP_PA_b0 --datain=' path6(1).folder '/acqparams.txt --config=b02b0.cnf --out=' path6(1).folder '/topup_AP_PA_b0 --iout=' path6(1).folder '/topup_AP_PA_b0_iout --fout=' path6(1).folder '/topup_AP_PA_b0_fout"']; %apply topup to every subject
        system(cmd)
        disp(gg)
    end
end
%cmd = 'topup --imain=AP_PA_b0 --datain=acqparams.txt --config=b02b0.cnf --out=topup_AP_PA_b0 --iout=topup_AP_PA_bo_iout --fout=topup_AP_PA_bo_fout';
%topup command
%DTI data after fslmerg (--imain=AP_PA_b0)
%acquisition parameters (--datain=acqparams.txt)
%configuration file specifying command line arguments (--datain=acqparams.txt); this file is already in the FSL directory and it contains some (default) specifications
%output file (--out=topup_AP_PA_b0)
%output file with unwarped images (--iout=topup_AP_PA_bo_iout)
%output file with field (Hz) (--fout=topup_AP_PA_bo_fout)
%% Generating a brain mask from the corrected b0 (fslmaths, bet) - fast
path5=dir([path_out '0*']);

for gg=2:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        cd([path5(gg).folder '/' path5(gg).name]); %change directory for each subject, so that a nodif_PA.nii.gz is created for each subject
        cmd = 'fslmaths topup_AP_PA_b0_iout.nii.gz -Tmean hifi_nodif'; %apply fslmaths to create a mask from the corrected b0
        system(cmd)
        disp(gg)
    end
end
% cmd = 'fslmaths topup_AP_PA_b0_iout.nii.gz -Tmean hifi_nodif';
%fslmaths command
%input (topup_AP_PA_b0_iout.nii.gz)
%specification -Tmean
%output (hifi_nodif)
% system(cmd)
%%
%extracting the brain from b0 - fast
for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        cd([path5(gg).folder '/' path5(gg).name]); %change directory for each subject, so that a nodif_PA.nii.gz is created for each subject
        cmd = 'bet hifi_nodif hifi_nodif_brain -m -f 0.2'; %apply BET to extract the brain from b0
        system(cmd)
        disp(gg)
    end
end
% cmd = 'bet hifi_nodif hifi_nodif_brain -m -f 0.2';
%bet command
%input brain mask (hifi_nodif)
%output (hifi_nodif_brain)
%parameters (-m -f 0.2); -f means "fraction intensity threshold"
% system(cmd)


%% New automatized way of creating index files

raw_path='/projects/MINDLAB2023_MEG-AuditMemDement/raw'; %path to the folder containing all the raw data
%mkdir('/scratch7/MINDLAB2022_MEG-EncodingMusicSeq/Nikita/DTI/') %creates the directory for storign analysis results
path1=dir([raw_path '/0*']); %query the content of the folder "raw": it contains one folder for each subject
for ii=1:length(path1) %loop across all subjects
    path2=dir([path1(ii).folder '/' path1(ii).name '/20*']); %for each subject, search folders named '2018*'. There should be 1 or 2 folders per participant: one contains the folder 'MR' and 'SR' the other one contains the folder 'MEG'. We need the folder 'MR'.
    for dd=1:length(path2)
        path3=dir([path2(dd).folder '/' path2(dd).name '/MR/*AP/files/PRO*']); %choose 'MR' rather than 'SR' - AP
        if ~isempty(path3)
            thispath = dir([path_out '/' path1(ii).name]);
            path_string = thispath(1).folder;
            volume_count = length(path3);
            ind = ones(volume_count,1);
            dlmwrite([path_string '/index.txt'], ind,'delimiter',' ')
        end
    end
end


%% Correcting for EDDY currents (eddy) - slow: about 10-15 hours/subject
%(currents generated in the MRI machine because of a rapid change of the magnetic field direction during the acquisition: echo planar images are acquired rapidly in different orientations)

path5=dir([path_out '0*']);

for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        cmd = ['submit_to_cluster -q long.q -n 1 -p MINDLAB2023_MEG-AuditMemDement "/usr/local/fsl/bin/eddy --imain=' path6(3).folder '/' path6(3).name ' --mask=' path6(1).folder '/hifi_nodif_brain_mask --index=' path6(1).folder '/index.txt --acqp=' path6(1).folder '/acqparams.txt --bvecs=' path6(2).folder '/' path6(2).name ' --bvals=' path6(1).folder '/' path6(1).name ' --fwhm=0 --topup=' path6(1).folder '/topup_AP_PA_b0 --flm=quadratic --out=' path6(3).folder '/eddy_unwarped_images"'];
        system(cmd)
        disp(gg)
    end
end
% cmd = 'eddy --imain=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.nii.gz --mask=hifi_nodif_brain_mask --index=index.txt --acqp=acqparams.txt --bvecs=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bvec --bvals=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bval --fwhm=0 --topup=topup_AP_PA_b0 --flm=quadratic --out=eddy_unwarped_images';
%command eddy
%main input DTI data from initial nifti convertion (--imain=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.nii.gz)
%mask outputted from BET (hifi_nodif_brain_mask); NOTE that the name of the output that we use is "hifi_nodif_brain_mask" and not "hifi_nodif_brain"
%file with indices of the volumes (--index=index.txt)
%acquisition parameters (--acqp=acqparams.txt)
%bvecs outputted from the initial nifti convertion (--bvecs=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bvec); bvecs indicate the direction of diffusion
%bvals outputted from the initial nifti convertion (--bvals=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bval); bvals indicate the amount of the diffusion
%parameter telling that no smoothing should be applied (--fwhm=0)
%output from TOPUP algorithm (--topup=topup_AP_PA_b0)
%parameter that assumes a quadratic model for the EC-fields (--flm=quadratic)
%output name (--out=eddy_unwarped_images)
% system(cmd)

%% TENSOR FITTING - fast (for later TBSS) Can run for the first subject to look how it looks in FSLeyes
%fitting the tensor into the DTI data (this is to get fractional anisotropy (FA) to see whether there are microstructural changes or differences between groups, etc. this is usually used in connection with TBSS)

path5=dir([path_out '0*']);
for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    if ~isempty(path6)
        cd([path5(gg).folder '/' path5(gg).name]); %change directory for each subject, so that a nodif_PA.nii.gz is created for each subject
        cmd = ['dtifit --data=' path6(1).folder '/eddy_unwarped_images.nii.gz --mask=' path6(1).folder '/hifi_nodif_brain_mask.nii.gz --bvecs=' path6(2).folder '/' path6(2).name ' --bvals=' path6(1).folder '/' path6(1).name ' --out=' path6(1).folder '/dti_fitted_tensors'];
        system(cmd)
        disp(gg)
    end
end
% cmd = ['dtifit --data=eddy_unwarped_images.nii.gz --mask=hifi_nodif_brain_mask.nii.gz --bvecs=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bvec --bvals=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bval --out=dti_fitted_tensors'];
%dtifit command
%input data from eddy (--data=eddy_unwarped_images.nii.gz)
%mask that you got before (--mask=hifi_nodif_brain_mask.nii.gz)
%bvecs (same as for eddy)
%bvals (same as for eddy)
%output file (--out=dti_fitted_tensors)
% system(cmd)
%returns several outputs (starting with "dti_fitted_tensors..", NOTE in particular the ones ending with "_V1" "_V2" "_V3" and open them in fsleyes to conduct a further inspection
