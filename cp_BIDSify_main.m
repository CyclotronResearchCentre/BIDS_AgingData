function fn_out = cp_BIDSify_main(pth_dat,pth_out,opt)
% Function to BIDSify the processed aging data from Callaghan et al. 2014
% FORMAT
%   fn_out = cp_BIDSify_main(pth_raw,pth_out)
% 
% INPUT
%   pth_dat : path to folder with all the data, see Readme
%   pth_out : path where to write the BIDSified data, see Readme
%   opt     : option structure flag
%       .gzip : zip all the NIfTI files (1) or not (0, default)
% 
% OUTPUT
%   fn_out : whole list of files in the BIDS folder
% 
% EXAMPLE
%   pth_dat = 'C:\Dox\2_Data\qMRI_MPM\Data4ChrisPhilips'
%   pth_out = 'C:\Dox\2_Data\qMRI_MPM\BIDS_AgingData'
%   opt = struct('gzip', true); % -> gzip all .nii files at the end
%   fn_out = cp_BIDSify_main(pth_dat,pth_out, opt)
% 
% REFERENCE
% Callaghan et al. 2014, https://doi.org/10.1016/j.neurobiolaging.2014.02.008
% 
% PROCESS
% - start top level duties
%       1. randomizing the subjects list, produce 'participants.tsv' file 
%       2. add the generic .json files
%       3. gather mean and mask images
% - then deal with the subjects images
% 
% STILL MISSING
% - data licence 
% - full description of how the data were spatially processed
% - JSON file describing the tissue-weighted smoothed normalized 
%   quantitative maps, globally for all the subjects.
%_______________________________________________________________________
% Copyright (C) 2023 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% Input check and default values
if nargin==0
    pth_dat = pwd;
end
if nargin<2
    pth_out = pth_dat;
end
if nargin<3
    opt = struct('gzip', false);
end

% Deal with pathes for the BIDSified data
if ~exist(pth_out,'dir'), mkdir(pth_out); end
pth_deriv = fullfile(pth_out,'derivatives','SPM12dartel');
if ~exist(pth_deriv,'dir'), mkdir(pth_deriv); end

%% Deal with top level files
% 1. Labels and regressors -> participants.tsv file
%==================================================
participant_id = cp_prepLabelsRegr(pth_dat,pth_out);
Nsubj = numel(participant_id);

% 2. Add the top-level .JSON files
%=================================
fn_dataset_desription_json = cp_prepTopJSON(pth_out); %#ok<*NASGU>

% 3. Arrange mean and mask images
%================================
% Deal with mask images, in the top data folder, and name them
% - atlas-GM_space-MNI_mask.nii/.json
% - atlas-WM_space-MNI_mask.nii/.json
fn_GMmask_orig = fullfile(pth_dat,'WinnerTakesAllMask_GM.nii');
fn_GMmask = fullfile(pth_deriv,'atlas-GM_space-MNI_mask.nii');
copyfile(fn_GMmask_orig,fn_GMmask)

fn_WMmask_orig = fullfile(pth_dat,'WinnerTakesAllMask_WM.nii');
fn_WMmask = fullfile(pth_deriv,'atlas-WM_space-MNI_mask.nii');
copyfile(fn_WMmask_orig,fn_WMmask)

fn_GMmask_json = spm_file(fn_GMmask,'ext','json');
GMmask_json = struct( ...
    'Name', 'Grey matter mask', ...
    'Description', ['Voxel with higher probability, on average over the ' ...
        'subjects in the study, to be grey matter than white matter ' ...
        'or CSF and with a priori probability higher than .2 to be GM.'], ...
    'DerivedFrom', 'Segmented GM, WM and CSF maps from all subjects', ...
    'ReferencesAndLinks', 'https://doi.org/10.1016%2Fj.neurobiolaging.2014.02.008');
spm_save(fn_GMmask_json, GMmask_json, 'indent', '\t')

fn_WMmask_json = spm_file(fn_WMmask,'ext','json');
WMmask_json = struct( ...
    'Name', 'White matter mask', ...
    'Description', ['Voxel with higher probability, on average over the ' ...
        'subjects in the study, to be white matter than frey matter ' ...
        'or CSF and with a priori probability higher than .2 to be GM.'], ...
    'DerivedFrom', 'Segmented GM, WM and CSF maps from all subjects', ...
    'ReferencesAndLinks', 'https://doi.org/10.1016%2Fj.neurobiolaging.2014.02.008');
spm_save(fn_WMmask_json, WMmask_json, 'indent', '\t')

% deal with averaged map image(s), only some mean MTsat maps available but
% with 3 different masks applied. :-(
fn_meanMTsat_orig1 = fullfile(pth_dat,'Average_MTsat_Maps','AveragedMTMap_DARTEL.nii');
fn_meanMTsat_orig2 = fullfile(pth_dat,'Average_MTsat_Maps','AveragedMT_MNI_BrainMask.nii');
fn_meanMTsat_orig3 = fullfile(pth_dat,'Average_MTsat_Maps','Averaged_MT_HeadMask.nii');

% orig1 = DARTEL mean?           -> atlas-MTsat_space-MNI_desc-meanFull.nii/json, lower resolution
fn_meanMTsat_1 = fullfile(pth_deriv,'atlas-MTsat_space-MNI_res-low_desc-meanFull.nii');
copyfile(fn_meanMTsat_orig1,fn_meanMTsat_1)
% orig2 = mean with tight mask   -> atlas-MTsat_space-MNI_desc-meanICV.nii/json, high resolution
fn_meanMTsat_2 = fullfile(pth_deriv,'atlas-MTsat_space-MNI_res-high_desc-meanICV.nii');
copyfile(fn_meanMTsat_orig2,fn_meanMTsat_2)
% orig3 = mean with broader mask -> atlas-MTsat_space-MNI_desc-mean.nii/json, high resolution
fn_meanMTsat_3 = fullfile(pth_deriv,'atlas-MTsat_space-MNI_res-high_desc-mean.nii');
copyfile(fn_meanMTsat_orig3,fn_meanMTsat_3)

% deal with JSONs
fn_meanMTsat_1_json = spm_file(fn_meanMTsat_1,'ext','json');
meanMTsat_1_json = struct( ...
    'Name', 'mean MTsat map', ...
    'Description', ['Mean image, over all the subjects, of the warped MTsat ' ...
        'maps obtained directly from DARTEL without any masking.'], ...
    'DerivedFrom', 'MTsat images from all subjects', ...
    'Resolution', 'Matching the Dartel template: [1.5 1.5 1.5] mm^3', ...
    'ReferencesAndLinks', 'https://doi.org/10.1016%2Fj.neurobiolaging.2014.02.008');
spm_save(fn_meanMTsat_1_json, meanMTsat_1_json, 'indent', '\t')

fn_meanMTsat_2_json = spm_file(fn_meanMTsat_2,'ext','json');
meanMTsat_2_json = struct( ...
    'Name', 'mean MTsat map ICV masked', ...
    'Description', ['Mean image, over all the subjects, of the warped MTsat ' ...
        'masked with SPM''s "intra-cranial volume" mask.'], ...
    'DerivedFrom', 'MTsat images from all subjects', ...
    'Resolution', 'Matching the resulting image: [1 1 1] mm^3', ...
    'ReferencesAndLinks', 'https://doi.org/10.1016%2Fj.neurobiolaging.2014.02.008');
spm_save(fn_meanMTsat_2_json, meanMTsat_2_json, 'indent', '\t')

fn_meanMTsat_3_json = spm_file(fn_meanMTsat_3,'ext','json');
meanMTsat_3_json = struct( ...
    'Name', 'mean MTsat map skull masked', ...
    'Description', ['Mean image, over all the subjects, of the warped MTsat ' ...
        'masked, including bits of the skull.'], ...
    'DerivedFrom', 'MTsat images from all subjects', ...
    'Resolution', 'Matching the resulting image: [1 1 1] mm^3', ...
    'ReferencesAndLinks', 'https://doi.org/10.1016%2Fj.neurobiolaging.2014.02.008');
spm_save(fn_meanMTsat_3_json, meanMTsat_3_json, 'indent', '\t')

%% Deal with individual subjects data
% 1. Define the path to all the images, 2 x 4 sets: [GM WM] x [A MTsat R1 R2*]
%    + the different types of maps & tissues
imgTypes_orig = {'A','MT','R1','R2s'};
imgTypes = {'PDmap','MTsat','R1map','R2starmap'}; % BIDS suffixes
tissueTypes = {'GM', 'WM'};
pth_qMRIs = cell(2,4);
for ii=1:2 % tissue types
    for jj=1:4 % maps types
        pth_qMRIs{ii,jj} = fullfile(pth_dat, ...
            sprintf('Fin_dart_p%d',ii),sprintf('Imgs_%s',imgTypes_orig{jj}));
    end
end

% 2. Deal with each subject one by one
for isub = 1:Nsubj
    % Create subject's folders
    pth_isub_anat = fullfile(pth_deriv,sprintf('sub-%s',participant_id{isub}),'anat');
    if ~exist(pth_isub_anat,'dir'), mkdir(pth_isub_anat); end
    
    % Deal with all 8 images
    for ii=1:2
        for jj=1:4
            fn_isub_orig = fullfile(pth_qMRIs{ii,jj}, ...
                sprintf('fin_dart_p%d%s_%s.nii',ii,participant_orig{isub},imgTypes_orig{jj}) );
            fn_isub = fullfile(pth_isub_anat, ...
                sprintf('sub-%s_space-MNI_desc-%ssmo_%s.nii', ...
                    participant_id{isub}, ...
                    tissueTypes{ii}, ...
                    imgTypes{jj} ) );
            if ~exist(fn_isub_orig,'file')
                fprintf('\nERROR. Could not find file :\n\t%s\n', fn_isub_orig);
            else
                copyfile(fn_isub_orig,fn_isub)
            end
        end
    end
end

%% GZIP all the .nii files to save some space
if opt.gzip
    fn_nii = spm_select('FPListRec',pth_out,'^.*\.nii$');
    for ii=1:size(fn_nii,1)
        gzip(deblank(fn_nii(ii,:))); % Gzip in situ
        delete(deblank(fn_nii(ii,:))); % Delete orginal .nii file
    end
end

%% Collect output -> whole list of files in the BIDS folder
fn_out = spm_select('FPListRec',pth_out,'.*');
    
end
%%

% gzip(cellstr(fn_nii)) % -> puts all gzip file into top folder! :-(
