function fn_out = cp_BIDSify_main(pth_dat,pth_out)
% Function to BIDSify the processed aging data from Callaghan et al. 2014
% FORMAT
%   fn_out = cp_BIDSify_main(pth_raw,pth_out)
% 
% INPUT
%   pth_dat : path to folder with all the data, see Readme
%   pth_out : path where to write the BIDSified data, see Readme
% 
% OUTPUT
%   fn_out : structure with full path file names
% 
% EXAMPLE
%   pth_dat = 'C:\Dox\2_Data\qMRI_MPM\Data4ChrisPhilips'
%   pth_out = 'C:\Dox\2_Data\qMRI_MPM\BIDS_AgingData'
%   fn_out = cp_BIDSify_main(pth_raw,pth_out)
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
%_______________________________________________________________________
% Copyright (C) 2023 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% Input check and defautl values
if nargin==0
    pth_dat = 'C:\Dox\2_Data\qMRI_MPM\Data4ChrisPhilips';
end
if nargin<2
    pth_out = pth_dat;
end
if ~exist(pth_out,'dir'), mkdir(pth_out); end
pth_deriv = fullfile(pth_out,'derivatives','SPM12dartel');
if ~exist(pth_deriv,'dir'), mkdir(pth_deriv); end

%% Deal with top level files
% 1. Labels and regressors -> participants.tsv file
%==================================================
% Load labels and regressors
fn_labels_regressors = fullfile(pth_dat,'Subjects4Chris.mat');
load(fn_labels_regressors)
Nsubj = numel(Subjects4Chris.ID);

% Randomization of subjects
SubjPerm = randperm(Nsubj);

% Create all arrays with participants informations
participant_id = cell(Nsubj,1);
participant_orig = cell(Nsubj,1);
age = cell(Nsubj,1);
sex = cell(Nsubj,1);
TIV = cell(Nsubj,1);
scanner = cell(Nsubj,1);
for isub = 1:Nsubj
    participant_id{isub} = sprintf('S%03d',isub);
    participant_orig{isub} = Subjects4Chris.ID{SubjPerm(isub)} ;
    age{isub} = Subjects4Chris.Age(SubjPerm(isub));
    if Subjects4Chris.Gender(SubjPerm(isub))
        sex{isub} = 'M';
    else
        sex{isub} = 'F';
    end
    TIV{isub} = Subjects4Chris.TIV{SubjPerm(isub)};
    if Subjects4Chris.Trio(SubjPerm(isub))
        scanner{isub} = 'trio';
    else
        scanner{isub} = 'quatro';
    end
end
% write out the table with original id keys
fn_IdKeys = fullfile(pth_out,'IdKeys.tsv');
spm_save(fn_IdKeys,table(participant_id,participant_orig))

% write out the participants.tsv/.json files
fn_participants_tsv = fullfile(pth_out,'participants.tsv');
spm_save(fn_participants_tsv,table(participant_id,age,sex,TIV,scanner))

fn_participants_json = spm_file(fn_participants_tsv,'ext','json');
participants_json = struct( ...
    'age' , struct( ...
        'Description', 'age of the participant', ...
        'Units', 'years'), ...
    'sex', struct( ...
        'Description', 'sex of the participant', ...
        'Levels', struct('M', 'male', 'F', 'female') ), ...
    'TIV', struct( ...
        'Description', 'total intracranial volume', ...
        'Units', 'liters'), ...
    'scanner', struct( ...
        'Description', 'acquisition scanner', ...
        'Levels', struct('quatro', 'quatro machine', 'trio', 'trio machine') ) ...
    );
spm_save(fn_participants_json, participants_json, 'indent', '\t')

% 2. Add the top-level .JSON files
%=================================
% write out the dataset_desription.json files
fn_dataset_desription_json = fullfile(pth_out,'dataset_desription.json');
Tool = struct( ...
    'Name', {'SPM' , 'VBQ'}, ...
    'Version', {'12' , 'original'});
GenBy = struct('GeneratedBy',Tool);
dataset_desription_json = struct( ...
    'name' , 'Processed MPM qMRI aging data', ...
    'BIDSversion', 'v1.8.0 + BEP011 + BEP038', ...
    'DatasetType', 'derivative', ...
    'Authors', 'Martina F. Callaghan & Christophe Phillips', ...
    'Acknowledgements', ['Elaine Anderson, Marinella Cappelletti, Rumana ', ...
        'Chowdhury, Joern Diedirchsen, Thomas H. B. Fitzgerald, and Peter ',...
        'Smittenaar as part of multiple cognitive neuroimaging studies ',...
        'performed at the Wellcome Centre for Human Neuroimaging'], ...
    'HowToAcknowledge', ['Please cite this paper: ', ...
        'https://doi.org/10.1016/j.neurobiolaging.2014.02.008'], ...
    'ReferencesAndLinks', 'https://doi.org/10.1016/j.neurobiolaging.2014.02.008', ...
    'GeneratedBy', GenBy, ...
    'Licence', '???' ...
    );
spm_save(fn_dataset_desription_json, dataset_desription_json, 'indent', '\t')

% TODO!
% Add the description tissue-weighted smmoothed warped quantitative maps: 
% A/PD, MTsat, R1 and R2* for the GM & WM.


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
imgTypes = {'A','MT','R1','R2s'};
pth_qMRIs = cell(2,4);
for ii=1:2
    for jj=1:4
        pth_qMRIs{ii,jj} = fullfile(pth_dat, ...
            sprintf('Fin_dart_p%d',ii),sprintf('Imgs_%s',imgTypes{jj}));
    end
end

for i_sub = 1:Nsubj
    % Deal with each subject one by one
    
end

end
%%

