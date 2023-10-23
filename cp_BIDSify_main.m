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
%       * randomizing the subjects list, produce 'participants.tsv' file 
%       * gather mean and mask images
%       * add the generic .json files
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

% 2. Arrange mean and mask images
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

% deal with averaged maps images, only some MTsat maps available...

% 3. Add the top-level .JSON files
%=================================



%% Deal with individual subjects data
for i_sub = 1:Nsubj
    % Deal with each subject one by one
end

end
%%
% Writing .tsv file ?
% - create a table
% - use writetable -> writetable(T,'participants.txt','Delimiter','\t')
% - turn .txt into .tsv file -> movefile(SOURCE,DESTINATION)

