function [participant_id,participant_orig,fn_out] = cp_prepLabelsRegr(pth_dat,pth_out)
% Function to prepare the BIDS labels and create the participants.tsv file
% with the 4 variables provided (age, sex, TIV, scanner). Moreover the 
% subjects order/label is randomized for improved anonymization.
% This function is called at the beginning of the main BIDS-ification
% function.
% 
% NOTE:
% If the 'IdKeys' file already exist, then it is loaded and the SAME order
% of the subject is re-used.
% 
% FORMAT
%   [participant_id,participant_orig,fn_out] = cp_prepLabelsRegr(pth_dat,pth_out)
% 
% INPUT
%   pth_dat : path to folder with all the data, see Readme
%   pth_out : path where to write the BIDSified data, see Readme
% 
% OUTPUT
%   participant_id   : BIDS labels of the subjects
%   participant_orig : original labels of the subjects
%   fn_out           : list of files created
% 
% EXAMPLE
% 
%   pth_dat = 'D:\ccc_DATA\qMRI_Ageing_MPM\Data4ChrisPhilips'
%   pth_out = 'D:\ccc_DATA\qMRI_Ageing_MPM\BIDS_AgingData'
%   [label_id,label_orig,fn_out] = cp_prepLabelsRegr(pth_dat,pth_out)
%_______________________________________________________________________
% Copyright (C) 2024 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

% 1. Labels and regressors -> participants.tsv file
%==================================================
% Load labels and regressors
fn_labels_regressors = fullfile(pth_dat,'Subjects4Chris.mat');
load(fn_labels_regressors)
Nsubj = numel(Subjects4Chris.ID);
% Filename of table
fn_IdKeys = fullfile(pth_out,'IdKeys.tsv');

if exist(fn_IdKeys,'file')
    IdKeys = spm_load(fn_IdKeys);
    SubjPerm = IdKeys.Index_orig;
else
    % Randomization of subjects
    SubjPerm = randperm(Nsubj);
end

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
% write out the table with original id keys and permutation
participants_tb = table(participant_id,participant_orig,SubjPerm');
participants_tb.Properties.VariableNames{3} = 'Index_orig';
spm_save(fn_IdKeys,participants_tb)

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

% List all generatred files
fn_out = char(...
    fn_IdKeys, ... % Key table for the BIDS-id and original-id match
    fn_participants_tsv,... % Participants file for BIDS
    fn_participants_json) ; % Description of particpants.tsv file

end