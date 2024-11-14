function [participant_id,fn_out] = cp_prepLabelsRegr(pth_dat,pth_out)
% Function to prepare the BIDS labels and create the participants.tsv file
% with the 4 variables provided (age, sex, TIV, scanner). Moreover the 
% subjects order/label is randomized for improved anonymization.
% This function is called at the beginning of the main BIDS-ification
% function
% 
% FORMAT
%   fn_out = cp_prepLabelsRegr(pth_dat,pth_out)
% 
% INPUT
%   pth_dat : path to folder with all the data, see Readme
%   pth_out : path where to write the BIDSified data, see Readme
% 
% OUTPUT
%   participant_id : BIDS lables of the subjects
%   fn_out         : list of files created
% 
% EXAMPLE
% 
%   pth_dat = 'C:\Dox\2_Data\qMRI_MPM\Data4ChrisPhilips'
%   pth_out = 'C:\Dox\2_Data\qMRI_MPM\BIDS_AgingData'
%   [participant_id,fn_out] = cp_prepLabelsRegr(pth_dat,pth_out)
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

% List all generatred files
fn_out = char(...
    fn_IdKeys, ... % Key table for the BIDS-id and original-id match
    fn_participants_tsv,... % Participants file for BIDS
    fn_participants_json) ; % Description of particpants.tsv file

end