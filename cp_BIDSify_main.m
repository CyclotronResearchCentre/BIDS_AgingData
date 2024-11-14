function [fn_out, fn_out_nii] = cp_BIDSify_main(pth_dat,pth_out,opt)
% Function to BIDSify the processed aging data from Callaghan et al. 2014
% 
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
%   fn_out     : whole list of files in the BIDS folder
%   fn_out_nii : whole list of (gzipped) Nifti files
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
fn_MaskMean = cp_prepMeanMask(pth_dat,pth_deriv);

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
    flag_gz = struct(...
        'filt','^.*\.nii$',... % all .nii files
        'rec', true, ...       % act recursively
        'delOrig', true);      % delete original file after gzipping 
    fn_out_nii = cp_gzip(pth_out, flag_gz);
else
    fn_out_nii = spm_select('FPListRec',pth_out,'^.*\.nii$');
end

%% Collect output -> whole list of files in the BIDS folder
fn_out = spm_select('FPListRec',pth_out,'.*');
    
end
%%

% gzip(cellstr(fn_nii)) % -> puts all gzip file into top folder! :-(
