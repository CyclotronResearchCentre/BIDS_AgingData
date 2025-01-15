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
%       .ops  : flag the data to be dealt with, a 1x3 vector, def. [1 1 1]
%           * TW-smoothed warped quantiative maps
%           * warped quantitative and tissue maps
%           * subject space tissue maps
%       .mask : apply a whole-brain mask on the quantitative maps (1,
%               default) or not (0)
%
% OUTPUT
%   fn_out     : whole list of files in the BIDS folder
%   fn_out_nii : whole list of (gzipped) Nifti files
%
% EXAMPLE
%   pth_dat = 'D:\ccc_DATA\qMRI_Ageing_MPM\Data4ChrisPhilips'
%   pth_out = 'D:\ccc_DATA\qMRI_Ageing_MPM\BIDS_AgingData'
%   opt = struct( ...
%       'gzip', true, ... % -> gzip all .nii files at the end
%       'mask', true, ... % -> ICV-mask the qunatitive maps in MNI space
%       'ops', [1 1 0]);  % -> deal with warped and TW-smoothed data but
%                         %    not subject-space tissue maps
%   fn_out = cp_BIDSify_main(pth_dat,pth_out, opt)
%
% REFERENCE
% Callaghan et al. 2014, 
% https://doi.org/10.1016/j.neurobiolaging.2014.02.008
%
% PROCESS
% - start top level duties
%       1. randomizing the subjects list, produce 'participants.tsv' file
%       2. add the generic .json files
%       3. gather mean and mask images
% - then deal with the subjects images
%       1. the warped images, q-maps and modulated tissue class images
%          -> "SPM8_dartel" derivative
%       2. the tissue-weighted (GM & WM) smoothed q-maps
%          -> "VBQ_TWsmooth" derivative
%       3. the native space tissue class images
%          -> "SPM8_preproc" derivative
%
% STILL MISSING
% - data licence
% - full description of how the data were spatially processed before hand
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
    opt = struct( ...
        'gzip', false, ...
        'mask', true, ...
        'ops', [1 1 1]);
end

%% Deal with pathes for the BIDSified data, in separate derivative folders
if ~exist(pth_out,'dir'), mkdir(pth_out); end
% Main derivatieve folder
pth_deriv = fullfile(pth_out,'derivatives');
if ~exist(pth_deriv,'dir'), mkdir(pth_deriv); end
% Specific derivative folders: 'dartel', 'TWsmooth' and 'preproc'
pth_drv_TWsmooth = fullfile(pth_deriv,'VBQ_TWsmooth');
pth_drv_dartel = fullfile(pth_deriv,'SPM8_dartel');
pth_drv_preproc = fullfile(pth_deriv,'SPM8_preproc');

%% Deal with top level files
% 1. Labels and regressors -> participants.tsv file
%==================================================
% Original and BIDS labels needed to cath, copy & rename all data...
[participant_id,participant_orig] = cp_prepLabelsRegr(pth_dat,pth_out);
Nsubj = numel(participant_id);

% 2. Add the top-level .JSON files
%=================================
fn_dataset_desription_json = cp_prepTopJSON(pth_out); %#ok<*NASGU>

% 3. Arrange mean and mask images
%================================
fn_MaskMean = cp_prepMeanMask(pth_dat,pth_deriv);

% % 4. Create empty top level folder for all the subjects
% %======================================================
% for isub = 1:Nsubj
%     % Create subject's target folders
%     pth_isub_anat = fullfile(pth_out, ...
%         sprintf('sub-%s',participant_id{isub}),'anat');
%     if ~exist(pth_isub_anat,'dir'), mkdir(pth_isub_anat); end
%     fn_emptyT1w = fullfile( pth_isub_anat , ...
%         sprintf('sub-%s_T1w.nii',participant_id{isub}) );
%     fid = fopen(fn_emptyT1w,'wb'); 
%     fwrite(fid,0,'uint8'); 
%     fclose(fid);
% end

% 5. Create the .bidsignore files to pass validator
%==============================================
% It's a text file named '.bidsignore' with just this in body
%   /participants.tsv
%   /participants.json
%   *_T1w.nii


%% Deal with TW-smoothed individual subjects data
if opt.ops(1)
    if ~exist(pth_drv_TWsmooth,'dir'), mkdir(pth_drv_TWsmooth); end
    % 1. Define the path to all the images, 2 x 4 sets: [GM WM] x [A MTsat R1 R2*]
    %    + the different types of maps & tissues
    imgTypes_orig = {'A','MT','R1','R2s'};
    imgTypes = {'PDmap','MTsat','R1map','R2starmap'}; % BIDS suffixes
    tissueTypes = {'GM', 'WM'};
    pth_swqMRIs = cell(2,4);
    for ii=1:2 % tissue types
        for jj=1:4 % maps types
            pth_swqMRIs{ii,jj} = fullfile(pth_dat, ...
                sprintf('Fin_dart_p%d',ii), ...
                sprintf('Imgs_%s',imgTypes_orig{jj}));
        end
    end
    
    % 2. Deal with each subject one by one
    for isub = 1:Nsubj
        % Create subject's target folders
        pth_isub_anat = fullfile(pth_drv_TWsmooth, ...
            sprintf('sub-%s',participant_id{isub}),'anat');
        if ~exist(pth_isub_anat,'dir'), mkdir(pth_isub_anat); end
        
        % Deal with all 8 images
        for ii=1:2
            for jj=1:4
                % Source image full-filename
                fn_isub_orig = fullfile(pth_swqMRIs{ii,jj}, ...
                    sprintf('fin_dart_p%d%s_%s.nii', ...
                    ii,participant_orig{isub},imgTypes_orig{jj}) );
                % BIDS image full-filename
                fn_isub = fullfile(pth_isub_anat, ...
                    sprintf('sub-%s_space-MNI_desc-%ssmo_%s.nii', ...
                    participant_id{isub}, ...
                    tissueTypes{ii}, ...
                    imgTypes{jj} ) );
                if ~exist(fn_isub_orig,'file')
                    fprintf('\nERROR. Could not find file :\n\t%s\n', ...
                        fn_isub_orig);
                else
                    copyfile(fn_isub_orig,fn_isub)
                end
            end
        end
    end
end

%% Deal with warped individual subjects data: quantitative and tissue maps
if opt.mask
    % Pick up mask generated at top level
    fn_mask = fn_MaskMean(end-1,:);
    fl_imCalc.dtype = 16; % Force a float32 format as original data
end
if opt.ops(2)
    if ~exist(pth_drv_dartel,'dir'), mkdir(pth_drv_dartel); end
    % 1. Define the path to all the images:
    % - modulated warped tissue maps, mwc1/2/3
    % - warped quantitative maps, w*[A MTsat R1 R2*]
    % - warp images, u*MT
    % plus Dartel template #6
    pth_wMaps = fullfile(pth_dat, 'MPM_Processing');
    tissueTypes = {'GM', 'WM', 'CSF'};
    imgTypes_orig = {'A','MT','R1','R2s'};
    imgTypes = {'PDmap','MTsat','R1map','R2starmap'}; % BIDS suffixes
    
    % 2. Deal with each subject one by one
    for isub = 1:Nsubj
        % Create subject's target folders
        pth_isub_anat = fullfile(pth_drv_dartel, ...
            sprintf('sub-%s',participant_id{isub}),'anat');
        if ~exist(pth_isub_anat,'dir'), mkdir(pth_isub_anat); end
        
        % Deal with modulated warped tissue maps, mwc1/2/3
        for ii=1:3
            % Source image full-filename
            fn_isub_orig = fullfile(pth_wMaps, ...
                sprintf('mwc%d%s_MT.nii', ii,participant_orig{isub}) );
            % BIDS image full-filename
            fn_isub = fullfile(pth_isub_anat, ...
                sprintf('sub-%s_MTsat_space-MNI_desc-mod_label-%s_probseg.nii', ...
                participant_id{isub}, ...
                tissueTypes{ii} ) );
            if ~exist(fn_isub_orig,'file')
                fprintf('\nERROR. Could not find file :\n\t%s\n', ...
                    fn_isub_orig);
            else
                copyfile(fn_isub_orig,fn_isub)
            end
        end
        
        % Deal with warped quantitative maps, w*[A MTsat R1 R2*]
        for ii=1:4
            % Source image full-filename
            fn_isub_orig = fullfile(pth_wMaps, ...
                sprintf('w%s_%s.nii', participant_orig{isub}, imgTypes_orig{ii}) );
            % BIDS image full-filename
            fn_isub = fullfile(pth_isub_anat, ...
                sprintf('sub-%s_space-MNI_%s.nii', ...
                participant_id{isub}, ...
                imgTypes{ii} ) );
            if ~exist(fn_isub_orig,'file')
                fprintf('\nERROR. Could not find file :\n\t%s\n', ...
                    fn_isub_orig);
            else
                % Proceed with data copy, incl. masking
                if opt.mask
                    spm_imcalc(char(fn_isub_orig, fn_mask), fn_isub, ...
                        'i1.*i2', fl_imCalc);
                else
                    copyfile(fn_isub_orig,fn_isub)
                end
            end
        end
        
        % Deal with warp images, u*MT
        % Original image full-filename
        fn_isub_orig = fullfile(pth_wMaps, ...
            sprintf('u%s_MT.nii', participant_orig{isub}) );
        % BIDS image full-filename
        fn_isub = fullfile(pth_isub_anat, ...
            sprintf('sub-%s_MTsat_desc-dartelwarps.nii', ...
            participant_id{isub} ) );
        if ~exist(fn_isub_orig,'file')
            fprintf('\nERROR. Could not find file :\n\t%s\n', ...
                fn_isub_orig);
        else
            copyfile(fn_isub_orig,fn_isub)
        end
    end
    
    % 3. Deal with Dartel template #6
    fn_template_orig = spm_select('FPList',pth_wMaps,'^Template_6');
    fn_template = spm_file(fn_template_orig,'path',pth_drv_dartel);
    for ii=1:size(fn_template_orig,1)
        if ~exist(fn_template_orig(ii,:),'file')
            fprintf('\nERROR. Could not find file :\n\t%s\n', ...
                fn_template_orig(ii,:));
        else
            copyfile(fn_template_orig(ii,:),fn_template(ii,:))
        end
    end
end

%% Deal with the tissue maps still in subject space, c1/2/3*
if opt.ops(3)
    if ~exist(pth_drv_preproc,'dir'), mkdir(pth_drv_preproc); end
    % 1. Define the path to all the images, c1/2/3
    pth_cMaps = fullfile(pth_dat, 'MPM_Processing');
    tissueTypes = {'GM', 'WM', 'CSF'};
    
    % 2. Deal with each subject one by one
    for isub = 1:Nsubj
        % Create subject's target folders
        pth_isub_anat = fullfile(pth_drv_preproc, ...
            sprintf('sub-%s',participant_id{isub}),'anat');
        if ~exist(pth_isub_anat,'dir'), mkdir(pth_isub_anat); end
        
        % Deal with modulated warped tissue maps, mwc1/2/3
        for ii=1:3
            % Source image full-filename
            fn_isub_orig = fullfile(pth_wMaps, ...
                sprintf('c%d%s_MT.nii', ii,participant_orig{isub}) );
            % BIDS image full-filename
            fn_isub = fullfile(pth_isub_anat, ...
                sprintf('sub-%s_MTsat_space-SUBJ_label-%s_probseg.nii', ...
                participant_id{isub}, ...
                tissueTypes{ii} ) );
            if ~exist(fn_isub_orig,'file')
                fprintf('\nERROR. Could not find file :\n\t%s\n', ...
                    fn_isub_orig);
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
