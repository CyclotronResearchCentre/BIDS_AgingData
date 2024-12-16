function fn_out = cp_prepTopJSON(pth_out)
% Function to prepare the other top-level JSON files, mainly the dataset
% description file plus the global JSON files describing the different
% types of image files.
% 
% FORMAT
%   fn_out = cp_prepTopJSON(pth_out)
% 
% INPUT
%   pth_out : path where to write the BIDSified data, see Readme
% 
% OUTPUT
%   fn_out         : list of files created
% 
% EXAMPLE
%   pth_out = 'C:\Dox\2_Data\qMRI_MPM\BIDS_AgingData'
%   fn_out = cp_prepTopJSON(pth_out)
% 
% TO DO!
% Add the description of
%   - warped quantitative maps
%   - modulate warped tissue maps
%   - tissue-weighted smmoothed warped quantitative maps
%   - ...
% 
%_______________________________________________________________________
% Copyright (C) 2024 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

% Write out the dataset_desription.json files
fn_dataset_desription_json = fullfile(pth_out,'dataset_description.json');
Tool = struct( ...
    'Name', {'SPM' , 'VBQ'}, ...
    'Version', {'8' , 'original'});
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
    'Licence', 'CC BY' ...
    );
spm_save(fn_dataset_desription_json, dataset_desription_json, ...
    'indent', '\t')

% TODO!
% Add the description tissue-weighted smmoothed warped quantitative maps: 
% A/PD, MTsat, R1 and R2* for the GM & WM.

fn_out = fn_dataset_desription_json;
end
