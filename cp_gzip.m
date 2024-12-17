function fn_out = cp_gzip(fn_in, flag)

% Function to gzip a list of files, either explicit of picked up with a
% regular expression filter, potentially recursively.
% 
% FORMAT
% 
%   fn_out = cp_gzip(fn_in, flag)
% 
% INPUT
% 
% fn_in     : single file, list of files, or path to a folder where files
%             will be selected to be gzipped
% flag      : structure with some options
%   .filt       : file filtering regular expression (Def. '^.*\.nii$')
%                 (see spm_select)
%   .rec        : act recursively (Def.) or not on folder
%   .delOrig    : delete the original file(s) or not (Def.)
% 
% OUPTUT
% 
% fn_out    : file(s) that were gzipped
% 
% EXAMPLE
% To gzip all nifti files in current folder and lower ones, not keeping the
% original ones
%   flag = struct(...
%     'filt','^.*\.nii$',... % all .nii files
%     'rec', true, ... % act recursively
%     'delOrig', true) % delete originals
%   pth_gzip = 'D:\ccc_DATA\qMRI_Ageing_MPM\BIDS_AgingData'
%   cp_gzip(pth_gzip, flag)
% 
% NOTE:
% When multiple files are passed, Gzip will put all the compressed files
% into the folder of the *1st* filer.
% Therefore, it is safer to loop over the files and ensure their
% compressed version is created next to the original one!
%_______________________________________________________________________
% Copyright (C) 2024 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

% Defaults & input
flag_ref = struct(...
    'filt','^.*\.nii$',... % all .nii files
    'rec', true, ... % act recursively
    'delOrig', false);  % do not delete original file after gzipping 
pth_ref = '.'; % use local location by default

if nargin<2, flag = flag_ref; end
if nargin<1, fn_in = pth_ref; end

% Find out if dealing with (list of) file(s) or a single path, from which
% to pick some files
if size(fn_in,1)>1 || exist(fn_in,'file')==2
    % assume it's single/list of file/s
    search_dir = false;
    fn_2gz = fn_in;
elseif size(fn_in,1)==1 && exist(fn_in,'dir')==7
    % assume single path to deal with
    search_dir = true;
    pth_seach = fn_in;
else
    % Unknown case
    beep
    fprintf('\nUnknow input format\n\t%s',fn_in')
end

% If searching a directory, then pick up all the target files
if search_dir
    % Do we go recursive or not ?
    if flag.rec
        listSelect = 'FPListRec';
    else
        listSelect = 'FPList';
    end
    fn_2gz = spm_select(listSelect,pth_seach,flag.filt);
end

% Gzipping things, one by one to keep thigns clean
for ii=1:size(fn_2gz,1)
    gzip(deblank(fn_2gz(ii,:)));
    % if requested delete original file
    if flag.delOrig
        delete(deblank(fn_2gz(ii,:)));
    end
end

% return list of gzipped files
fn_out = fn_2gz;

end

% pth_test = 'D:\ccc_DATA\P2021_MS7T_test\MS7T_Data\Test_gzip'
% fn_gz = spm_select('FPListRec',pth_test,'^.*\.nii$')
% gzip(cellstr(fn_gz))
% % 
% % When multiple files are passed, Gzip will put all the compressed files
% % into the folder of the *1st* file.
% % Therefore, it is safer to loop over the files and ensure their
% % compressed version is created next to the original one!
% tmp = fn_gz(3,:)
% fn_gz(3,:) = fn_gz(1,:);
% fn_gz(1,:) = tmp

