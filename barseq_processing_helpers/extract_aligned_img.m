function extract_aligned_img(folders)
%extract zipped aligned files from specified folders. if no folder is
%specified, then extrac tfrom all folders.
if ~exist('folders','var')
    folders=get_folders(); % if no folder is specified, then extract from all FOVs
elseif ~iscell(folders)
    folders={folders};
end

parfor i=1:numel(folders)
    zip_name=fullfile(folders{i},'aligned.zip');
    try
        unzip(zip_name)
    catch
    
    end

end

end