

function organize_hyb_files(hybfolder)

% move files
    if ~exist('hybfolder','var')
        hybfolder='hyb01';
    end
    
    folders=dir('MAX*');
    folders=sort_nat({folders.name});
    
    
    cd(['../',hybfolder]);
    
    % copy sequential round images to the correct folders
    for i=1:length(folders)
        copyfile([folders{i},'.tif'],['../processed/',folders{i},'/','nuclear.tif']);
    end
    cd ../processed
end