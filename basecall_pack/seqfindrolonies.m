function [rol, sig]=seqfindrolonies(filename,thresh)
%find rolonies and signal intensities in all channels in all tif images whose filenames containe filename.
%thresh is an n x 1 vector with intensity threshold values for finding rolonies.

files=dir(fullfile(['*',filename,'*.tif']));
files=sort_nat({files.name});
rol=cell(length(files),1);
sig=rol;
for i=1:length(files)
    im1=[];
    for n=1:4
        im1(:,:,n)=imread(files{i},n);
    end
    if i<=length(thresh)
        rol{i}=findrolonies(im1,thresh(i));
    else
        warning(['threshold value missing for cycle ', num2str(i),'. Using last threshold value.']);
        rol{i}=findrolonies(im1,thresh(end));
    end
    sig{i}=readintensities(rol{i},im1,[1 1 1;1 1 1;1 1 1]);%read intensities using the original rol0
    
end

end