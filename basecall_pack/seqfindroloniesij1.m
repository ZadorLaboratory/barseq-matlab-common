function [rol, sig,rol0]=seqfindroloniesij1(filename,thresh)
%find rolonies and signal intensities in all channels in all tif images whose filenames containe filename.
%thresh is an n x 4 vector with intensity threshold values for finding rolonies.
%Uses the function "Find Maxima" in Fiji, need Fiji installed with Miji and
%need MATLAB running JAVA JRE 1.8

%this version added rol0, an optional output that records rolony positions
%as they were initially identified in ImageJ. This version calls rolony
%signals using only 


files=dir(fullfile(['*',filename,'*.tif']));
files=sort_nat({files.name});
rol=cell(length(files),1);
sig=rol;
rol0=rol;

% find rolonies using Miji and find intensities in matlab
Miji(false);
for i=1:length(files)
    im1=uint16([]);
    currol=[];
    %find rolonies
    if i<=size(thresh,1)
        currthresh=thresh(i,:);
    else
        warning(['threshold value missing for cycle ', num2str(i),'. Using last threshold value.']);
        currthresh=thresh(end,:);
    end
    
    for n=1:4
        im1(:,:,n)=imread(files{i},n);
        MIJ.createImage(im1(:,:,n));
        MIJ.run('Find Maxima...', ['noise=',num2str(currthresh(n)),' output=List']);
        currol=[currol;MIJ.getResultsTable];
        %cleanup
        MIJ.run('Clear Results');
        MIJ.closeAllWindows
    end
    rol0{i}=currol;
       
    %combine rolonies within 1.5 pixels of each other
    allpeaks=zeros(size(im1,1),size(im1,2));
    allpeaks(sub2ind([size(im1,1),size(im1,2)],currol(:,2)+1,currol(:,1)+1))=1;
    allpeaks=allpeaks&~imdilate(allpeaks,[1 1 1;1 0 0; 0 0 0]);
    [idx1,idx2]=ind2sub(size(allpeaks),find(allpeaks));
    rol{i}=[idx1,idx2];
    %figure;scatter(idx2,idx1);
    sig{i}=readintensities(rol{i},double(im1),[1 1 1;1 1 1;1 1 1]);%read intensities using the original rol0
    
end
MIJ.exit;
%save('rolonies.mat','rol','sig');
end