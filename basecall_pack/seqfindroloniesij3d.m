function rol=seqfindroloniesij3d(maxprojfilename,stackname,thresh)
%find rolonies and signal intensities in all channels in all tif images whose filenames containe filename.
%thresh is an n x 4 vector with intensity threshold values for finding rolonies.
%Uses the function "Find Maxima" in Fiji, need Fiji installed with Miji and
%need MATLAB running JAVA JRE 1.8
%This version finds rolonies in the maxproj image, but additionally locate the
%rolonies in z-stack. run this in the maxproj folder. This version also
%does not read out rolony intensities.

files=dir(['*',maxprojfilename,'*.tif']);
filesstack=dir(['../*',stackname,'*']);
if length(files)~=length(filesstack)
    error('unequal number of stack and maxproj files.');
end

rol=cell(length(files),1);

% find rolonies using Miji and find intensities in matlab
Miji(false);
for i=1:length(files)
    allrol=[];
    chidx=[];
    sig=[];
    %find rolonies,the channel the rolonies are identified in, and the
    %signals
    if i<=size(thresh,1)
        currthresh=thresh(i,:);
    else
        warning(['threshold value missing for cycle ', num2str(i),'. Using last threshold value.']);
        currthresh=thresh(end,:);
    end
    
    for n=1:4
        im1=uint16(imread(files(i).name,n));
        MIJ.createImage(im1);
        MIJ.run('Find Maxima...', ['noise=',num2str(currthresh(n)),' output=List']);
        currol=MIJ.getResultsTable;
        allrol=[allrol;currol];
        chidx=[chidx;n*ones(size(currol,1),1)];
        sig=[sig;im1(sub2ind([size(im1,1),size(im1,2)],currol(:,2)+1,currol(:,1)+1))];
        
        %cleanup
        MIJ.run('Clear Results');
        MIJ.closeAllWindows
    end
       
    %combine rolonies within 1.5 pixels of each other
    labelim=zeros(size(im1,1),size(im1,2));
    labelim(sub2ind([size(im1,1),size(im1,2)],allrol(:,2)+1,allrol(:,1)+1))=1:length(allrol);
    
    sigim=labelim;
    sigim(sub2ind([size(im1,1),size(im1,2)],allrol(:,2)+1,allrol(:,1)+1))=sig;
    
    %remove rolonies whose signals are weaker than another rolony within
    %1.5 pixels distance
    filt=zeros(3,3,4);
    filt([1 2 3 4 6 7 8 9])=1;
    for n=1:length(filt)
        sigim1=conv2(sigim,filt(:,:,1),'same');
        labelim(sigim<sigim1)=0;
    end
    rol{i}=allrol(unique(labelim(labelim>0)),[2 1])+1;
    chidx=chidx(unique(labelim(labelim>0)));
    sig=sig(unique(labelim(labelim>0)));
    
    %read stack image, assign z position for each rolony by finding z with max signal in the channel the rolony was found in.
    olddir=cd(['../',filesstack(i).name]);
    stackfiles=dir('*.tif');
    tcz=zeros(length(stackfiles),3);
    prefixend=regexp(stackfiles(1).name,'T[01234567890]+C[01234567890]+Z');
    for k=1:length(stackfiles)
        tcz(k,:)=cell2mat(textscan(stackfiles(k).name(prefixend:end),'T%uC%uZ%u')); %change this if file names change
    end
    stackfiles(tcz(:,2)>4)=[];%remove non-sequencing channels
    tcz(tcz(:,2)>4,:)=[];%remove non-sequencing channels
    
    imstack=cell(1,1,length(stackfiles));
    for k=1:length(stackfiles)
        imstack{k}=imread(stackfiles(k).name);
    end
    imstack=cell2mat(imstack);
    
    rol{i}(:,3)=0;
    for k=1:size(rol{i},1)
        [~,I]=max(imstack(rol{i}(k,1),rol{i}(k,2),tcz(:,2)==chidx(k)));
        z=tcz(tcz(:,2)==chidx(k),3);
        rol{i}(k,3)=z(I);
    end
    
    cd(olddir);
    
    %[rol{i}, sig{i}=readintensitiesed(rol{i},imstack,[1 1 1;1 1 1;1 1 1]);%read intensities, not working
    
end
MIJ.exit;
save('rolonies.mat','rol');
end