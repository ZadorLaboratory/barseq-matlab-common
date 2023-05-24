function [seq,seqC,sig,scores]=pixelbasecall(rol,foldername)
%locally align image stacks in 3d of each cycle  to first cycle (and rolony coordinates in
%first cycle). basecall using pixel values in the correct z-plane.

%need to fix registration for edge tiles.



%% basecall first cycle: balance channels, read rolony signals, 
%read first cycle image.
alignedfolders=dir(['*',foldername,'*']);

olddir=cd(alignedfolders(1).name);
stackfiles=dir('*.tif');
tcz=zeros(length(stackfiles),3);
prefixend=regexp(stackfiles(1).name,'T[01234567890]+C[01234567890]+Z');
for k=1:length(stackfiles)
    tcz(k,:)=cell2mat(textscan(stackfiles(k).name(prefixend:end),'T%uC%uZ%u')); %change this if file names change
end
stackfiles(tcz(:,2)>4)=[];%remove non-sequencing channels
tcz(tcz(:,2)>4,:)=[];%remove non-sequencing channels

imstack1=cell(1,1,length(stackfiles));
for k=1:length(stackfiles)
    imstack1{k}=imread(stackfiles(k).name);
end
imstack1=cell2mat(imstack1);

cd(olddir);

%make 3d rolony label image 
labelim=zeros(size(imstack1,1),size(imstack1,2),max(tcz(:,3)));
labelim(sub2ind([size(imstack1,1),size(imstack1,2),max(tcz(:,3))],rol(:,1),rol(:,2),rol(:,3)))=1:size(rol,1);

%initialize signals

%initial call first base using imstack1 filtered by gaussian for
%calibration
sig=zeros(size(rol,1),length(alignedfolders),4);%sequencing signals rol-cycle-channel
for i=1:4
    chim=imgaussfilt(imstack1(:,:,tcz(:,2)==i),1);
    sig(labelim(labelim>0),1,i)=chim(labelim>0);
end
[maxsig,I]=max(sig(:,1,:),[],3);

%optional: balance bases using median of maxsig in each channel
for i=1:4
    sigmedian(i)=median(maxsig(I==i));
end
sigmedian=sigmedian/max(sigmedian);

for i=1:4
    imstack1(:,:,tcz(:,2)==i)=uint16(imstack1(:,:,tcz(:,2)==i)/sigmedian(i));
end
%recall first base using the calibrated signals.
sig=zeros(size(rol,1),length(alignedfolders),4);%sequencing signals rol-cycle-channel

for i=1:4
    chim=imgaussfilt(imstack1(:,:,tcz(:,2)==i),1);
    sig(labelim(labelim>0),1,i)=chim(labelim>0);
end

%make a channel-sum z-stack for alignment later
chsum1=zeros(size(imstack1,1),size(imstack1,2),max(tcz(:,3)));
for i=1:max(tcz(:,3))
    chsum1(:,:,i)=sum(imstack1(:,:,tcz(:,3)==i),3);
end

%z-profile for each rolony location.
labelimmaxp=max(labelim,[],3);
rolzprofile=zeros(size(rol,1),max(tcz(:,3)));
for i=1:max(tcz(:,3))
    chsum1z=chsum1(:,:,i);
    rolzprofile(labelimmaxp(labelimmaxp>0),i)=imgaussfilt(chsum1z(labelimmaxp>0),1);
end

%remove image stack to save memory




%cut cycle 1 image into 400px by 400px tiles to prepare for alignment and
%basecalling of subsequent cycles. tiles overlap each other by 50%,so that
%every pixel is covered twice except at the edge. In each tile, only call
%rolonies at the center 300px by 300px range.



tile1=1:300:floor(size(imstack1,1)/300-1)*300+1;
tile2=1:300:floor(size(imstack1,2)/300-1)*300+1;
tile1=repmat(tile1,length(tile2),1);
tile2=repmat(tile2',1,length(tile1));
tile=reshape(cat(3,tile1,tile2),[],2);%this has the starting index of each tile.

clearvars imstack1 chim chsum1z ;


%% for subsequent cycles:
for n=2:length(alignedfolders)
    %read cycle n
    
    olddir=cd(alignedfolders(n).name);
    stackfiles2=dir('*.tif');
    tcz2=zeros(length(stackfiles2),3);
    prefixend=regexp(stackfiles2(1).name,'T[01234567890]+C[01234567890]+Z');
    for k=1:length(stackfiles2)
        tcz2(k,:)=cell2mat(textscan(stackfiles2(k).name(prefixend:end),'T%uC%uZ%u')); %change this if file names change
    end
    stackfiles2(tcz2(:,2)>4)=[];%remove non-sequencing channels
    tcz2(tcz2(:,2)>4,:)=[];%remove non-sequencing channels
    
    imstack2=cell(1,1,length(stackfiles2));
    for k=1:length(stackfiles2)
        imstack2{k}=imread(stackfiles2(k).name);
    end
    imstack2=cell2mat(imstack2);
    
    
    
    
    %scale channels by calibration obtained from cycle 1.
    for i=1:4
        imstack2(:,:,tcz2(:,2)==i)=uint16(imstack2(:,:,tcz2(:,2)==i)/sigmedian(i));
    end
    
    %make a channel-sum z-stack for alignment later
    chsum2=zeros(size(imstack2,1),size(imstack2,2),max(tcz2(:,3)));
    for i=1:max(tcz2(:,3))
        chsum2(:,:,i)=sum(imstack2(:,:,tcz2(:,3)==i),3);
    end
    cd(olddir);
    
    %cut images into 500px x 500 px images, and align each small stack with cycle 1.
    for i=1:size(tile,1)
        
        if tile(i,1)+798>size(chsum1,1)
            e1=size(chsum1,1);
        else
            e1=tile(i,1)+399;
        end
        
        if tile(i,2)+798>size(chsum1,2)
            e2=size(chsum1,2);
        else
            e2=tile(i,2)+399;
        end
        
        chsum1tile=chsum1(tile(i,1):e1,tile(i,2):e2,:);
        chsum2tile=chsum2(tile(i,1):e1,tile(i,2):e2,:);
        imstack2tile=imstack2(tile(i,1):e1,tile(i,2):e2,:);
        labelimtile=labelim(tile(i,1):e1,tile(i,2):e2,:);
        labelimtile(1:50,:,:)=0;
        labelimtile(:,1:50,:)=0;
        labelimtile(end-49:end,:,:)=0;
        labelimtile(:,end-49:end,:)=0;
        
        
        %align tiles in xy plane allowing rotation and translation.
        %Then align a 3x3 area around each rolony along z using xcorr
        
        %make z-proj, calculate tform, and apply to the whole stack.
        chsum1tilemax=max(chsum1tile,[],3);
        chsum2tilemax=max(chsum2tile,[],3);
        
        optimizer = registration.optimizer.OnePlusOneEvolutionary;
        metric = registration.metric.MattesMutualInformation;
        tform1 = imregtform(chsum2tilemax,chsum1tilemax,'rigid',optimizer,metric);%doesn't work for boundary slices due to the bounary being very bright. Need to detect boundary and exclude from alignment
        Rfixed=imref2d(size(chsum1tilemax));
        warpedchsum2tile=imwarp(chsum2tile,tform1,'OutputView',Rfixed);%channel-summed stack warped in xy, used for z-alignment
        imstack2tile=imwarp(imstack2tile,tform1,'OutputView',Rfixed);%original stack warped in xy, used for basecalling
        %for each rolony, align in the z-direction by maximizing cross
        %correlation
        labelimtilemaxp=max(labelimtile,[],3);
        tilerollist=labelimtilemaxp(labelimtilemaxp>0);%in the sequence of linear indexing in a single plane
        %z-profile for each rolony location in this cycle.
        rolzprofiletile=zeros(length(tilerollist),max(tcz(:,3)));
        for k=1:max(tcz(:,3))
            chsum2z=warpedchsum2tile(:,:,k);
            rolzprofiletile(:,k)=imgaussfilt(chsum2z(labelimtilemaxp>0),1); %note that rolzprofiletile follows the sequence of linear indexing of labelimtile, not rolony index as rolzprofile
        end
        
        %shift labelimtile in z to maximize xcorr
        for k=1:length(tilerollist)
            %tic
            [c,l] = xcorr(rolzprofiletile(k,:),rolzprofile(tilerollist(k),:));
            [~,I] = max(c);
            t = l(I);
            [idx1,idx2]=ind2sub(size(labelimtilemaxp),find(labelimtilemaxp==tilerollist(k)));
            labelimtile(idx1,idx2,:)=circshift(labelimtile(idx1,idx2,:),t);
            %toc
        end
        %call base using imstack2
        for k=1:4
            chim=imstack2tile(:,:,tcz2(:,2)==k);
            sig(labelimtile(labelimtile>0),n,k)=chim(labelimtile>0);
        end
       
    end
end

%% read out sequence and scores
%sequence
[~,seq]=max(sig,[],3);
%sequence in ATCG
seqC(seq==1)='G';
seqC(seq==2)='T';
seqC(seq==3)='A';
seqC(seq==4)='C';
seqC=reshape(seqC,length(rol),[]);
%scores
scores=max(sig,[],3)./sqrt(sum(sig.^2,3));


%% save
save('seq.mat','seq','seqC','scores','sig');



