function [seq,seqC,sig,scores]=pixelbasecall2d1(rol,foldername)
%locally align maxproj images  to first cycle (and rolony coordinates in
%first cycle). basecall using pixel values. This version uses tiny tiles
%for local alignment (50x50)

tilesize=150;
edgesize=30;

%% basecall first cycle: balance channels, read rolony signals, 
%read first cycle image.
olddir=cd(foldername);
files=dir('*.tif');
im=cell(length(files),1);
for i=1:length(files)
    im{i}=imread(files(i).name,1);
    for n=2:4
        im{i}(:,:,n)=imread(files(i).name,n);
    end
end
cd ..
labelim=zeros(size(im{1},1),size(im{1},2));
labelim(sub2ind([size(im{1},1),size(im{1},2)],rol(:,1),rol(:,2)))=1:size(rol,1);

%initialize signals

%initial call first base using imstack1 filtered by gaussian for
%calibration
sig=zeros(size(rol,1),length(files),4);%sequencing signals rol-cycle-channel
for i=1:4
    chim=imgaussfilt(im{1}(:,:,i),1);
    sig(labelim(labelim>0),1,i)=chim(labelim>0);
end
[maxsig,I]=max(sig(:,1,:),[],3);

%optional: balance bases using median of maxsig in each channel
for i=1:4
    sigmedian(i)=median(maxsig(I==i));
end
sigmedian=sigmedian/max(sigmedian);
for i=1:length(files)
    im{i}=double(im{i})./repmat(reshape(sigmedian,1,1,4),size(im{i},1),size(im{i},2),1);
end

%recall first base using the calibrated signals.
sig=zeros(size(rol,1),length(files),4);%sequencing signals rol-cycle-channel

for i=1:4
    chim=imgaussfilt(im{1}(:,:,i),1);
    sig(labelim(labelim>0),1,i)=chim(labelim>0);
end

%cut cycle 1 image into 400px by 400px tiles to prepare for alignment and
%basecalling of subsequent cycles. tiles overlap each other by 50%,so that
%every pixel is covered twice except at the edge. In each tile, only call
%rolonies at the center 300px by 300px range.



tile1=1:(tilesize-2*edgesize):floor(size(im{1},1)/(tilesize-2*edgesize)-1)*(tilesize-2*edgesize)+1;
tile2=1:(tilesize-2*edgesize):floor(size(im{1},2)/(tilesize-2*edgesize)-1)*(tilesize-2*edgesize)+1;
tile1=repmat(tile1,length(tile2),1);
tile2=repmat(tile2',1,length(tile1));
tile=reshape(cat(3,tile1,tile2),[],2);%this has the starting index of each tile.


%% for subsequent cycles:
for n=2:length(files)
    %cut images into 500px x 500 px images, and align each small stack with cycle 1.
    for i=1:size(tile,1)
        
        if tile(i,1)+tilesize*2-2>size(im{1},1)
            e1=size(im{1},1);
        else
            e1=tile(i,1)+tilesize-1;
        end
        
        if tile(i,2)+tilesize*2-2>size(im{1},2)
            e2=size(im{1},2);
        else
            e2=tile(i,2)+tilesize-1;
        end
        
        im1tile=im{1}(tile(i,1):e1,tile(i,2):e2,:);
        im2tile=im{n}(tile(i,1):e1,tile(i,2):e2,:);
        labelimtile=labelim(tile(i,1):e1,tile(i,2):e2);
        labelimtile(1:edgesize,:)=0;
        labelimtile(:,1:edgesize)=0;
        labelimtile(end-edgesize+1:end,:)=0;
        labelimtile(:,end-edgesize+1:end)=0;
        
        
        %align tiles in xy plane allowing rotation and translation.
        %Then align a 3x3 area around each rolony along z using xcorr
        
        %make z-proj, calculate tform, and apply to the whole stack.
        %optimizer = registration.optimizer.OnePlusOneEvolutionary;
        %metric = registration.metric.MattesMutualInformation;
        %tform1 = imregtform(sum(im2tile,3),sum(im1tile,3),'rigid',optimizer,metric);%doesn't work for boundary slices due to the bounary being very bright. Need to detect boundary and exclude from alignment
        tform1 = imregcorr(sum(im2tile,3),sum(im1tile,3),'rigid','Window',0);%doesn't work for boundary slices due to the bounary being very bright. Need to detect boundary and exclude from alignment
        
        Rfixed=imref2d(size(im1tile(:,:,1)));
        warpedim2tile=imwarp(im2tile,tform1,'OutputView',Rfixed);%channel-summed stack warped in xy, used for z-alignment
        %call base 
        for k=1:4
            chim=imgaussfilt(warpedim2tile(:,:,k),1);
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
seq2d=seq;seqC2d=seqC;scores2d=scores;sig2d=sig;

%% save
save('seq2d.mat','seq2d','seqC2d','scores2d','sig2d');



