function segim=mmsegmentbc2(currrol1,varargin)

%% parse inputs

p=inputParser;
addParameter(p,'bgn',50,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'lthresh',800,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'projthresh',0.2,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'rejthresh',0.8,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'cellmasksize',30,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'cellratio',0.5,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'maxcellrange',300,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'mincellarea',400,@(x)isnumeric(x)&&isscalar(x));

parse(p,varargin{:});
bgn  = p.Results.bgn;
lthresh = p.Results.lthresh;
projthresh  = p.Results.projthresh;
rejthresh  = p.Results.rejthresh;
cellmasksize  = p.Results.cellmasksize;
cellratio  = p.Results.cellratio;
maxcellrange  = p.Results.maxcellrange;
mincellarea  = p.Results.mincellarea;

% bgn=50;
% lthresh=3000;
% projthresh=0.2;
% rejthresh=0.8;
% cellmasksize=8;
% cellratio=0.5;
% maxcellrange=80;
% mincellarea=25;


%%
addpath('C:\Fiji.app\scripts');
javaaddpath 'C:\Program Files\MATLAB\R2018b\java\mij.jar'

%segment cells using barcodes
files=dir('*aligned*BC*.tif');
files=sort_nat({files.name});

info=imfinfo(files{1});
im=zeros(info(1).Height,info(1).Width,4,length(files));
for i=1:length(files)
    for n=1:4
        im(:,:,n,i)=imread(files{i},n);
    end
end
%img=gpuarray(im);
img=im;
%find cells in each cycle/channel
imblurg=imgaussfilt(img,3);
%imblur=reshape(gather(imblurg),size(imblurg,1),size(imblurg,2),[]);
imblur=reshape(imblurg,size(imblurg,1),size(imblurg,2),[]);




%% normalize each cycle,additionally get intensity profile.
im1g=(imblurg+bgn)./repmat(sqrt(sum((imblurg+bgn).^2,3)),1,1,4,1);
im1g(isnan(im1g))=0;
im1g=reshape(im1g,size(im1g,1),size(im1g,2),[]);
%im1=gather(im1g);
im1=(im1g);
%im1
intimg=max(max(img,[],3),[],4);

%% at each "cell", find all surrounding pixels with similar colors.
% %sort currol by intensity and get intensity and normalized signals
% intim=gather(intimg);
% currrolint=zeros(size(currrol,1),1);
% currrolsig=zeros(size(currrol,1),1,size(imblur,3));
% for i=1:size(currrol,1)
%     currrolint(i)=intim(currrol(i,2)+1,currrol(i,1)+1);
%     currrolsig(i,1,:)=im1(currrol(i,2)+1,currrol(i,1)+1,:);
% end
% currrolsig=squeeze(currrolsig);
%
%
% [currrolint,I]=sort(currrolint,'descend');
% currrol=currrol(I,:);
% currrolsig=currrolsig(I,:);
%
% %for each "cell", measure distance to cell for neighboring pixels in
% %101x101 area centered around the "cell".
% segim=zeros(size(im1,1),size(im1,2));
% for i=1:size(currrol,1)
%     if segim(currrol(i,2)+1,currrol(i,1)+1)==0
%         subim=im1(max(1,currrol(i,2)-149):min(currrol(i,2)+151,size(im1,1)),max(1,currrol(i,1)-149):min(currrol(i,1)+151,size(im1,1)),:);
%         d=reshape(pdist2(reshape(subim,[],size(subim,3)),currrolsig(i,:),'euclidean'),size(subim,1),size(subim,2));
%         subint=intim(max(1,currrol(i,2)-149):min(currrol(i,2)+151,size(im1,1)),max(1,currrol(i,1)-149):min(currrol(i,1)+151,size(im1,1)));
%         subseg=(segim(max(1,currrol(i,2)-149):min(currrol(i,2)+151,size(im1,1)),max(1,currrol(i,1)-149):min(currrol(i,1)+151,size(im1,1)))==0) ...
%             .*(d<=distthresh)*i;
%         segim(max(1,currrol(i,2)-149):min(currrol(i,2)+151,size(im1,1)), ...
%             max(1,currrol(i,1)-149):min(currrol(i,1)+151,size(im1,1)))= ...
%             segim(max(1,currrol(i,2)-149):min(currrol(i,2)+151,size(im1,1)), ...
%             max(1,currrol(i,1)-149):min(currrol(i,1)+151,size(im1,1)))+subseg;
%     end
% end
%
%
% rgb=label2rgb(segim);
% figure;imshow(rgb);


%% intensity based. For each "cell", project neighboring pixels onto the vector of the "cell", then use a intensity threshold segment
if ~isempty(currrol1)
    %intim=gather(max(max(img,[],3),[],4));
    intim=(max(max(img,[],3),[],4));
    currrolint=zeros(size(currrol1,1),1);
    currrolsig=zeros(size(currrol1,1),1,size(imblur,3));
    for i=1:size(currrol1,1)
        currrolint(i)=intim(currrol1(i,2)+1,currrol1(i,1)+1);
        currrolsig(i,1,:)=imblur(currrol1(i,2)+1,currrol1(i,1)+1,:);%non-normalized signal
        
    end
    currrolsig=squeeze(currrolsig);
    %basecall "cells"
    sig1=reshape(currrolsig,size(currrolsig,1),4,[]);
    [~,seq]=max(sig1,[],2);
    seq=squeeze(seq);
    
    %sort "cells" by intensity
    [currrolint,I]=sort(currrolint,'descend');
    currrol1=currrol1(I,:);
    currrolsig=currrolsig(I,:);
    
    
    
    
    %pad both image and segmentation
    segim=zeros(size(im1,1)+300,size(im1,2)+300);
    imblurpad=padarray(imblur,[150 150],'both');
    
    for i=1:size(currrol1,1)
        if segim(currrol1(i,2)+151,currrol1(i,1)+151)==0 %&& ~ismember(seq(i,:),seq(uniqid(uniqid>0),:),'rows')
            subim=imblurpad(currrol1(i,2)+1:currrol1(i,2)+301,currrol1(i,1)+1:currrol1(i,1)+301,:);
            subim1=reshape(subim,[],size(subim,3));
            %project signal vectors of each pixel onto the "cell" signal vector
            proj=sum(subim1.*repmat(currrolsig(i,:),size(subim1,1),1),2)/sqrt(sum(currrolsig(i,:).^2));
            rejv=subim1-repmat(proj,1,size(currrolsig,2)).*repmat(currrolsig(i,:),size(subim1,1),1)/sqrt(sum(currrolsig(i,:).^2));
            rej=sqrt(sum(rejv.^2,2));
            
            proj=reshape(proj,size(subim,1),size(subim,2));
            rej=reshape(rej,size(subim,1),size(subim,2));
            %group pixels to that "cell" if the projection is above a threshold
            %and if the rejection/projection ratio is small
            subseg=(segim(currrol1(i,2)+1:currrol1(i,2)+301,currrol1(i,1)+1:currrol1(i,1)+301)==0) ...
                .*(proj>=projthresh*currrolint(i)).*(rej./proj<=rejthresh);
            %remove pixels not connected to the "cell"
                cc=bwconncomp(subseg>0);
            if ~isempty(cc.PixelIdxList)
                idx=zeros(length(cc.PixelIdxList),1);
                for n=1:length(idx)
                    idx(n)=ismember(sub2ind(size(proj),151,151),cc.PixelIdxList{n});
                end
                subseg(:)=0;
                if sum(idx)>0
                    subseg(cc.PixelIdxList{find(idx)})=1;
                    
                    cellmask=padarray(fspecial('disk',cellmasksize)>0, [150-cellmasksize,150-cellmasksize],'both');
                    [y,x]=ind2sub(size(subseg),cc.PixelIdxList{find(idx)});
                    %filter by cell size
                    if 0
                    elseif length(cc.PixelIdxList{find(idx)})<mincellarea %cells that are too small
                        subseg(:)=0;
                    elseif  sum(sum((subseg>0).*cellmask))<cellratio*sum(cellmask(:)>0) %if not enough pixels surrounding the "cell" has the same barcode
                        subseg(:)=0;
                        %elseif sum(sum((subseg>0).*cellmask))/sum(subseg(:)>0)<cellratio1
                        subseg(:)=0;
                    elseif range(y)+range(x)>maxcellrange %cells that take up too big a range
                        subseg(:)=0;
                    end
                end
                %fix cell shape
                %subseg=imopen(subseg,strel('disk',10,0));
                subseg=imclose(subseg,strel('disk',10,0));
                %if touching, check sequence and eliminate
                
                %subseg=imdilate(subseg,strel('disk',3,0));
                
                subseg=subseg*i.*(segim(currrol1(i,2)+1:currrol1(i,2)+301,currrol1(i,1)+1:currrol1(i,1)+301)==0);
                segim(currrol1(i,2)+1:currrol1(i,2)+301,currrol1(i,1)+1:currrol1(i,1)+301)= ...
                    segim(currrol1(i,2)+1:currrol1(i,2)+301,currrol1(i,1)+1:currrol1(i,1)+301)+subseg;
            end
        end
    end
    
    %rgb=label2rgb(segim,'jet','w','shuffle');
    %figure;imshow(rgb);
    %imwrite(rgb,'segmentation.tif');
    
    %remove padding and save.
    segim=segim(151:end-150,151:end-150);
else
    segim=zeros(size(im1,1),size(im1,2));
end


save('segmentation.mat','segim');









