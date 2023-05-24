function mmbasecallgenesoptimize2lowres(codebookbin,rolthresh1,mapthresh,psf,bgn,nctrl)
%basecalling gene rolonies based on codebook using a soft basecalling based
%on rolonies identified in the first cycle based on intensity. works on
%single position images in the aligned folder. If psf is not provided,
%guess psf using deconvblind.


%% parse inputs

if ~exist('rolthresh1','var')
    rolthresh1=[15 15 40 30];
end

if ~exist('mapthresh','var')
    mapthresh=0.05:0.01:0.3;
elseif isempty(mapthresh)
    mapthresh=0.05:0.01:0.3;
end

if ~exist('bgn','var')
   bgn=rolthresh1;
end



%noise=uint16(rolthresh1/2);
noise=[1 1 1 1];

% %% start fiji
% addpath('C:\Fiji.app\scripts');
% javaaddpath 'C:\Program Files\MATLAB\R2018b\java\mij.jar'
% Miji(false);

%% read 1st cycle and find rolonies
seqfiles=dir(fullfile('*gene*seq*.tif'));
[seqfiles,~]=sort_nat({seqfiles.name});
 lim=imread(seqfiles{1},1);
 for n=2:4
     lim(:,:,n)=imread(seqfiles{1},n);
 end



im=lim;
if ~exist('psf','var')
    parfor n=1:4
        [im(:,:,n),psf{n}]=deconvblind(imread(seqfiles{1},n),ones(11),3,uint16(noise(n)));
    end
else
    parfor n=1:4
        im(:,:,n)=deconvlucy(imgaussfilt(imread(seqfiles{1},n),1),psf{n},10); %current deconv settings are on the aggresive side.
    end
end

%%
%find rolonies after gaussian blur with sigma=1
r={};rsub={};
parfor n=1:4
    
    a=im(:,:,n);
    CC = bwconncomp(imregionalmax(imreconstruct(max(a-rolthresh1(n),0),a)));
    r{n}=zeros(length(CC.PixelIdxList),1);
    for i=1:length(CC.PixelIdxList)
        [~,I]=max(a(CC.PixelIdxList{i}));
        r{n}(i)=CC.PixelIdxList{i}(I); %linear indexed peak positions
    end
    [y,x]=ind2sub(size(im(:,:,n)),r{n});
    rsub{n}=[x,y];%rsub is in x, y, consistent with lroi1
    
    
end

for n=1:4
    %deconv amplifies noise. Remove by amplitude thresholding.
    b=lim(:,:,n);
    rsub{n}(b(r{n})<bgn(n),:)=[];
    r{n}(b(r{n})<bgn(n))=[];
end
% %% visualize rolonies
% i=4;
% figure;subplot(1,2,1);imagesc(lim(:,:,i),[0 300]);
% hold on;scatter(rsub{i}(:,1),rsub{i}(:,2),'+r');
% set(gca,'xlim',[1300 1600],'ylim',[1300 1600]);
% subplot(1,2,2);imagesc(lim(:,:,i),[0 300]);
% set(gca,'xlim',[1300 1600],'ylim',[1300 1600]);
% %Rolonies are mostly real,perhaps 10-20% are noise, dependent on channel.

%%
%combine rolonies
lpeaks=zeros(size(lim,1),size(lim,2));
%lpeaks(sub2ind([size(lim,1),size(lim,2)],currrol1(:,2)+1,currrol1(:,1)+1))=1;
for n=1:4
    lpeaks(r{n})=1;
end


%lpeaks=lpeaks&~imdilate(lpeaks,triu(ones(3))-diag([ones(1,2),zeros(1,1)]));
[lidxy,lidxx]=find(lpeaks);
lroi1=[lidxx,lidxy];

%figure out which channel these rolonies are from
lpeaks1=zeros(size(lim,1),size(lim,2));
for n=1:4
    lpeaks1(r{n})=n;
end

lpeaks1=lpeaks1.*(lpeaks>0);
rsub1={};       
for n=1:4
    [y,x]=ind2sub(size(lim(:,:,n)),find(lpeaks1(:)==n));
    rsub1{n}=[x,y];
end




%% readout signals after deconvolution
currsig1=ones(length(lroi1),length(seqfiles),4);

for m=1:length(seqfiles) %for each cycle, read images and convolute by 7x7
    if ~exist('psf','var')
        parfor n=1:4
            [lim(:,:,n),psf{n}]=deconvblind(imread(seqfiles{m},n),ones(11),10,uint16(noise(n)));
        end
    else
        parfor n=1:4
            lim(:,:,n)=imgaussfilt(deconvlucy(imread(seqfiles{m},n),psf{n},10),1);
        end
    end
    
    for n=1:length(lroi1) %for each rolony, readout signal
        currsig1(n,m,:)=lim(lroi1(n,2),lroi1(n,1),:);
    end
end

%figure;imshow(lim(:,:,1),[0 500]);




%%
%calculate score
currscore1=max(currsig1,[],3)./sqrt(sum(currsig1.^2,3));
currscore1(isnan(currscore1))=0.5;

currsignorm=currsig1./repmat(sum(currsig1,3),1,1,4);
currsignorm=reshape(permute(currsignorm,[1,3,2]),size(currsignorm,1),[]);
currsignorm(isnan(currsignorm))=0;

%% match to codebook by projecting to each barcode, then
codes=double(cell2mat(codebookbin(:,2)));
proj1=zeros(size(currsignorm,1),size(codes,1));
for i=1:size(codes,1)
    proj1(:,i)=sum(currsignorm.*repmat(codes(i,1:size(currsignorm,2)),size(currsignorm,1),1),2)/sqrt(sum(codes(i,1:size(currsignorm,2)).^2));
end
[sortedproj,id21]=sort(proj1,2,'descend');
qual1=1-sortedproj(:,2)./sortedproj(:,1); %range [0 1] with 0 being equaly matching the top two barcodes and 1 being only matching one barcode
id22={};
for i=1:numel(mapthresh)
    id22{i}=id21(:,1);
    id22{i}(qual1<mapthresh(i))=0;
end

figure;
subplot(1,3,1);
plot(mapthresh,cellfun(@(x) sum(ismember(x,nctrl))/numel(nctrl),id22));
title('Total errors per gene')
xlabel('Quality threshold');ylabel('Counts');
subplot(1,3,2);
plot(mapthresh,cellfun(@(x) sum(~ismember(x,nctrl)&x~=0),id22));
title('Total valid reads');
xlabel('Quality threshold');ylabel('Counts');
subplot(1,3,3);
plot(mapthresh,cellfun(@(x) sum(ismember(x,nctrl))/numel(nctrl)./sum(~ismember(x,nctrl)&x~=0)*(size(codebookbin,1)-numel(nctrl)),id22));
title('Predicted error rate')
xlabel('Quality threshold');ylabel('Fraction of errors');






end