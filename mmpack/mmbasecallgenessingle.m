function [id221,currsig1,currscore1,lroi1]=mmbasecallgenessingle(codebookbin,rolthresh1)
%basecalling gene rolonies by directly identifying rolonies from images. no
%intesity comparison between channels. Only works on non-combinatorial
%code.


%% parse inputs

if ~exist('rolthresh1')
    rolthresh1=[10 10 20 20];
end

rolrad=5;



%% read 1st cycle and find rolonies
seqfiles=dir(fullfile('*gene*seq*.tif'));
[seqfiles,~]=sort_nat({seqfiles.name});
lim=imread(seqfiles{1},1);
for n=2:4
    lim(:,:,n)=imread(seqfiles{1},n);
end

% %% try image processing
% %close all;
% %h=fspecial('gaussian',10,5);
% h = offsetstrel('ball',5,300);
% 
% 
% %blob=imfilter(double(lim(:,:,4)),h);
% blob=imgaussfilt(imopen(double(im(:,:,3)),h),1);
% figure;
% subplot(1,2,1);
% imagesc(double(im(:,:,3))./double(max(max(im(:,:,3)))),[0 0.2]);
% %set(gca,'xlim',[1250,1500],'ylim',[1050 1350]);
% 
% subplot(1,2,2);
% imagesc(blob./double(max(max(im(:,:,3)))),[0 0.2]);
% %set(gca,'xlim',[1250,1500],'ylim',[1050 1350]);
% 



%%
%find rolonies after gaussian blur with sigma=1
r={};
h = offsetstrel('ball',rolrad,300);

for n=1:4
    a=imgaussfilt(imopen(double(lim(:,:,n)),h),1);
    %a=imgaussfilt(lim(:,:,n),1);
    CC = bwconncomp(imregionalmax(imreconstruct(max(a-rolthresh1(n),0),a)));
    r{n}=zeros(length(CC.PixelIdxList),1); 
    for i=1:length(CC.PixelIdxList)
        [~,I]=max(a(CC.PixelIdxList{i}));
        r{n}(i)=CC.PixelIdxList{i}(I); %linear indexed peak positions
    end
end




%lpeaks(sub2ind([size(lim,1),size(lim,2)],currrol1(:,2)+1,currrol1(:,1)+1))=1;
for n=1:4
    lpeaks{n}=zeros(size(lim,1),size(lim,2));

    lpeaks{n}(r{n})=1;
    [lidxy{n},lidxx{n}]=find(lpeaks{n});
end

lroi1=[];
id221=[];
codes=double(cell2mat(codebookbin(:,2)));
[~,I]=max(codes(:,1:4),[],2);
for n=1:4
    lroi1=[lroi1;lidxx{n},lidxy{n}];
    id221=[id221;ones(length(lidxx{n}),1)*find(I==n)];
end

%% readout signals after deconvolution
currsig1=ones(length(lroi1),length(seqfiles),4);

for m=1:length(seqfiles) %for each cycle, read images and convolute by 7x7
    for n=1:4
        lim(:,:,n)=imread(seqfiles{m},n);
    end
    for n=1:length(lroi1) %for each rolony, readout signal
        currsig1(n,m,:)=lim(lroi1(n,2),lroi1(n,1),:);
    end
end

%calculate score
currscore1=max(currsig1,[],3)./sqrt(sum(currsig1.^2,3));
currscore1(isnan(currscore1))=0.5;


%figure;imshow(lim(:,:,1),[0 500]);

save('basecalls.mat','id221','currsig1','currscore1','lroi1');
end