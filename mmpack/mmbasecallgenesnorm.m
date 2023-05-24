function [id21,currsig1,currscore1,lroi1,id31]=mmbasecallgenesnorm(codebookbin,rolthresh1,mapthresh,psf,bgn)
%basecalling gene rolonies based on codebook using a soft basecalling based
%on rolonies identified in the first cycle based on intensity. works on
%single position images in the aligned folder. If psf is not provided,
%guess psf using deconvblind. This version normalizes channel signals
%within a cycle.


%% parse inputs

if ~exist('rolthresh1')
    rolthresh1=[15 15 40 30];
end

if ~exist('mapthresh')
    mapthresh=0.25;
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


lpeaks=lpeaks&~imdilate(lpeaks,triu(ones(3))-diag([ones(1,2),zeros(1,1)]));
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

% balance the channels in each cycle
currsig1=currsig1./repmat(mean(currsig1,1),size(currsig1,1),1,1).*repmat(mean(currsig1(:,:,1),1),size(currsig1,1),1,size(currsig1,3));


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
id21(qual1<mapthresh,:)=0;




% match to code book by hard code 
codes=double(cell2mat(codebookbin(:,2)));
codesmax=reshape(codes,size(codes,1),4,[]);
[~,codesmax]=max(permute(codesmax,[1 3 2]),[],3);

[~,currsigmax]=max(currsig1,[],3);
[~,id31]=ismember(currsigmax,codesmax(:,1:size(currsigmax,2)),'rows');
id31(min(currscore1,[],2)<0.8)=0;





% 
% currsig1=reshape(permute(currsig1,[1,3,2]),size(currsig1,1),[]);
% d=squareform(pdist(lroi1,'euclidean'));
% 
% % for each rolony with low basecalling quality, find neighboring rolonies
% [d2,I]=sort(d,2,'ascend');
% 
% u=fspecial('gaussian', [1 19], 7);%make a 1d gaussian filter to simulate a rolony
% u=u(10:end)/max(u);
% 
% d1=d;
% d1(d1>10)=0;
% d1(d1>0)=u(ceil(d1(d1>0)));
% 
% bgnsig=d1*currsig1;
% %subtract background for low quality rolonies and redo basecall
% currsig2=currsig1-bgnsig;
% currsig2(currsig2<0)=0;
% 
% currsig2=permute(reshape(currsig2,size(currsig2,1),4,[]),[1 3 2]);
% 
% currsignorm2=currsig2./repmat(sum(currsig2,3),1,1,4);
% currsignorm2=reshape(permute(currsignorm2,[1,3,2]),size(currsignorm2,1),[]);
% currsignorm2(isnan(currsignorm2))=0;
% 
% proj2=zeros(size(currsignorm2,1),size(codes,1));
% for i=1:size(codes,1)
%     proj2(:,i)=sum(currsignorm2.*repmat(codes(i,1:size(currsignorm2,2)),size(currsignorm2,1),1),2)/sqrt(sum(codes(i,1:size(currsignorm2,2)).^2));
% end
% [sortedproj2,id221]=sort(proj2,2,'descend');
% %id22(sortedproj2(:,1)*mapthresh<sortedproj(:,2),:)=0;
% qual2=1-sortedproj2(:,2)./sortedproj2(:,1); %range [0 1] with 0 being equaly matching the top two barcodes and 1 being only matching one barcode
% id221(qual1>qual2,:)=id21(qual1>qual2,:);
% id21(qual2<mapthresh,:)=0;
% 
% % figure;
% % subplot(1,3,1);
% % scatter(qual,qual2,'.');title('sigma 7 ');
% % subplot(1,3,2);scatter(d2(:,2),max([qual2,qual],[],2),'.');
% % subplot(1,3,3);scatter(d2(:,2),qual,'.'); %trying to clean up signal of the ones with nearby rolonies
% 
% 
% %filter out low quality barcodes
% qual1=max([qual2,qual1],[],2);
% id221(qual1<mapthresh,:)=0; %this correspond to roughly a quality score of 0.8
% 




% %% visualize rolonies
% sum(id21(:,1)==18)
% sum(id21(:,1)==17)
% sum(id21(:,1)>0)
% sum(id21(:,1)==0)
% 
% sum(id31(:,1)==18)
% sum(id31(:,1)==17)
% sum(id31(:,1)>0)
% sum(id31(:,1)==0)
% 
% sum(id31(:,1)==18|id21(:,1)==18)
% sum(id31(:,1)==17|id21(:,1)==17)
% sum(id31(:,1)>0|id21(:,1)>0)
% sum(id31(:,1)==0&id21(:,1)==0)
% 
% 
% cd ../RGB
% rgbfiles=dir('*gene*.tif');
% rgbfiles=sort_nat({rgbfiles.name});
% rgbim=imread(rgbfiles{1});
% defaultcolors=get(groot,'DefaultAxesColorOrder');
% figure;
% imshow(rgbim);hold on;
% i=0;
% scatter(lroi1(id21(:,1)==i,1),lroi1(id21(:,1)==i,2),'+r');
% scatter(lroi1(id21(:,1)>0,1),lroi1(id21(:,1)>0,2),'+g');
% scatter(lroi1(id21(:,1)>16,1),lroi1(id21(:,1)>16,2),'or');
% set(gca,'ydir','reverse');
% figure;
% imshow(rgbim);hold on;
% i=0;
% scatter(lroi1(id31(:,1)==i,1),lroi1(id31(:,1)==i,2),'+r');
% scatter(lroi1(id31(:,1)>0,1),lroi1(id31(:,1)>0,2),'+g');
% scatter(lroi1(id31(:,1)>16,1),lroi1(id31(:,1)>16,2),'or');
% set(gca,'ydir','reverse');

% 
% for i=1:4
%     subplot(1,4,i)
%     scatter(lroi1(id221(:,1)==i,1),lroi1(id221(:,1)==i,2),3,defaultcolors(i,:),'filled');
%     set(gca,'ydir','reverse');
% 
% end
% %as comparison call rolonies by first base intensity only
% l1=currsignorm(:,1:4)/double(codes(:,1:4));
% [~,id1]=sort(l1,2,'descend');
% figure;
% for i=1:4
%     subplot(1,4,i)
%     scatter(lroi1(id1(:,1)==i,1),lroi1(id1(:,1)==i,2),3,defaultcolors(i,:),'filled');
%     set(gca,'ydir','reverse');
%
% end
%
% figure;
% for i=1:4
%     subplot(1,4,i)
%     scatter(lroi(id22(:,1)==i,1),lroi(id22(:,1)==i,2),3,defaultcolors(i,:),'filled');
%     set(gca,'ydir','reverse');
%
% end





% %%
% mappedx=tsne(currsignorm);
% %%
% defaultcolors=get(groot,'DefaultAxesColorOrder');
% figure;hold on;
% for i=1:5
%     scatter(mappedx(id221(:,1)==i-1,1),mappedx(id221(:,1)==i-1,2),3,defaultcolors(i,:));
% end
% legend({'null','Slc17a7','Gad1','Slc30a3','Cdh13'});
% figure;hold on;
% for i=1:5
%     scatter(mappedx(id1(:,1)==i-1,1),mappedx(id1(:,1)==i-1,2),3,defaultcolors(i,:));
% end
% legend({'null','Slc17a7','Gad1','Slc30a3','Cdh13'});
%
% figure;hold on;
% for i=1:5
%     scatter(mappedx(:,1),mappedx(:,2),3,qual,'filled');
% end
% legend({'null','Slc17a7','Gad1','Slc30a3','Cdh13'});

%% combine basecalling and save
id21(id21(:,1)==0&id31>0,1)=id31(id21(:,1)==0&id31>0);

% %% check intensity of four channels
% [maxsig,currsigmax]=max(currsig1,[],3);
% ratio(1)=mean(maxsig(currsigmax(:,1)==1&id21(:,1)>0,1))
% ratio(2)=mean(maxsig(currsigmax(:,1)==2&id21(:,1)>0,1))
% ratio(3)=mean(maxsig(currsigmax(:,1)==3&id21(:,1)>0,1))
% ratio(4)=mean(maxsig(currsigmax(:,1)==4&id21(:,1)>0,1))
% ratio=ratio/max(ratio);
% 
% n1=histcounts(maxsig(currsigmax(:,4)==1&id31(:,1)>0,1),0:50:2000);
% n2=histcounts(maxsig(currsigmax(:,4)==2&id31(:,1)>0,1),0:50:2000);
% n3=histcounts(maxsig(currsigmax(:,4)==3&id31(:,1)>0,1),0:50:2000);
% n4=histcounts(maxsig(currsigmax(:,4)==4&id31(:,1)>0,1),0:50:2000);
% figure;plot(n1/max(n1))
% hold on;
% plot(n2/max(n2));plot(n3/max(n3));plot(n4/max(n4));
% legend({'1','2','3','4'},'Location','eastoutside');
% normlize currsig by mean intensity of channels and recall bases
% currsig2=currsig1./repmat(reshape(ratio,1,1,[]),size(currsig1,1),size(currsig1,2),1);
% currscore1=max(currsig2,[],3)./sqrt(sum(currsig2.^2,3));
% currscore1(isnan(currscore1))=0.5;
% 
% currsignorm=currsig2./repmat(sum(currsig2,3),1,1,4);
% currsignorm=reshape(permute(currsignorm,[1,3,2]),size(currsignorm,1),[]);
% currsignorm(isnan(currsignorm))=0;
% 
% match to codebook by projecting to each barcode, then
% codes=double(cell2mat(codebookbin(:,2)));
% proj1=zeros(size(currsignorm,1),size(codes,1));
% for i=1:size(codes,1)
%     proj1(:,i)=sum(currsignorm.*repmat(codes(i,1:size(currsignorm,2)),size(currsignorm,1),1),2)/sqrt(sum(codes(i,1:size(currsignorm,2)).^2));
% end
% [sortedproj,id21]=sort(proj1,2,'descend');
% qual1=1-sortedproj(:,2)./sortedproj(:,1); %range [0 1] with 0 being equaly matching the top two barcodes and 1 being only matching one barcode
% id21(qual1<mapthresh,:)=0;
% 
% match to code book by hard code 
% codes=double(cell2mat(codebookbin(:,2)));
% codesmax=reshape(codes,size(codes,1),4,[]);
% [~,codesmax]=max(permute(codesmax,[1 3 2]),[],3);
% 
% [~,currsigmax]=max(currsig2,[],3);
% [~,id31]=ismember(currsigmax,codesmax(:,1:size(currsigmax,2)),'rows');
% id31(min(currscore1,[],2)<0.7)=0;
% 

% %% check rolonies against image
% cd ../RGB
% rgbfiles=dir('*gene*.tif');
% rgbfiles=sort_nat({rgbfiles.name});
% rgbim=imread(rgbfiles{1});
% defaultcolors=get(groot,'DefaultAxesColorOrder');
% figure;
% imshow(rgbim);hold on;
% i=3;
% scatter(lroi1(id21(:,1)==i,1),lroi1(id21(:,1)==i,2),'+r');
% %scatter(lroi1(id21(:,1)==0,1),lroi1(id21(:,1)==0,2),'+g');
% %scatter(lroi1(id21(:,1)>16,1),lroi1(id21(:,1)>16,2),'or');
% set(gca,'ydir','reverse');
% figure;
% imshow(rgbim);hold on;
% i=3;
% scatter(lroi1(id31(:,1)==i,1),lroi1(id31(:,1)==i,2),'+r');
% %scatter(lroi1(id31(:,1)==0,1),lroi1(id31(:,1)==0,2),'+g');
% %scatter(lroi1(id31(:,1)>16,1),lroi1(id31(:,1)>16,2),'or');
% set(gca,'ydir','reverse');
% 




%%
save('basecalls.mat','id21','currsig1','currscore1','id31','lroi1');
%
% MIJ.exit;
end