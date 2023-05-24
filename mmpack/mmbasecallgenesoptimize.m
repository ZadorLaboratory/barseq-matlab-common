function mmbasecallgenesoptimize(rolthreshall,psf,bgn)
%optimize rolthresh for mmbasecallgenes


%% parse inputs

if ~exist('rolthreshall')
    rolthreshall=repmat([2 4 6 8 10 20 50 100 150 200 300]',1,4);
elseif isempty(rolthreshall)
    rolthreshall=repmat([2 4 6 8 10 20 50 100 150 200 300]',1,4);
end

if ~exist('bgn','var')
    bgn=rolthreshall;
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
        im(:,:,n)=deconvlucy(imgaussfilt(imread(seqfiles{1},n),1),psf{n},5); %current deconv settings are on the aggresive side.
    end
end


%%
%find rolonies after gaussian blur with sigma=1
rsuball={};
parfor nn=1:size(rolthreshall,1)
    rolthresh1=rolthreshall(nn,:);
    r={};rsub={};
    for n=1:4
        
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
    rsuball{nn}=rsub1;
end
%% Generate rolony plots on top of images.
for n=1:4
    a=im(:,:,n);
    maxthresh=prctile(a(:),99.9);
    figure;
    s=zeros(size(rolthreshall,1),1);
    for i=size(rolthreshall,1)-4:size(rolthreshall,1)
        subplot(2,size(rolthreshall,1),i)
        imshow(a,[min(a(:)),maxthresh])
        hold on;
        scatter(rsuball{i}{n}(:,1),rsuball{i}{n}(:,2),'+r');
        set(gca,'xlim',[500 800],'ylim',[500 800])
        title(['Channel ',num2str(n),' thresh ',num2str(rolthreshall(i,n))]);
        subplot(2,size(rolthreshall,1),i+size(rolthreshall,1))
        imshow(a,[min(a(:)),maxthresh])
        set(gca,'xlim',[500 800],'ylim',[500 800])
         s(i)=size(rsuball{i}{n},1);
    end
    figure;plot(rolthreshall(:,n),s);
    title(['Channel ',num2str(n)])
end
%
% MIJ.exit;
end