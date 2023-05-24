function mmbasecallhyboptimize(codebookhyb,hybthresh1,bgn,psf,chprofile)
%basecalling gene rolonies based on codebook using a soft basecalling based
%on rolonies identified in the first cycle based on intensity. works on
%single position images in the aligned folder. If psf is not provided,
%guess psf using deconvblind.chprofile provides channel bleedthrough
%profile if images have not been fixed completely.



%% parse inputs
codes=num2str(cell2mat(codebookhyb(:,2)));

if ~exist('hybthresh1','var')
    hybthresh1=repmat([5 10 20 30 50 80 100 200 300 400 500 600 700 800 900 1000 ]',1,4);
elseif isempty(hybthresh1)
    hybthresh1=repmat([5 10 20 30 50 80 100 200 300 400 500 600 700 800 900 1000 ]',1,4);
end

if ~exist('bgn','var')
    bgn=hybthresh1;
elseif isempty(bgn)
    bgn=hybthresh1;
end


if ~exist('chprofile')
    chprofile=eye(4);
end
%noise=uint16(rolthresh1/2);
noise=[1 1 1 1];

%% read 1st cycle and find rolonies
seqfiles=dir(fullfile('*hyb*.tif'));
[seqfiles,~]=sort_nat({seqfiles.name});
%%

m=1;

imoriginal=[];


for n=1:4
    imoriginal(:,:,n)=imread(seqfiles{m},n);
end
imfixed=reshape(uint16(double(reshape(imoriginal,[],4))/chprofile),size(imoriginal,1),size(imoriginal,2),4);%subtract camera baseline and correct for bleeding
%%
parfor n=1:4
    lim(:,:,n)=deconvlucy(imgaussfilt(imfixed(:,:,n),1),psf{n},10); %current deconv settings are on the aggresive side.
    im(:,:,n)=imgaussfilt(imfixed(:,:,n),1);
end

%%
lroiall={};idall={};
parfor nn=1:size(hybthresh1,1)
    lroi1=[];id1=[];
    r={};rsub=r;
    for n=1:4
        a=lim(:,:,n);
        CC = bwconncomp(imregionalmax(imreconstruct(max(a-hybthresh1(nn,n),0),a)));
        r{n}=zeros(length(CC.PixelIdxList),1);
        for i=1:length(CC.PixelIdxList)
            [~,I]=max(a(CC.PixelIdxList{i}));
            r{n}(i)=CC.PixelIdxList{i}(I); %linear indexed peak positions
        end
        [y,x]=ind2sub(size(lim(:,:,n)),r{n});
        rsub{n}=[x,y];%rsub is in x, y, consistent with lroi1
    end
    for n=1:4
        c=find(codes(:,m)==num2str(n));
        if ~isempty(c)
            a=im(:,:,n);
            lroi1=[lroi1;rsub{n}(a(r{n})>bgn(nn,n),:)];
            id1=[id1;repmat(c,sum(a(r{n})>bgn(nn,n)),1)];
        end
    end
    lroiall{nn}=lroi1;
    idall{nn}=id1;
end


%% Generate rolony plots on top of images.
iplot=ceil(size(hybthresh1,1)/4*[1 2 3 4]);
for nn=1:size(codes,1)
    n=str2double(codes(nn));
    a=im(:,:,n);
    maxthresh=prctile(a(:),99.9);
    figure;
    s=zeros(size(hybthresh1,1),1);
    for i=1:size(hybthresh1,1)
        if ismember(i,iplot)
            
            subplot(2,4,find(iplot==i,1))
            imshow(a,[min(a(:)),maxthresh])
            hold on;
            scatter(lroiall{i}(idall{i}==nn,1),lroiall{i}(idall{i}==nn,2),'+r');
            set(gca,'xlim',[500 800],'ylim',[500 800])
            title(['Channel ',num2str(n),' thresh ',num2str(hybthresh1(i,n))]);
            subplot(2,4,find(iplot==i,1)+4)
            imshow(a,[min(a(:)),maxthresh])
            set(gca,'xlim',[500 800],'ylim',[500 800])
        end
        s(i)=sum(idall{i}==nn);
    end
    figure;plot(hybthresh1(:,n),s);
    title(['Channel ',num2str(n)])
end

end