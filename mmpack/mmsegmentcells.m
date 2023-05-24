function segmentbc()
bgn=50;
lthresh=800;



addpath('C:\Fiji.app\scripts');
javaaddpath 'C:\Program Files\MATLAB\R2018b\java\mij.jar'

%segment cells using barcodes
files=dir('*aligned*BC*.tif');
files=sort_nat({files.name});

info=imfinfo(files{1});
im=zeros(info(1).Height,info(1).Width,size(info,1),numel(files));
for i=1:length(files)
    for n=1:size(info,1)
        im(:,:,n,i)=imread(files{i},n);
    end
end
img=gpuArray(im);
%find cells in each cycle/channel
imblurg=imgaussfilt(img,3);
imblur=reshape(gather(imblurg),size(imblurg,1),size(imblurg,2),[]);
currrol=[];
%find rolonies after gaussian blur with sigma=3
Miji(false);
for n=1:size(img,3)
    MIJ.createImage(imblur(:,:,n));
    MIJ.run('Find Maxima...', ['noise=',num2str(lthresh),' output=List']);
    currrol=[currrol;MIJ.getResultsTable];
    %cleanup
    MIJ.run('Clear Results');
    MIJ.closeAllWindows
end

MIJ.exit;



%% normalize each cycle,additionally get intensity profile.
imblurg1=reshape(imblurg,size(imblurg,1),size(imblurg,2),[]);
imgradg=imblurg1;
img1=reshape(img,size(img,1),size(img,2),[]);
for i=1:size(imgradg,3)
    imgradg(:,:,1)=imgradient(img1(:,:,i));
end

imgradg1=imblur;
for i=1:size(imblur,3)
    imgradg(:,:,1)=imgradient(imblur(:,:,i));
end




im1g=(img+bgn)./repmat(sqrt(sum((img+bgn).^2,3)),1,1,size(info,1),1);
im1g(isnan(im1g))=0;
intimg=sum(sum(img,3),4);

%% at each "cell", find all surrounding pixels with similar colors.
%sort currol by intensity
intim=gather(intimg);
currrolint=zeros(size(currrol,1),1);
for i=1:size(currrol,1)
    currrolint(i)=intim(currrol(i,2)+1,currrol(i,1)+1);
end

[currrolint,I]=sort(currrolint,'descend');
currrol=currrol(I,:);





