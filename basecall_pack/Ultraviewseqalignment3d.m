function [tforms,chprofiles]=Ultraviewseqalignment3d(filename,chprofile)
%fix bleedthrough in all tif images in the current folder and align them.
%updated with 20x bleedthrough profile and channel shift correction
%4/8/2018. This version uses FFT phase correlation for image alignment and
%cuts off the border for alignment. This version takes 3d stacks as inputs 
%but only aligns in the xy plane and uses. The z-alignment is done later
%locally during basecalling. This version also uses sparse nmf to find 
%channel bleedthrough. 

%fix and preprocess tif files:
seqfolders=dir(['*',filename,'*']);
mkdir original;

%check if a channel bleedthrough profile is provided. If not, initialize
%chprofile
if ~exist('chprofile','var')
    chprofile=zeros(4);
end

% Fixes bleedthrough for all files, and store output filename list
fixedfolders=cell(length(seqfolders),1);
for k=1:length(seqfolders)
    [fixedfolders{k},chprofiles{k}]=fixbleed3d(seqfolders(k).name,chprofile);
    movefile(seqfolders(k).name,['original/',seqfolders(k).name]);
end
%copy the first cycle fixed images directly as aligned images
copyfile(fixedfolders{1},['aligned',fixedfolders{1}]);
cd(['aligned',fixedfolders{1}]);
files=dir('*.tif');
for i=1:length(files)
    movefile(files(i).name,['aligned',files(i).name]);
end
cd ..

%chprofiles=cell(1);
%fixedfolders=dir('fixed*');
%fixedfolders={fixedfolders.name};

%align subsequent cycles to cycle 1
for k=2:length(fixedfolders)
    tforms{k}=alignseq(fixedfolders{k},fixedfolders{1});
end


%movefile('aligned*.tif','aligned/');
%delete('fixed*tif');

end


% This subfunction fixes bleed through between channels in illumina seq tiff file.
function [fixedfolder,chprofile]=fixbleed3d(seqfolders,chprofile)
cd(seqfolders);
files=dir('*.tif');
tcz=zeros(length(files),3);
for i=1:length(files)
    tcz(i,:)=cell2mat(textscan(files(i).name,'T%uC%uZ%u'));
end

im=cell(1,1,length(files));
for i=1:length(files)
    im{i}=imread(files(i).name);
end
im=cell2mat(im);


%correct for channel alignment
im(:,:,tcz(:,2)==1)=imtranslate(im(:,:,tcz(:,2)==1),[1.9 -0.4]);
%median filter
for i=1:length(files)
    im(:,:,i)=medfilt2(im(:,:,i));
end

% if chprofile is empty, find channel profile by sparse nmf on the center of the z-proj image
if sum(sum(chprofile))==0
    %find z-proj
    immax=zeros(size(im,1),size(im,2),4);
    for i=1:4
        immax(:,:,i)=max(im(:,:,tcz(:,2)==i),[],3);
    end
    %bgn subtraction
    radius=20;
    ball=strel('ball', radius, radius);
    bgn=imopen(immax, ball);
    immax=immax-bgn;
    %sparse-nmf
    [A,~]=sparsenmfnnls(reshape(double(immax(round(size(immax,1)/2)-500:round(size(immax,1)/2)+500, ...
        round(size(immax,2)/2)-500:round(size(immax,2)/2)+500,:)),[],4)',4); %Y should be correspond to the four channels. Can tune the sparseness constraint if necessary (doesn't seem to affect much)
    [~,I]=max(A,[],1);
    [~,I]=sort(I);
    chprofile=A(:,I)';%This assumes that for each color, the strongest channel is still the intended channel. (can't have two channels with max signal in the same imaging channel)
    chprofile=chprofile./repmat(max(chprofile,[],2),1,4);
end

%bgn subtract of stacked images and fix bleedthrough using channel
%profile
radius=20;
ball=strel('ball', radius, radius);
for i=1:max(tcz(:,3))
    bgn=imopen(im(:,:,tcz(:,2)<=4&tcz(:,3)==i), ball);
    im(:,:,tcz(:,2)<=4&tcz(:,3)==i)=im(:,:,tcz(:,2)<=4&tcz(:,3)==i)-bgn;
    im(:,:,tcz(:,3)==i & tcz(:,2)<=4)=reshape(uint16(double(reshape(im(:,:,tcz(:,3)==i & tcz(:,2)<=4),[],4))/chprofile),size(im,1),size(im,2),4);
    i
end

%channel profile measured manually. These are not used in this version
%CtoY=0.05;
%YtoC=0.35;
%MtoY=0.02;
%MtoW=0.84;
%WtoM=0.05;

mkdir(['../fixed',seqfolders]);
fixedfolder=['fixed',seqfolders];
cd ..
for i=1:length(files)
    imwrite(im(:,:,i),[fixedfolder,'/fixed',files(i).name]);
end
end

function tform=alignseq(imagefolder,templatefolder)
%read template sum across seq channels
cd(templatefolder);
filestemplate=dir('*.tif');
tcztemplate=zeros(length(filestemplate),3);
for i=1:length(filestemplate)
    tcztemplate(i,:)=cell2mat(textscan(filestemplate(i).name,'fixedT%uC%uZ%u'));
end

filestemplate(tcztemplate(:,2)>4)=[];%remove non-sequencing channels
tcztemplate(tcztemplate(:,2)>4,:)=[];%remove non-sequencing channels

imtemplate=cell(1,1,length(filestemplate));
for i=1:length(filestemplate)
    imtemplate{i}=imread(filestemplate(i).name);
end
imtemplate=cell2mat(imtemplate);
imsumtemplate=zeros(size(imtemplate,1),size(imtemplate,2),4);
for i=1:4
    imsumtemplate(:,:,i)=max(imtemplate(:,:,tcztemplate(:,2)==i),[],3); %max proj
end
imsumtemplate=sum(imsumtemplate,3);%z-proj of channel-summed template
clearvars imtemplate filestemplate tcztemplate %the non-summed template images are not used after this point

%read moving image and sum
cd(['../',imagefolder]);
files=dir('*.tif');
tcz=zeros(length(files),3);
for i=1:length(files)
    tcz(i,:)=cell2mat(textscan(files(i).name,'fixedT%uC%uZ%u'));
end

im=cell(1,1,length(files));
for i=1:length(files)
    im{i}=imread(files(i).name);
end
im=cell2mat(im);

imsum=zeros(size(im,1),size(im,2),4);
for i=1:4
    imsum(:,:,i)=max(im(:,:,tcz(:,2)==i),[],3); 
end
imsum=sum(imsum,3); %z-proj of channel-summed moving

%align moving sum to template sum
tform = imregcorr(imsum(200:end-200,200:end-200),imsumtemplate(200:end-200,200:end-200),'rigid','Window',0);
Rfixed=imref2d(size(imsumtemplate));

%write aligned image to file.
mkdir(['../aligned',imagefolder])
cd(['../aligned',imagefolder]);
for i=1:length(files)
    imwrite(uint16(imwarp(im(:,:,i),tform,'OutputView',Rfixed)),['aligned',files(i).name]);
end
cd ..
end







