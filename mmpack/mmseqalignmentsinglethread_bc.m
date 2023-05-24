function [chprofile,tforms,warnmsg]=mmseqalignmentsinglethread_bc(filename,chprofile, bgnradius,shiftmat)
%fix bleedthrough in all tif images in the current folder and align them.
%This version works for the xlight seq-o-matic. This version requires bleedthrough profile 6/10/2019. This version uses a
%single core. This version preserves non-seq frames and allows the use of
%dic image to pre-align. This version adds the option of shifting channels
%to correct for dichroic thickness. shiftmat is n x 2 for n channels. Each
%row is [Tx Ty].
%This version registers on the original images, works better for sparse
%signals


lastwarn('');

%check if channels need to be shifted
if ~exist('shiftmat','var')
    shiftmat=zeros(4,2);
end
    


%read all tif file names in dirname folder
files=dir(fullfile([filename,'*.tif']));
files=sort_nat({files.name});
L=size(files);
fixedfilelist=cell(L(2));
mkdir original;
% Fixes bleedthrough for all files, and store output filename list
%k=1;currentfile=files{k};
%[chprofile,fixedfilelist{k}]=fixbleed(currentfile,chprofile,bgnradius);
%movefile(files{k},['original/',files{k}]);
parfor k=1:L(2)
    currentfile=files{k};
    [~,fixedfilelist{k}]=fixbleed(currentfile,chprofile,bgnradius,shiftmat);
    movefile(files{k},['original/',files{k}]);
end

copyfile(fixedfilelist{1},['aligned',fixedfilelist{1}]);
parfor k=1:L(2)
    tforms{k}=alignseq(fixedfilelist{k},['original/',files{k}],['original/',files{1}]);
end

movefile('aligned*.tif','aligned/');
delete('fixed*tif');


warnmsg=lastwarn;
end

% align all images in fixedfilelist to fixedfilelist(1).





% This subfunction fixes bleed through between channels in illumina seq tiff file.
function [chprofile,fixedimage]=fixbleed(filename,chprofile,radius,shiftmat)
%read file and median filter
info=imfinfo(filename);
im=[];
for i=1:4
    %im(:,:,i)=medfilt2(imread(filename,i));
    im(:,:,i)=imread(filename,i);
end
if length(info)>4
    for i=5:length(info)
        im1(:,:,i-4)=imread(filename,i);
    end
end

%shift chanels
for i=1:size(shiftmat,1)
    im(:,:,i)=imtranslate(im(:,:,i),shiftmat(i,:));
end

%fix n2v artifacts in corners
im(im>60000)=0;%camera is 12 bit, so anything close to the 16-bit range is artifact


%bgn subtraction
%radius=20;
ball=strel('ball', radius, radius);
im=im-imopen(im,ball);

% Assign bleed through coefficients and bgn values
%chprofile=[1 0.76 0 0;0.19 1 0 0;0 0 1 0.37;0 0 0.48 1]';

im=reshape(uint16(double(reshape(im,[],4))/chprofile),size(im,1),size(im,2),4);%subtract camera baseline and correct for bleeding


fixedimage=strcat('fixed',filename);
outputfile=fixedimage;
imwrite(uint16(im(:,:,1)),outputfile);
for i=2:4
    imwrite(uint16(im(:,:,i)), outputfile, 'WriteMode','append');
end
if length(info)>4
    for i=5:length(info)
        imwrite(uint16(im1(:,:,i-4)),outputfile,'WriteMode','append');
    end
end

end

function tform=alignseq(refimagename,imagename,templatename)
% dicframe=5;
% imagedic=imread(imagename,dicframe);
% templatedic=imread(templatename,dicframe);
% 
% [optimizer,metric] = imregconfig('multimodal');
% optimizer.InitialRadius=optimizer.InitialRadius/5;
% tform = imregtform(imagedic(100:end-100,100:end-100), templatedic(100:end-100,100:end-100), 'rigid', optimizer, metric);
% Rfixed=imref2d(size(templatedic));

%read images and calculate sum of seq channels. 
info=imfinfo(imagename);
im=zeros(info(1).Height,info(1).Width,length(info));
for n=1:length(info)
    %im(:,:,n)=imwarp(imread(imagename,n),tform,'OutputView',Rfixed);
    im(:,:,n)=imread(imagename,n);
end



imagesum=double(sum(im(:,:,1:4),3))./double(max(max(sum(im(:,:,1:4),3))));

template1=imread(templatename,1);
template2=imread(templatename,2);
template3=imread(templatename,3);
template4=imread(templatename,4);
templatesum=double(template1+template2+template3+template4)./max(max(double(template1+template2+template3+template4)));

% %alignment using ECC
% par.transform = 'translation';
% par.levels = 6;
% par.iterations = 200;
% ransacWarp=iat_ecc(imagesum(200:end-200,200:end-200),templatesum(200:end-200,200:end-200),par);
% [M,N]=size(template4);


%align using imregtform or imregcorr
[optimizer,metric] = imregconfig('multimodal');
%optimizer.InitialRadius = optimizer.InitialRadius/5; %original conditions
%optimizer.MaximumIterations=optimizer.MaximumIterations*5;
optimizer.InitialRadius = optimizer.InitialRadius/5; %trial 1 conditions
%optimizer.GrowthFactor=1.01;
optimizer.Epsilon=optimizer.Epsilon;
optimizer.MaximumIterations=optimizer.MaximumIterations*5;
tform = imregtform(imagesum(100:end-100,100:end-100), templatesum(100:end-100,100:end-100), 'translation', optimizer, metric);
Rfixed=imref2d(size(templatesum));


info=imfinfo(refimagename);
im2=zeros(info(1).Height,info(1).Width,length(info));
for n=1:length(info)
    im2(:,:,n)=imread(refimagename,n);
end
alignedim=im2;

% for n=1:length(info)
%     [alignedim(:,:,n),~]=iat_inverse_warping(im1(:,:,n),ransacWarp,par.transform,1:N, 1:M);
% end
% tform=[];

for n=1:length(info)
    alignedim(:,:,n)=imwarp(im2(:,:,n),tform,'OutputView',Rfixed);
end


%write aligned image to file.
alignedfile=strcat('aligned',refimagename);
imwrite(uint16(alignedim(:,:,1)),alignedfile);
for i=2:size(im,3)
    imwrite(uint16(alignedim(:,:,i)),alignedfile, 'WriteMode','append');
end
end







