%fix bleedthrough in all tif images in the current folder and align them.
%This version works for the xlight seq-o-matic. bleed through profile is
%preliminary and may change in the future. 3/1/2019. This version uses
%parfor

function [chprofile,tforms]=xlightseqalignment(filename,chprofile)
%read all tif file names in dirname folder
files=dir(fullfile(['*',filename,'*.tif']));
files=sort_nat({files.name});
L=size(files);
fixedfilelist=cell(L(2));
mkdir original;
% Fixes bleedthrough for all files, and store output filename list
k=1;currentfile=files{k};
[chprofile,fixedfilelist{k}]=fixbleed(currentfile,chprofile);
movefile(files{k},['original/',files{k}]);
parfor k=2:L(2)
    currentfile=files{k};
    [~,fixedfilelist{k}]=fixbleed(currentfile,chprofile);
    movefile(files{k},['original/',files{k}]);
end

copyfile(fixedfilelist{1},['aligned',fixedfilelist{1}]);
parfor k=2:L(2)
    currentfile=fixedfilelist{k};
    tforms{k}=alignseq(currentfile,fixedfilelist{1});
end
%rgbout('bgnsub');
%movefile('RGB*.tif','rgb/');
%movefile('p*.tif','original/');

movefile('aligned*.tif','aligned/');
delete('fixed*tif');

%rgbout('aligned');
%movefile('RGB*.tif','rgb/');
%movefile('rough*.tif','rough/');
%movefile('aligned*.tif','aligned/');
end

% align all images in fixedfilelist to fixedfilelist(1).





% This subfunction fixes bleed through between channels in illumina seq tiff file.
function [chprofile,fixedimage]=fixbleed(filename,chprofile)
%read file and median filter
for i=1:4
    im(:,:,i)=medfilt2(imread(filename,i));
end


%bgn subtraction
radius=20;
ball=strel('ball', radius, radius);
im=im-imopen(im,ball);

% Assign bleed through coefficients and bgn values
%chprofile=[1 0.76 0 0;0.19 1 0 0;0 0 1 0.37;0 0 0.48 1]';

% if chprofile is empty, find channel profile by sparse nmf on the center of the z-proj image
%sparse-nmf
if isempty(chprofile)
    [A,~]=sparsenmfnnls(reshape(double(im(200:end-200,200:end-200,:)),[],4)',4); %Y should be correspond to the four channels. Can tune the sparseness constraint if necessary (doesn't seem to affect much)
    [~,I]=max(A,[],1);
    [~,I]=sort(I);
    chprofile=A(:,I)';%This assumes that for each color, the strongest channel is still the intended channel. (can't have two channels with max signal in the same imaging channel)
    chprofile=chprofile./repmat(max(chprofile,[],2),1,4);
end

im=reshape(uint16(double(reshape(im,[],4))/chprofile),size(im,1),size(im,2),4);%subtract camera baseline and correct for bleeding


fixedimage=strcat('fixed',filename);
outputfile=fixedimage;
imwrite(uint16(im(:,:,1)),outputfile);
for i=2:4
    imwrite(uint16(im(:,:,i)), outputfile, 'WriteMode','append');
end
end

function tform=alignseq(imagename,templatename)

%read images and calculate sum of seq channels. 
image1=imread(imagename,1);
image2=imread(imagename,2);
image3=imread(imagename,3);
image4=imread(imagename,4);
%correctedimage4=max(0,image4-1000)*1.2;
imagesum=uint16(65000*double(image1+image2+image3+image4)./double(max(max(image1+image2+image3+image4))));
template1=imread(templatename,1);
template2=imread(templatename,2);
template3=imread(templatename,3);
template4=imread(templatename,4);
templatesum=uint16(65000*double(template1+template2+template3+template4)./double(max(max(template1+template2+template3+template4))));

%alignment
par.transform = 'euclidean';
par.levels = 3;
par.iterations = 100;
ransacWarp=iat_ecc(imagesum(200:end-200,200:end-200),templatesum(200:end-200,200:end-200),par);

[M,N]=size(template4);

%transform images of nth cycle. 
[alignedimage1,~]=iat_inverse_warping(image1,ransacWarp,par.transform,1:N, 1:M);
[alignedimage2,~]=iat_inverse_warping(image2,ransacWarp,par.transform,1:N, 1:M);
[alignedimage3,~]=iat_inverse_warping(image3,ransacWarp,par.transform,1:N, 1:M);
[alignedimage4,~]=iat_inverse_warping(image4,ransacWarp,par.transform,1:N, 1:M);

tform=[];
%write aligned image to file.
alignedfile=strcat('aligned',imagename);
imwrite(uint16(alignedimage1), alignedfile);
imwrite(uint16(alignedimage2),alignedfile, 'WriteMode','append');
imwrite(uint16(alignedimage3),alignedfile, 'WriteMode','append');
imwrite(uint16(alignedimage4),alignedfile, 'WriteMode','append');
end







