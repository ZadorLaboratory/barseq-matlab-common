function [chprofile,tforms,warnmsg]=mmseqalignmentsinglethread_addcycles(filename,chprofile, bgnradius,shiftmat,ch)
%fix bleedthrough in all tif images in the current folder and align them.
%This version preserves non-seq frames and allows the use of
%dic image to pre-align. This version adds the option of shifting channels
%to correct for dichroic thickness. shiftmat is n x 2 for n channels. Each
%row is [Tx Ty].
%this version adds cycles and align them to already processed first-cycle
%image. Can be used during imaging?
%ch chooses which channels to use for registration, array of 1:4 by default


lastwarn('');

%check if channels need to be shifted
if ~exist('shiftmat','var')
    shiftmat=zeros(4,2);
elseif isempty(shiftmat)
    shiftmat=zeros(4,2);
end

if ~exist('ch','var')
    ch=1:4;
end
    


%read all tif file names in dirname folder
files=dir(fullfile([filename,'*.tif']));
files=sort_nat({files.name});
L=size(files);
fixedfilelist=cell(L(2));

parfor k=1:L(2)
    currentfile=files{k};
    [~,fixedfilelist{k}]=fixbleed(currentfile,chprofile,bgnradius,shiftmat);
    movefile(files{k},['original/',files{k}]);
end
cd aligned
imfiles=dir(['aligned*',filename,'*.tif']);
imfiles={imfiles.name};
imfiles=imfiles{1};
cd ..
parfor k=1:L(2)
    currentfile=fixedfilelist{k};
    tforms{k}=alignseq(currentfile,['aligned/',imfiles],ch);
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

function tform=alignseq(imagename,templatename,ch)
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
alignedim=im;



imagesum=double(sum(im(:,:,ch),3))./double(max(max(sum(im(:,:,ch),3))));




template1=imread(templatename,1);
template2=imread(templatename,2);
template3=imread(templatename,3);
template4=imread(templatename,4);
templatesum=double(template1+template2+template3+template4)./max(max(double(template1+template2+template3+template4)));

% %pre-register using FFT
% tform=imregcorr(imagesum(100:end-100,100:end-100), templatesum(100:end-100,100:end-100),'rigid','Window',0);
% [optimizer,metric] = imregconfig('multimodal');
% tform = imregtform(imagesum,templatesum, 'translation', optimizer, metric);
% 
% Rfixed=imref2d(size(templatesum));
% for n=1:length(info)
%     im(:,:,n)=imwarp(im(:,:,n),tform,'OutputView',Rfixed);
% end
%  imagesum=imwarp(imagesum,tform,'OutputView',Rfixed);

%alignment using ECC
% par.transform = 'translation';
% par.levels = 6;
% par.iterations = 200;
% ransacWarp=iat_ecc(imagesum(50:end-50,50:end-50),templatesum(50:end-50,50:end-50),par);
% [M,N]=size(template4);
% for n=1:length(info)
%     [alignedim(:,:,n),~]=iat_inverse_warping(im(:,:,n),ransacWarp,par.transform,1:N, 1:M);
% end
% tform=[];

%align using imregtform or imregcorr
[optimizer,metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/5;
optimizer.MaximumIterations=optimizer.MaximumIterations*5;
tform = imregtform(imagesum(100:end-100,100:end-100), templatesum(100:end-100,100:end-100), 'translation', optimizer, metric);
% tform=imregcorr(imagesum(100:end-100,100:end-100), templatesum(100:end-100,100:end-100),'rigid','Window',0);
Rfixed=imref2d(size(templatesum));
for n=1:length(info)
    alignedim(:,:,n)=imwarp(im(:,:,n),tform,'OutputView',Rfixed);
end


%write aligned image to file.
alignedfile=strcat('aligned',imagename);
imwrite(uint16(alignedim(:,:,1)),alignedfile);
for i=2:size(im,3)
    imwrite(uint16(alignedim(:,:,i)),alignedfile, 'WriteMode','append');
end
end







