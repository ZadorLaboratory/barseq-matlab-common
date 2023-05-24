%fix bleedthrough in all tif images in the current folder.
%This version works for the xlight seq-o-matic. bleed through profile is
%preliminary and may change in the future. 3/1/2019.

function chprofile=xlightfiximages(filename)
%read all tif file names in dirname folder
files=dir(fullfile(['*',filename,'*.tif']));
files=sort_nat({files.name});
L=size(files);
fixedfilelist=cell(L(2));
mkdir original;
% Fixes bleedthrough for all files, and store output filename list
for k=1:L(2)
    currentfile=files{k};
    if k==1
        [chprofile,fixedfilelist{k}]=fixbleed(currentfile,[]);
    else
        [chprofile,fixedfilelist{k}]=fixbleed(currentfile,chprofile);
    end
    
    movefile(files{k},['original/',files{k}]);
end

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




