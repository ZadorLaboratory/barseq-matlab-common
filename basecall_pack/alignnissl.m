function alignnissl(nisslfilename,cyclenum,channelnum,seqpath)
%align nissl image to the sequencing image using one of the sequencing 
%channels. nisslfilename=name of the nissl
%image (tif file with two channels, first channel being the nissl image and
%second channel being a sequencing channel), cyclenum is the sequencing 
%cycle number to align to,chanelnum is the sequencing channel also imaged 
%in the nissl image. seqpath is the file path to the aligned sequencing 
%images. Saves the aligned nissl image in the same folder and output 
%paired image of the two aligned sequencing cycles for visual inspection 

%% align nissl image to last sequencing cycle.
imnissl=imread(nisslfilename);
imnissl(:,:,2)=imread(nisslfilename,2);
olddir=cd(seqpath);
files=dir(fullfile('*.tif'));
files={files.name};
imseq=imread(files{cyclenum},channelnum);

%preprocessing seq channel in nissl image to match sequencing image
radius=10;
ball=strel('ball', radius, radius);
imnissl1=max(medfilt2(imnissl(:,:,2))-imopen(medfilt2(imnissl(:,:,2)),ball),0);

%calculate alignment
par.transform = 'translation';
par.levels = 8;
par.iterations = 50;
ransacWarp=iat_ecc(imnissl1,imseq,par);

[M,N]=size(imseq);
[alignedimnissl,~]=iat_inverse_warping(imnissl(:,:,1),ransacWarp,par.transform,1:N, 1:M);

[alignedimniss2,~]=iat_inverse_warping(imnissl1,ransacWarp,par.transform,1:N, 1:M);
figure;imshowpair(double(alignedimniss2)/max(max(double(alignedimniss2))),double(imseq)/max(max(double(imseq))));

cd(olddir);
imwrite(uint16(alignedimnissl),'nissl/alignednissl.tif');

