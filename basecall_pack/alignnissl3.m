function alignnissl3(nisslnum,origpath)
%Produce aligned nissl image. This version directly gets the nissl channel
%from the first sequencing channel
%% align nissl image to first sequencing cycle.
olddir=cd(origpath);
files=dir(fullfile('*.tif'));
files={files.name};
imseq=imread(files{1},nisslnum);

%calculate alignment

cd(olddir);
mkdir('nissl');
imwrite(imseq,'nissl/alignednissl.tif');
%imwrite(uint16(alignedimniss2),'nissl/alignednissl.tif','WriteMode','Append');


