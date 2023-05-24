function alignnissl1(nisslfilename,cyclenum,dicnum,origpath)
%align nissl image to the sequencing image using one of dic image
%. nisslfilename=name of the nissl image (tif file with two channels, 
%first channel being the nissl image and second channel being dic
%cyclenum is the sequencing 
%cycle number to align to,chanelnum is the sequencing channel also imaged 
%in the nissl image. seqpath is the file path to the aligned sequencing 
%images. Saves the aligned nissl image in the same folder and output 
%paired image of the two aligned sequencing cycles for visual inspection 

%% align nissl image to first sequencing cycle.
imnissl=imread(nisslfilename);
imnissl(:,:,2)=imread(nisslfilename,2);
olddir=cd(origpath);
files=dir(fullfile('*.tif'));
files={files.name};
imseq=imread(files{cyclenum},dicnum);

%calculate alignment


        tform = imregcorr(imnissl(200:end-200,200:end-200,2),imseq(200:end-200,200:end-200),'rigid','Window',0);
        Rfixed=imref2d(size(imseq));
        alignedimnissl=imwarp(imnissl(:,:,1),tform,'OutputView',Rfixed);
        alignedimniss2=imwarp(imnissl(:,:,2),tform,'OutputView',Rfixed);
        
cd(olddir);
mkdir('nissl');
imwrite(uint16(alignedimnissl),'nissl/alignednissl.tif');
imwrite(uint16(alignedimniss2),'nissl/alignednissl.tif','WriteMode','Append');


