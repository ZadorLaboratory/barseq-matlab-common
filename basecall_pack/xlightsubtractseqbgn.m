function xlightsubtractseqbgn(filename,f)
%subtract cleaved images from sequencing images. To fix for residual
%fluorescence. Assumes that odd number images are the sequencing images and
%even number images are after cleavage.f is multiplied to the background
%before subtraction to adjust intensity.


%read all tif file names in dirname folder
files=dir(fullfile(['*',filename,'*.tif']));
files=sort_nat({files.name});
mkdir original;

idx=1:2:length(files);
parfor k=1:length(idx)
    if idx(k)+1<=length(files) %if a cleavage image is available, subtract seq image. Otherwise, subtract the previous cleavage image
        bgnsub(files{idx(k)},files{idx(k)+1},f);
    else
        bgnsub(files{idx(k)},files{idx(k)-1},f);
    end
end
for k=1:length(files)
    movefile(files{k},['original/',files{k}]);
end
end

% This subfunction aligns cleavage image to seq image and subtract.
function bgnsub(seqfile,bgnfile,f)
%read file and median filter
for i=1:4
    imseq(:,:,i)=double(medfilt2(imread(seqfile,i)));
end
for i=1:4
    imbgn(:,:,i)=double(medfilt2(imread(bgnfile,i)));
end

imseqsum=sum(imseq,3);
imbgnsum=sum(imbgn,3);

par.transform = 'euclidean';
par.levels = 3;
par.iterations = 100;
ransacWarp=iat_ecc(imbgnsum,imseqsum,par);

[M,N]=size(imseq(:,:,1));
aligned=cell(4,1);
imseqbgnsub=imseq;
for i=1:4
    [aligned{i},~]=iat_inverse_warping(imbgn(:,:,i),ransacWarp,par.transform,1:N, 1:M);
    imseqbgnsub(:,:,i)=max(double(imseq(:,:,i))-double(aligned{i})*f,0);
end
% [optimizer,metric] = imregconfig('multimodal');
% tform = imregtform(imbgnsum(100:end-100,100:end-100), imseqsum(100:end-100,100:end-100), 'rigid', optimizer, metric);
% 
% Rfixed=imref2d(size(imseqsum));
% for i=1:4
%     aligned{i}=imwarp(imbgn(:,:,i),tform,'OutputView',Rfixed);
%     imseqbgnsub(:,:,i)=max(double(imseq(:,:,i))-double(aligned{i})*f,0);
% end


fixedimage=strcat('bgnsub',seqfile);
outputfile=fixedimage;
imwrite(uint16(imseqbgnsub(:,:,1)),outputfile);
for i=2:4
    imwrite(uint16(imseqbgnsub(:,:,i)), outputfile, 'WriteMode','append');
end
end


