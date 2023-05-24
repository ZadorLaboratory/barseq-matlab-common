function [imnuc,imcell]=findnuclei1(alignednisslfilename,blackpoint,bgnthresh,nucleisizethresh,cellsize)
%find nuclei and estimte cell location using aligned nuclei stain. Output
%nuc and cell with nuclei/cell labels. background is labeled 0. 
%   alignednisslfilename: file name of the aligned nissl image produced by
%   alignnissl.m
%   blackpoint: blackpoint of the image for initial background subtraction.
%   Can be set pretty aggresively
%   bgnthresh: additional threshold for background after blackpoint
%   subtraction, for finding backgrounds.
%   nucleisizethresh: max size of nuclei for filtering out multiple nuclei,
%   in # of pixels.
%   cellsize: the number of pixels to expand from each side of nuclei.


%% segment nuclei and expand to find cells

im=imread(alignednisslfilename,1);
%figure; imagesc(im);


%remove background
im1=max(im-blackpoint,0);

%find foreground (nuclei) by opening + erosion
se=strel('disk',5);
Io=imopen(im1,se);
%figure;
%imshow(Io(2000:2500,2000:2500),[0 100]);

se2 = strel('disk',3);
se21=strel('disk',4);
%fgm1 = imerode(imclose(imregionalmax(Iobr),se2),se21);
fgm1 = imerode(imregionalmax(Io),se21);

figure;imshowpair(double(im1(2000:2500,2000:2500))/max(double(im1(:))),fgm1(2000:2500,2000:2500));

im2=medfilt2(double(im1))-imopen(double(im1),strel('ball',50,50));

%find background by thresholding
%bgn=imbinarize(im1,0.1); %might need manual adjustment here
%figure;imshow(bgn);
bgn=(im2-bgnthresh)>0;
%thin background
D = bwdist(bgn);
DL = watershed(D);
bgm = DL == 0;
figure; imshow(bgm); title('Watershed ridge lines (bgm)')

%marker-controlled watershed on the gradient magnitude
L=watershed(imimposemin(imgradient(im1), bgm | fgm1));
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
%figure
%imshow(Lrgb)
%title('Nuclei(unfiltered)')

% filter out cells that are too large (likely multiple cells merged).
n=histcounts(L,[-0.5;unique(L)+0.5]); %pixel # for each nuclei
idx=unique(L);
[~,b]=max(n);
%idx(n>500)=b;
Lfilt=L;
Lfilt(ismember(Lfilt(:),idx(n>nucleisizethresh)))=b;
Lfiltrgb = label2rgb(Lfilt, 'jet', 'w', 'shuffle');
%figure
%imshow(Lfiltrgb)
%title('Nuclei(filtered)')

% switch background to 0;
Lfilt=Lfilt+1;
Lfilt(Lfilt==(b+1))=0;
%expand nuclei on each side by gradually dilate image to find cells
Lnew=Lfilt;
for i=1:cellsize
    se=strel('disk',i);
    L1=imdilate(Lfilt,se);
    Lnew(Lnew==0)=L1(Lnew==0);
end
Lnew1rgb = label2rgb(Lnew, 'jet', 'w', 'shuffle');
figure
imshow(Lnew1rgb)
title('Cells (size filtered)')

imnuc=Lfilt;
imcell=Lnew;


save('segmentation.mat','imnuc','imcell');
