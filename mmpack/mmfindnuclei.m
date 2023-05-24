function [imnuc1,imcell1,cellpos1,celllist1]=mmfindnuclei(nucleifile,nucleich,genefile,genech,hybfile,hybch,blackpoint,bgnthresh,maxcellsize,expandsize)
%find nuclei as markers and segment cells using either nuclei or rolony density using
%marker-based watershed.Output nuc and cell with nuclei/cell labels. background is labeled 0. 
%This version is optimized for 40x (scale all filters accordingly for other
%mags)
%   nucleifile: partial file name of the nuclei stain files
%   genefile: partial filename for sequencing files (this has cell-body
%   background)
%   hybfile: partial file name for sequential round files. (this only has
%   rolonies)
%   nucleich: nuclear signal ch
%   genech: sequencing signal ch
%   hybch: sequential signal ch to be included (i.e. exclude the ones already in gene ch)
%   blackpoint: blackpoint of the image for initial background subtraction.
%   Can be set pretty aggresively
%   bgnthresh: additional threshold for background after blackpoint
%   subtraction, for finding backgrounds.
%   nucleisizethresh: max size of nuclei for filtering out multiple nuclei,
%   in # of pixels.
%   cellsize: the number of pixels to expand from each side of nuclei.

%%
%find nuclei as markers for watershed
nucfile=dir(['*',nucleifile,'*.tif']);
nucfile=sort_nat({nucfile.name});

%read image minus blackpoint
im=imread(nucfile{1},nucleich);
im1=max(im-blackpoint,0);
im11=imgaussfilt(max(im-blackpoint,0),10); %reduces granularity within a nucleus, better than pill box


%find foreground (nuclei) by opening + erosion
se=strel('disk',12); %This determines how far the nuclei must be a part or they'll get merged
Io=imopen(im11,se);
%figure; imshow(min(Io*2,1)),title('Opening');
se2 = strel('disk',5); %these are not very important for tuning which to call nuclei, but determines the spot size for each nuclei marker.
se21=strel('disk',7);%these are not very important for tuning which to call nuclei, but determines the spot size for each nuclei marker.
fgm1 = imerode(imclose(imregionalmax(Io),se2),se21);
%figure;
%subplot(1,3,1);
%imshowpair(min(double(im11)/max(double(im11(:)))*3,1),fgm1);


%%
%find background by thresholding
%bgn=imbinarize(im1,0.1); %might need manual adjustment here
%figure;imshow(bgn);
im2=medfilt2(double(im1))-imopen(double(im1),strel('ball',100,100));
bgn=(im2-bgnthresh)>0;
%thin background
D = bwdist(bgn);
DL = watershed(D);
bgm = DL == 0;
% figure; 
% %imshow(bgm); 
% %subplot(1,3,2);
% imshowpair(im2,bgm+fgm1); title('Watershed ridge lines (bgm)')


%%
%marker-controlled watershed on the gradient magnitude

%do a first watershed on nuclei to find the whole nuclei
nuc=watershed(imimposemin(imgradient(im), bgm | fgm1));
 n=histcounts(nuc,[-0.5;unique(nuc)+0.5]); %pixel # for each nuclei
[~,b]=max(n);
nuc(nuc==(b-1))=0;
nuc(nuc>0)=1;

%watershed with rolony density. This is not very satisfactory probably
%because rolony density is not high enough.
% if isempty(rol1)
%     L=watershed(imimposemin(imgradient(im1), bgm | fgm1));
% else
%     %make rolony density
%     roldensity=zeros(size(im1));
%     for n=1:length(rol1)
%         roldensity(rol1(n,2),rol1(n,1))=roldensity(rol1(n,2),rol1(n,1))+1;
%     end
%     roldensity=imgaussfilt(conv2(roldensity,fspecial('disk',40),'same'),10);
%     %roldensity=imgaussfilt(roldensity,30);
%     L=watershed(imimposemin(imgradient(im), bgm | fgm1));
% end

%Watershed with aligned sequencing images. These images have
%cell body backgrounds may help when rolonies are not dense enough.
genefile=dir(['*',genefile,'*.tif']);
genefile=sort_nat({genefile.name});

geneim=imread(genefile{1},genech(1));
for n=2:length(genech)
    geneim(:,:,n)=imread(genefile{1},genech(n));
end

if ~isempty(hybch)
    hybfile=dir(['*',hybfile,'*.tif']);
    hybfile=sort_nat({hybfile.name});
    hybim=imread(hybfile{1},hybch(1));
    for n=2:length(hybch)
        hybim=imread(hybfile{1},hybch(n));
    end
    sumim=cat(3,geneim,hybim);
else
    sumim=geneim;
end

for n=1:size(sumim,3)
    sumim(:,:,n)=imgaussfilt(sumim(:,:,n),10);%blur rolonies
end
sumim=sum(double(sumim),3);
%figure;imagesc(sumim);


%figure;imagesc(imgradient(sumim));

L=watershed(imimposemin(imgradient(sumim), bgm | nuc)); 


%%
%[fgm1y,fgm1x]=ind2sub(size(im1),find(fgm1));
% figure;imagesc(roldensity);hold on;
% imagesc(im1);
% 
% scatter(rol1(:,1),rol1(:,2),'.g');
% scatter(fgm1x,fgm1y,'.r');

Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
% figure
% imshow(Lrgb);hold on;%scatter(rol1(:,1),rol1(:,2),0.5,'g','.');scatter(fgm1x,fgm1y,'.r');
% % title('Nuclei(unfiltered) image gauss=10 nuc')

%%
% filter out cells that are too large (likely multiple cells merged).
n=histcounts(L,[-0.5;unique(L)+0.5]); %pixel # for each nuclei
idx=unique(L);
[~,b]=max(n);
%idx(n>500)=b;
Lfilt=L;
Lfilt(ismember(Lfilt(:),idx(n>maxcellsize)))=b;
%Lfiltrgb = label2rgb(Lfilt, 'jet', 'w', 'shuffle');
%figure
%imshow(Lfiltrgb)
%title('Nuclei(filtered)')

% switch background to 0;
Lfilt=Lfilt+1;
Lfilt(Lfilt==(b+1))=0;
%expand nuclei on each side by gradually dilate image to find cells
Lnew=Lfilt;
for i=1:expandsize
    se=strel('disk',i);
    L1=imdilate(Lfilt,se);
    Lnew(Lnew==0)=L1(Lnew==0);
end
Lnew1rgb = label2rgb(Lnew, 'jet', 'w', 'shuffle');
% figure;
% %subplot(1,3,3)
% imshow(Lnew1rgb)
% title('Cells (size filtered)')

imnuc1=nuc;
imcell1=Lnew;
%% calculate cell positions
celllist1=unique(Lnew);
celllist1(celllist1==0)=[];
cellpos1=zeros(length(celllist1),2);
for i=1:length(celllist1)
    [y,x]=ind2sub(size(Lnew),find(Lnew==celllist1(i)));
    cellpos1(i,:)=[median(x),median(y)];
end





%%

save('segmentation.mat','imnuc1','imcell1','cellpos1','celllist1');
