function mmsegimgcrop(segim1,scorewhole1,cellid10x1,qualthresh)
%crop and concatenate 1st cycle gene sequence, barcode sequene, and
%segmentation image. run from the RGB folder.
if ~exist('qualthresh','var')
    qualthresh=0.8;
end


BCfiles=dir(fullfile('*BC*.tif'));
[BCfiles,~]=sort_nat({BCfiles.name});
genefiles=dir(fullfile('*hyb*.tif'));
[genefiles,~]=sort_nat({genefiles.name});

BCim=imread(BCfiles{1});
geneim=imread(genefiles{1});
geneim=padarray(geneim,[160,160],0,'both');
BCim=padarray(BCim,[160,160],0,'both');
segim2=padarray(segim1,[160,160],0,'both');

for i=1:length(cellid10x1)
    if min(scorewhole1(i,:))>qualthresh
    %find center of mass for cell
    [I1,I2]=ind2sub(size(segim1),find(segim1==cellid10x1(i)));
    
    cim=[repmat(imcrop(uint8((segim2==cellid10x1(i))*255),[mean(I2),mean(I1),320,320]),1,1,3),imcrop(BCim,[mean(I2),mean(I1),320,320]),imcrop(geneim,[mean(I2),mean(I1),320,320])];
    cim(:,[321 642],:)=255;
    imwrite(cim,[num2str(cellid10x1(i)),'_segcrop.tif']);
    
    end
end



