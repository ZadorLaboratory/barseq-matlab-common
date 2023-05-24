function nuccellid1=mmsegimgcrop2(bcsegim1,segim1,scorewhole1,cellid10x1,qualthresh)
%crop and concatenate 1st cycle gene sequence, barcode sequene, and
%segmentation image from both BC segmentatino and nuclei segmentation. run from the RGB folder.
%also calculate the correspondance from bcsegim1 to segim1
if ~exist('qualthresh','var')
    qualthresh=0.8;
end
%find correponding number in segim1 from bcsegim1
nuccellid1=zeros(length(cellid10x1),1);
for i=1:length(cellid10x1)
    u=unique([0;segim1(bcsegim1(:)==cellid10x1(i))]);%add 0 in case all pixels overlap with segim1.
    %if no nuclei segmentatino found, then nuccellid=0. otherwise find the
    %nuc segmentation with the largest overlap
    if sum(u~=0)>0
        c = histcounts(segim1(bcsegim1==cellid10x1(i)),u+0.5);
        [~,I]=max(c);
        nuccellid1(i)=u(I+1);
    else
        nuccellid1(i)=0;
    end
    
end



BCfiles=dir(fullfile('*BC*.tif'));
[BCfiles,~]=sort_nat({BCfiles.name});
genefiles=dir(fullfile('*hyb*.tif'));
[genefiles,~]=sort_nat({genefiles.name});

BCim=imread(BCfiles{1});
geneim=imread(genefiles{1});
geneim=padarray(geneim,[160,160],0,'both');
BCim=padarray(BCim,[160,160],0,'both');
bcsegim2=padarray(bcsegim1,[160,160],0,'both');
segim2=padarray(segim1,[160,160],0,'both');

for i=1:length(cellid10x1)
    if min(scorewhole1(i,:))>qualthresh
    %find center of mass for cell
    [I1,I2]=ind2sub(size(bcsegim1),find(bcsegim1==cellid10x1(i)));
    
    cim=[repmat(imcrop(uint8((bcsegim2==cellid10x1(i))*255),[mean(I2),mean(I1),320,320]),1,1,3), ...
        repmat(imcrop(uint8((segim2==nuccellid1(i))*255),[mean(I2),mean(I1),320,320]),1,1,3), ...
        imcrop(BCim,[mean(I2),mean(I1),320,320]), ...
        imcrop(geneim,[mean(I2),mean(I1),320,320])];
    cim(:,[321 642 963],:)=255;
    imwrite(cim,[num2str(cellid10x1(i)),'_segcropnuclei.tif']);
    
    end
end



