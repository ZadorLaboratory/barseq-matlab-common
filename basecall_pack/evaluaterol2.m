function genes=evaluaterol2(match,ref,seqC,rol)
%compare rolonies to first cycle image for sequencing against genes
%run this after seqfindroloniesij and nricpbasecall in the aligned folder.
%This version takes rol of first cycle only and call on maxproj.

%assign genes to rolonies
ref1=char(ref);
d=pdist2(seqC,ref1(:,1:size(seqC,2)),'hamming')*size(seqC,2);
[mind,genes]=min(d,[],2);
genes(sum(d==repmat(mind,1,length(ref)),2)>1)=0;%filter out ambiguous ones.
genes(mind>(size(seqC,2)-match))=0; %filter out genes not matching enough bases.



%read the rgb image of the first sequencing cycle
olddir=cd('RGB/');
files=dir('*.tif');
im=imread(files(1).name);

cd(olddir);

%sort genes by abundance in ascending order so to plot the most abundant genes first.
c=zeros(length(ref),1);
for i=1:length(ref)
	c(i)=sum(genes==i);
end
[~,I]=sort(c,'descend');

%plot rolonies side by side with first cycle sequencing rgb image
figure('Position',[100 100 1000 700]);
subplot('Position',[0.05 0.1 0.4 0.8]); hold on;
for i=1:length(ref)
    scatter(rol(genes==I(i),2),rol(genes==I(i),1),1,'filled');
    set(gca,'xlim',[0 size(im,2)],'ylim',[0 size(im,1)]);
    axis off
end
%legend({'Slc30a3','Ptprk','Fam19a9','Serpine2','Gad1','Gad2','Slc17a7'});
set(gca,'ydir','reverse');
title(['Rol']);

subplot('Position', [0.55 0.1 0.4 0.8]);
imshow(im);
title(['Img']);
savefig('rol.fig');


%plot rolonies superimposed on the rgb first cycle sequencing image 
figure;
imshow(im); hold on;
for i=1:length(ref)
    scatter(rol(genes==I(i),2),rol(genes==I(i),1),'+');
    set(gca,'xlim',[0 size(im,2)],'ylim',[0 size(im,1)]);
    axis off
end
set(gca,'ydir','reverse');
savefig('rol1.fig');

