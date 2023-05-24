function genes=evaluaterol(match,Nmatch)
%compare rolonies to first cycle image for sequencing against genes
%run this after seqfindroloniesij and nricpbasecall in the aligned folder
load rolonies
load seq1;
load ref;

%assign genes to rolonies
genes=zeros(length(seqC),1);
for i=1:length(ref)
    ref1=repmat(ref{i},length(seqC),1);
    I=sum(ref1(:,1:size(seqC,2))==seqC,2);
    genes(I>=match&(I+sum(seqC=='N',2))>=(Nmatch+match))=i;
end

%read the rgb image of the first sequencing cycle
olddir=cd('../RGB/');
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
    scatter(rol{1}(genes==I(i),2),rol{1}(genes==I(i),1),1,'filled');
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
    scatter(rol{1}(genes==I(i),2),rol{1}(genes==I(i),1),'+');
    set(gca,'xlim',[0 size(im,2)],'ylim',[0 size(im,1)]);
    axis off
end
set(gca,'ydir','reverse');
savefig('rol1.fig');

