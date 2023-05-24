function [seq,seqC,score,sig,ftlroi]=mmbasecallpoints(slicenames,lthresh)
%find corresponding points of lroi in the 10x sequencing files, basecall
%those positions. Used for calling 10x cell body barcode sequencing images.
addpath('C:\Fiji.app\scripts');
javaaddpath 'C:\Program Files\MATLAB\R2018b\java\mij.jar'
Miji(false);
sig=cell(length(slicenames),1);
score=cell(length(slicenames),1);
seq=cell(length(slicenames),1);
seqC=cell(length(slicenames),1);
lroi=cell(length(slicenames),1);
for i=1:length(slicenames)
    %fine-tune lroi positions
    rootdir=cd([slicenames{i},'/aligned/']);
    seqfiles=dir(fullfile('*seq*.tif'));
    [seqfiles,~]=sort_nat({seqfiles.name});
    
    lim=imread(seqfiles{1},1);
    for n=2:4
        lim(:,:,n)=imread(seqfiles{1},n);
    end
    
    currrol=[];
    %find rolonies after gaussian blur with sigma=1
    for n=1:4
        MIJ.createImage(imgaussfilt(lim(:,:,n),1));
        MIJ.run('Find Maxima...', ['noise=',num2str(lthresh),' output=List']);
        currrol=[currrol;MIJ.getResultsTable];
        %cleanup
        MIJ.run('Clear Results');
        MIJ.closeAllWindows
    end
    
    lpeaks=zeros(size(lim,1),size(lim,2));
    lpeaks(sub2ind([size(lim,1),size(lim,2)],currrol(:,2)+1,currrol(:,1)+1))=1;
    lpeaks=lpeaks&~imdilate(lpeaks,triu(ones(15))-diag([ones(1,8),zeros(1,7)]));
    [lidxy,lidxx]=find(lpeaks);
    lroi{i}=[lidxx,lidxy];
    
    %readout signals
    currsig=ones(length(lroi{i}),length(seqfiles),4);
    
    for m=1:length(seqfiles) %for each cycle, read images and convolute by 7x7
        lim=conv2(imread(seqfiles{m}),ones(7,7),'same');
        
        for n=2:4
            lim(:,:,n)=conv2(imread(seqfiles{m},n),ones(7,7),'same');
        end
        
        for n=1:length(lroi{i}) %for each rolony, readout signal
            currsig(n,m,:)=lim(lroi{i}(n,2),lroi{i}(n,1),:);
            
        end
        
    end
    %calculate score
    currscore=max(currsig,[],3)./sqrt(sum(currsig.^2,3));
    %convert to base
    [~,currseq]=max(currsig,[],3);
    currseqC=char(currseq);
    currseqC(currseq==1)='G';
    currseqC(currseq==2)='T';
    currseqC(currseq==3)='A';
    currseqC(currseq==4)='C';
    %assign to cells
    currroi=lroi{i};
    save('basecalls.mat','currseq','currseqC','currscore','currsig','currroi');
    sig{i}=currsig;
    seq{i}=currseq;
    score{i}=currscore;
    seqC{i}=currseqC;
    cd(rootdir);
end
ftlroi=lroi;
save('allbasecalls.mat','seq','seqC','score','sig','ftlroi');
MIJ.exit;
end