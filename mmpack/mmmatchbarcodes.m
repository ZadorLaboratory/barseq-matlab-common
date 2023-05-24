function [B,Bnorm,proj,projraw]=mmmatchbarcodes(barcodematrix,refbarcodes,spikes,seqCwhole,slicenames)
% given barcodematrix, individualdatawithdistancealigned, spikes, produce
% B, CoPNnorm, AsPNnorm,filtaligneddepth,fisseqCoPNnorm,fisseqAsPNnorm.
% slicenames is not used.

%filter MAPseq barcodes
lowercutoff=1;
uppercutoff=10000;

B=barcodematrix(max(barcodematrix,[],2)>lowercutoff & max(barcodematrix,[],2)<uppercutoff,:);
seq=refbarcodes(max(barcodematrix,[],2)>lowercutoff & max(barcodematrix,[],2)<uppercutoff,:);
x=[];for i=1:size(barcodematrix,2);x(i)=length(spikes(i).counts2u);end

Bnorm=B./repmat(x,length(B),1);


%in situ barcode length
for i=1:length(seqCwhole)
    if ~isempty(seqCwhole{i})
        bclength=size(seqCwhole{i},2);
        break;
    end
end


% examine min hamming distance and calculate false positive rate using 15 bases and allow one mismatch
seq15=seq(:,1:bclength);
% what is the hamming distance profile for these barcodes?
d15=min(squareform(pdist(seq15,'hamming'))+eye(size(seq15,1)),[],2)*bclength;

figure;
histogram(d15,-0.5:9.5,'FaceColor',[0 0.5 0.5]);xlabel('Min hamming distance');ylabel('Counts');

%match in situ reads to projections
proj=cell(length(seqCwhole),1);projraw=proj;
for i=1:length(seqCwhole)
    if ~isempty(seqCwhole{i})
        seqC1=int8(seqCwhole{i});
        d=pdist2(seqC1,seq15,'hamming')*bclength;
        [mind,idx]=min(d,[],2);
        c=sum(d==repmat(mind,1,size(d,2)),2);
        proj{i}=zeros(size(d,1),size(B,2));
        proj{i}(mind<=1&c==1,:)=Bnorm(idx(mind<=1&c==1),:);
        projraw{i}=zeros(size(d,1),size(B,2));
        projraw{i}(mind<=1&c==1,:)=B(idx(mind<=1&c==1),:);
    else
        proj{i}=[];projraw{i}=[];
    end
end

