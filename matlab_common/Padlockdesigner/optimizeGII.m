function [idx,mind]=optimizeGII(oligpool,curridx,seqlength,numoligs,replicates)
%given a pool of oligos and a starting pool of oligo indices, find numoligs
%number of additional oligos to maximize oligo to oligo distance in the
%expected seqlength.

oligpool=char(oligpool);
mind=zeros(numoligs,size(oligpool,2)-seqlength+1,replicates);
idx0=zeros(numoligs+length(curridx),replicates);
parfor m=1:replicates
    idx=reshape(curridx,numel(curridx),1);
    mind1=zeros(numoligs,size(oligpool,2)-seqlength+1);
    for i=1:numoligs
        currolig=oligpool(idx,:);
        %idx1=1:length(oligpool);
        d=min(pdist2(double(oligpool),double(currolig),'hamming')*size(oligpool,2),[],2);
        dshort=min(pdist2(double(oligpool(:,1:seqlength)),double(currolig(:,1:seqlength)),'hamming')*seqlength,[],2);
        
        %for n=0:size(oligpool,2)-seqlength
        %    d=min(pdist2(oligpool(:,1:seqlength),currolig(:,1:seqlength),'hamming')*seqlength,[],2);
        %    idx1=idx1(d(idx1)==max(d(idx1)));
        %end
        idx1=find(d==max(d)&dshort==max(dshort(d==max(d))));
        idx=[idx;idx1(randperm(length(idx1),1))];
        for n=1:size(oligpool,2)-seqlength+1
            mind1(i,n)=min(pdist(double(oligpool(idx,1:seqlength+n-1)),'hamming')*(seqlength+n-1));
        end
        
    end
    mind(:,:,m)=mind1;
   idx0(:,m)=idx;
end
[~,I]=max(sum(sum(mind(:,1,:))));
mind=mind(:,:,I);
idx=idx0(:,I);
