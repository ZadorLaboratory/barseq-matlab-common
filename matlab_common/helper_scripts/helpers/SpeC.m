function [idx,pval]=SpeC(x,dist)
% kmean cluster x into two groups, check significance. x rows are
% observations and columns are features.
%if size(x,1)>269
            W=squareform(pdist(x,dist));
            W=max(max(W))-W;
            
            [idx, ~, ~] = SpectralClustering(W, 2,3);
            idx=full(idx);
            idx=idx(:,1)*1+idx(:,2)*2;    
    paramstruct.iCovEst=2;
    paramstruct.vclass=idx';
    %normx=x+(x>0).*repmat(repmat(max(max(x)),size(x,1),1)-max(x,[],2),1,size(x,2)); %normalize x
    %normx=x./repmat(sum(x,2),1,size(x,2));    %normalize x
    if sum(idx==1)<64||sum(idx==2)<64 % if either cluster is smaller than a predefinned value, remove the clusters
        
%else
        idx=ones(size(x,1),1);
        pval=1;
        C=0;
        z=0;
    else
        [pval,~]=SigClustSM(x',paramstruct);
    end
    
end


