function validx=RFvalidate(x,idx)
%validation of clusters using random forest. Same way as Tasic et al., 2016
% pairwise validation of clusters using random forest to find the
% probability of a cell in a particular cluster.

Treenum=50;
clusters=unique(idx);
validx=zeros(size(x,1),size(clusters,1),size(clusters,1));


% do 10 times of 5x cross-validation
for m=1:10
    m
    Indices = crossvalind('Kfold', size(x,1), 5);
    validx=validx+xval1(x,idx,clusters,Indices,Treenum);
end
% fill the diagonal with 10
for i=1:size(clusters)
        validx(:,i,i)=10;
end


end

% for one round of xval
function validx=xval1(x,idx,clusters,Indices,Treenum)
validx=zeros(size(x,1),size(clusters,1),size(clusters,1));
for n=1:max(Indices)
    testx=x(Indices==n,:);
    trainx=x(Indices~=n,:);
    trainidx=idx(Indices~=n);
    n
    % train classifier for each pair ot clusters
    for i=1:(size(clusters,1)-1)
        for j=(i+1):size(clusters,1)
        % pick the two groups for validation
        subx=trainx(trainidx==clusters(i)|trainidx==clusters(j),:);
        subidx=trainidx(trainidx==clusters(i)|trainidx==clusters(j));
        Mdl=TreeBagger(Treenum,subx,subidx, 'MaxNumSplits',10);
        validx(Indices==n,i,j)=cellfun(@(y) str2double(y)==clusters(i),predict(Mdl,testx));
        validx(Indices==n,j,i)=cellfun(@(y) str2double(y)==clusters(j),predict(Mdl,testx));
        
        end
    end
end

end
    
        
        