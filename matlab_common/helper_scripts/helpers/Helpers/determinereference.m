function p=determinereference(input,reference,n,z,distance)


if z==1
    reference=zscore(reference');
    reference=reference';
end


dist=[];
for repeats=1:n;
    %sample from referencebarcodes
    sample=reference(randi(size(reference,1),size(input,1),1),:);
    %find closest neighbours
    [idx,d]=knnsearch(sample,input,'K',1,'Distance',distance);
    %keep track of distances
    dist=[dist;d];
end
p=hist(dist,0:0.001:10)/sum(hist(dist,0:0.001:10));
