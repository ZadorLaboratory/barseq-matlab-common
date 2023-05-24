function []=findclosestrandommatch(input,BCnumber,BCnumbertarget,reference,z,distance)


if z==1
    reference=zscore(reference');
    reference=reference';
end


dist=[];
% figure;
for targets=1:length(BCnumber);
    figure;
for repeats=1:5;
    %sample from referencebarcodes
    sample=reference(randi(size(reference,1),size(input,1),1),:);
    %find closest neighbours
    [idx,d]=knnsearch(sample,input(BCnumber(targets),:),'K',1,'Distance',distance);
    %keep track of distances
    dist=[dist;d];
    
%    subplot(length(BCnumber),5,(targets-1)*5+repeats);
    plot(sample(idx,:),'black');
    hold on;
    
end
plot(input(BCnumber(targets),:),'red');
plot(input(BCnumbertarget(targets),:),'blue');
xlim([0 23])
end
