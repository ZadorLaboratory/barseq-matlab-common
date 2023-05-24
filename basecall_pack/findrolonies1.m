function rol=findrolonies1(im,thresh)
%find rolonies in image im. The coordinates of all indices are in rol.
%using imhmax

filteredpeaks=zeros(size(im));
for i=1:4
    filteredpeaks(:,:,i)=imhmax(im(:,:,i),thresh);
end

%check rolonies
for i=1:4
figure;imshow(im(:,:,i),[0 300]);
hold on;
[idx1,idx2]=ind2sub(size(im(:,:,i)),find(filteredpeaks(:,:,i)));
scatter(idx2,idx1,'+');
end

%pool rolonies, combine rolonies right next to each other 
allpeaks=sum(filteredpeaks,3)>0;
%[idx1,idx2]=ind2sub(size(allpeaks),find(allpeaks));
%figure;scatter(idx2,idx1);

allpeaks=allpeaks&~imdilate(allpeaks,[1 1 1;1 0 0; 0 0 0]);
[idx1,idx2]=ind2sub(size(allpeaks),find(allpeaks));
%figure;scatter(idx2,idx1);
rol=[idx1,idx2];
end
    
        
        
        
    
    
    
    


