function rol=findrolonies(im,thresh)
%find rolonies in image im. The coordinates of all indices are in rol.
%for i=1:4
%    im1(:,:,i)=imread('max.tif',i);
%    im1(:,:,i)=medfilt2(im1(:,:,i));
%end
%im=im1(300:600,300:600,:);

%thresh=20;
bgn=0;
%find all maxima separated by thresh
peaks=im>=imdilate(im, [1 1 1; 1 0 1; 1 1 1]); %peak positions

filteredpeaks=zeros(size(im));
for i=1:4
    %sort peaks by values
    peakvalues=peaks(:,:,i).*im(:,:,i);
    [~,I]=sort(reshape(peakvalues,[],1),'descend');
    [idx1,idx2]=ind2sub(size(peaks(:,:,i)),I(1:sum(sum(peakvalues>(bgn+thresh)))));
    I2=zeros(size(peaks(:,:,i)))==1;
    for n=1:length(idx1)
        if ~I2(idx1(n),idx2(n))
            %flood-fill image using thresh
            %tic
            I1=im(:,:,i)<(peakvalues(I(n))-thresh);
            %tic
            I2=I2|(imfill(I1,[idx1(n),idx2(n)])&~I1);
            %toc
            filteredpeaks(idx1(n),idx2(n),i)=1;
            %toc
        end
        
    end
    
end

%check rolonies
%for i=1:4
%figure;imshow(im(:,:,i),[0 300]);
%hold on;
%[idx1,idx2]=ind2sub(size(im(:,:,i)),find(filteredpeaks(:,:,i)));
%scatter(idx2,idx1,'+');
%end

%pool rolonies, combine rolonies right next to each other 
allpeaks=sum(filteredpeaks,3)>0;
%[idx1,idx2]=ind2sub(size(allpeaks),find(allpeaks));
%figure;scatter(idx2,idx1);

allpeaks=allpeaks&~imdilate(allpeaks,[1 1 1;1 0 0; 0 0 0]);
[idx1,idx2]=ind2sub(size(allpeaks),find(allpeaks));
%figure;scatter(idx2,idx1);
rol=[idx1,idx2];
end
    
        
        
        
    
    
    
    


