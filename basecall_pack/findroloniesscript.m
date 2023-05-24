for i=1:4
    im1(:,:,i)=imread('max.tif',i);
    im1(:,:,i)=medfilt2(im1(:,:,i));
end
im=im1(300:600,300:600,:);


thresh=20;
bgn=100;
peaks=im>=imdilate(im, [1 1 1; 1 0 1; 1 1 1]); %peak positions
filteredpeaks=zeros(size(im));
for i=1:4
    %sort peaks by values
    
    peakvalues=uint16(peaks(:,:,i)).*im(:,:,i);
    [~,I]=sort(reshape(peakvalues,[],1),'descend');
    [idx1,idx2]=ind2sub(size(peaks(:,:,i)),I(1:sum(sum(peakvalues>(bgn+thresh)))));
    I2=zeros(size(peaks(:,:,i)),'logical');
    for n=1:1000
        if ~I2(idx1(n),idx2(n))
            %flood-fill image using thresh
            I1=im(:,:,i)<(peakvalues(I(n))-thresh);
            I2=I2|(imfill(I1,[idx1(n),idx2(n)])&~I1);
            filteredpeaks(idx1(n),idx2(n),i)=1;
        end
        
    end
    
end

%check rolonies
for i=1:4
figure;imshow(im(:,:,i),[100 300]);
hold on;
[idx1,idx2]=ind2sub(size(im(:,:,i)),find(filteredpeaks(:,:,i)));
scatter(idx2,idx1,'+');
end

%pool rolonies, combine rolonies right next to each other 
allpeaks=sum(filteredpeaks,3)>0;
[idx1,idx2]=ind2sub(size(allpeaks),find(allpeaks));
figure;scatter(idx2,idx1);

allpeaks=allpeaks&~imdilate(allpeaks,[1 1 1;1 0 0; 0 0 0]);
[idx1,idx2]=ind2sub(size(allpeaks),find(allpeaks));
figure;scatter(idx2,idx1);

%transform
allpeaks1=imtranslate(allpeaks,[15 25]);


%plot original peaks vs transformed peaks
[idx1,idx2]=ind2sub(size(allpeaks),find(allpeaks));
figure;scatter(idx2,idx1,'x');
hold on;
[idx3,idx4]=ind2sub(size(allpeaks1),find(allpeaks1));
scatter(idx4,idx3,'+');

%drop rolonies
idx3(randperm(length(idx3),30))=0;
idx4(idx3==0)=[];
idx3(idx3==0)=[];

figure;scatter(idx2,idx1,'x');
hold on;
scatter(idx4,idx3,'+');

%ICP
[tr,tt]=icp([idx1';idx2';ones(1,length(idx1))],[idx3';idx4';ones(1,length(idx3))],100);
figure;figure;scatter(idx2,idx1,'x');
hold on;scatter(idx4,idx3,'+');


save('rolonies.mat','filteredpeaks');

    
        
        
        
    
    
    
    


