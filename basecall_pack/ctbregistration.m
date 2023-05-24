%%
file10x=dir('MAX*.tif');
file10x=sort_nat({file10x.name});
file40x=dir('stitched*.tif');
file40x=sort_nat({file40x.name});


%%
for i=1:numel(file40x)
    %%
    im1=imread(file10x{i},1);
    im1(:,:,2)=imread(file10x{i},2);
    finfo=imfinfo(file40x{i});
    im2=[];
    for n=1:numel(finfo)
        im2(:,:,n)=imread(file40x{i},n);
    end
    im1=imresize(imrotate(im1,180),4);
    %%
    

    %%
    normim1=im1(:,:,1);
    normim1=double((normim1-min(normim1(:))))/double(prctile(normim1(:)-min(normim1(:)),99.9));
    normim2=im2(:,:,6);
    normim2=double((normim2-min(normim2(:))))/double(prctile(normim2(:)-min(normim2(:)),99.9));
    [selectedMovingPoints,selectedFixedPoints] = cpselect(normim1,normim2,'Wait',true);
    tform=fitgeotrans(selectedMovingPoints,selectedFixedPoints,'nonreflectivesimilarity');
    
    rfixed=imref2d(size(im2(:,:,6)));
    im1t=imwarp(im1,tform,'OutputView',rfixed);
    
    %%
    normim1t=imwarp(normim1,tform,'OutputView',rfixed);
        figure;
    subplot(2,2,1)
    imshow(normim1t);
    title(file10x{i});
    
    subplot(2,2,2)
    imshow(normim2);
    title(file40x{i});
    
    subplot(2,2,3)
    imshow(normim1t);
    hold on;
    plot([3250 3250],[3200 3300],'r','LineWidth',0.5);
    plot([3200 3300],[3250 3250],'r','LineWidth',0.5);
    
    set(gca,'xlim',[3000 3500],'ylim',[3000 3500]);
    title(file10x{i});
    
    subplot(2,2,4)
    imshow(normim2);
    hold on;
    plot([3250 3250],[3200 3300],'r','LineWidth',0.5);
    plot([3200 3300],[3250 3250],'r','LineWidth',0.5);
    set(gca,'xlim',[3000 3500],'ylim',[3000 3500]);
    title(file40x{i});
    
    
    
    %%
    im2t=im2;
    im2t(:,:,7)=im2t(:,:,4);
    im2t(:,:,4)=im1t(:,:,2);
    im2t=uint16(im2t);
    for n=1:size(im2t,3)
        if n==1
            imwrite(im2t(:,:,n),['tformed',file40x{i}]);
        else
            imwrite(im2t(:,:,n),['tformed',file40x{i}],'WriteMode','Append');
        end
    end
end