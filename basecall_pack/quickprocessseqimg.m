imfiles=dir('*.tif');
for i=1:length(imfiles)
    info=imfinfo(imfiles(i).name);
    for n=1:length(info)
        im(:,:,n)=medfilt2(imread(imfiles(i).name,n));
    end
    
    im=double(im-imopen(im,strel('ball',20,20)));
    im(:,:,1)=max(im(:,:,1)-im(:,:,2)*0.74,0);
    im(:,:,4)=max(im(:,:,4)-im(:,:,3)*0.42,0);
    imwrite(uint16(im(:,:,1)),['fixed',imfiles(i).name]);
    for n=1:length(info)
        imwrite(uint16(im(:,:,n)),['fixed',imfiles(i).name],'WriteMode','Append');
    end
    
    imrgb=repmat(im(:,:,1),1,1,3);
    imrgb(:,:,1)=sum(im(:,:,2:4),3);
    imrgb(:,:,2)=sum(im(:,:,[1 2 4]),3);
    imrgb(:,:,3)=sum(im(:,:,[1 3 4]),3);
    imrgb=min(imrgb./repmat(max(max(imrgb,[],1),[],2),size(imrgb,1),size(imrgb,2),1),1);
    imwrite(imrgb,['RGB',imfiles(i).name]);
end


    