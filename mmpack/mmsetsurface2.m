function contours=mmsetsurface2(poslist,idx)
%select points delineating the top and bottom of the cortex. idx indicates
%the channel to display. This version operates on folders containing  
%sequencing images of individual slices.
contours=struct;
for n=1:length(poslist)
    cd([poslist{n},'/original']);
    files=dir('*.tif');
    files=sort_nat({files.name});
    im=imread(files{1},idx);
    im1=sort(im(:),'descend');
    figure; imshow(im,[0 im1(floor(length(im1)/10))],'InitialMagnification',50);title([poslist{n},', set top, press enter; then set bottom']);
    hold on;
    topx=[];topy=[];
    while 1
        [x,y]=myginput(1,'topl');
        if numel(x)==1
            topx=[topx;x];
            topy=[topy;y];
            if numel(topx)>=2
                plot(topx(end-1:end),topy(end-1:end),'r','LineWidth',2);
            else
                scatter(x,y,10,'r','filled');
            end
        else
            break
        end
    end
    bottomx=[];bottomy=[];
    while 1
        [x,y]=myginput(1,'cross');
        if numel(x)==1
            bottomx=[bottomx;x];
            bottomy=[bottomy;y];
            if numel(bottomx)>=2
                plot(bottomx(end-1:end),bottomy(end-1:end),'g','LineWidth',2);
            else
                scatter(x,y,10,'g','filled');
            end
        else
            break
        end
    end
   
    close all;
    contours(n).slice=poslist{n};
    contours(n).top=[topx,topy];
    contours(n).bottom=[bottomx,bottomy];
    cd ../..
end
save('contours.mat','contours');


