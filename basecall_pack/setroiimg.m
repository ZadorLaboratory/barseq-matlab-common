function [region,filtpxl]=setroiimg(idx,region)
%draw an ROI and find all pixels within the region.idx is the index of 
%the channel used to draw ROI.

    files=dir('*.tif');
    files=sort_nat({files.name});
    
if ~exist('region','var')
    %draw area
    for n=1:length(files)
        im=imread(files{n},idx);
        im1=sort(im(:),'descend');
        figure; imshow(im,[0 im1(floor(length(im1)/10))],'InitialMagnification',50);title([files{n},', set points defining convex hull, then press enter']);
        hold on;
        topx=[];topy=[];
        while 1
            [x,y]=myginput(1,'cross');
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
        close all;
        region{n}=[topx(convhull(topx,topy)),topy(convhull(topx,topy))];
    end
end

%find all pixels within convex hull
for n=1:length(files)
    %find all pixels within convex hall
    im=imread(files{n},idx);
    [py,px]=ind2sub(size(im),1:numel(im));
    in = inpolygon(px,py,region{n}(:,1),region{n}(:,2));
    py=py(in);px=px(in);
    filtpxl{n}=[px',py'];
end

end
