function countfish(hybthresh1,fixbleedthrough,ctbstrength)
%manually select cells with and without CTB labeling, manually draw cell contours then count fish
%signals in the cells. hybthresh1 is a 1x4 vector with the prominence
%threshold for identifying local maxima in fish.

if ~exist('ctbstrength','var')
    ctbstrength=99.9;
end


%check if images have been processed before
folders=dir('processed*');
if isempty(folders)
    if ~exist('fixbleedthrough','var')
        fixbleedthrough=1;
    end

    %%
    files=dir('*.tif');
    files=sort_nat({files.name});

    %segmentation produced by cellpose

    for i=1:numel(files)
        iminfo=imfinfo(files{i});
        im=zeros(iminfo(1).Height,iminfo(1).Width,numel(iminfo));
        parfor n=1:size(im,3)
            im(:,:,n)=imtophat(medfilt2(imread(files{i},n)),offsetstrel('ball',200,200));

        end
        if fixbleedthrough~=0
            im(:,:,4)=im(:,:,4)-0.57*im(:,:,2);%fix bleedthrough to CTB.
        end

        mkdir(['processed',files{i}(1:end-4)]);
        imwrite(uint16(im(:,:,1)),['processed',files{i}(1:end-4),'/processed',files{i}]);
        for n=2:size(im,3)
            imwrite(uint16(im(:,:,n)),['processed',files{i}(1:end-4),'/processed',files{i}],'WriteMode','Append');
        end
    end
end

%% draw cell positions and contours

%the following assumes that the channels are Slc30a3, Cdh, Slc17a7,
%CTB, DAPI, DIC

%make slc17a7+CTB image

folders=dir('processed*');
folders=sort_nat({folders.name});

for o=1:numel(folders)
    %%
    cd(folders{o});
    f=dir('processed*.tif');
    f={f.name};
    iminfo=imfinfo(f{1});
    im=zeros(iminfo(1).Height,iminfo(1).Width,numel(iminfo));
    parfor n=1:size(im,3)
        im(:,:,n)=imread(f{1},n);
    end
    figure;imshow(cat(3,im(:,:,4)./prctile(im(:,:,4),ctbstrength,'all'),im(:,:,3)./prctile(im(:,:,3),99.99,'all'),zeros(size(im(:,:,1)))));
    title('Please click on cells with CTB, then press enter.');
    hold on;
    
    ctbcells=[];
    while 1
        [x,y]=myginput(1,'cross');
        if isempty(x)
            break
        else
            ctbcells=[ctbcells;x,y];
            scatter(x,y,10,[0 1 1],'filled');
        end
    end
    title('Please click on cells WITHOUT CTB, then press enter.');
    nonctbcells=[];
    while 1
        [x,y]=myginput(1,'crosshair');
        if isempty(x)
            break
        else
            nonctbcells=[nonctbcells;x,y];
            scatter(x,y,10,[1 0 1],'filled');
        end
    end
    
    %%
    title('Please mark top of the cortex, then press enter.')
    top=[];
    while 1
        [x,y]=myginput(1,'cross');
        if isempty(x)
            break
        else
            if isempty(top)
                scatter(x,y,3,[1 0 0],'filled');
            else
                plot([top(end,1),x],[top(end,2),y],'-r','Linewidth',2);
            end
            top=[top;x,y];
        end
    end
    
    title('Please mark bottom of the cortex, then press enter.')
    bottom=[];
    while 1
        [x,y]=myginput(1,'cross');
        if isempty(x)
            break
        else
            if isempty(bottom)
                scatter(x,y,3,[0 1 0],'filled');
            else
                plot([bottom(end,1),x],[bottom(end,2),y],'-r','Linewidth',2);
            end
            bottom=[bottom;x,y];
        end
    end
    close;
    %calculate depths
    contours.top=top;
    contours.bottom=bottom;
    pixelsize=6.5/40;
    cellpos=[ctbcells;nonctbcells];
    [depths,angles]=mmcalculatedepth(contours,{cellpos},pixelsize);
    

    ctb=[ones(size(ctbcells,1),1);zeros(size(nonctbcells,1),1)];
    save('cellpos.mat','cellpos','ctb','top','bottom','depths','angles');
    
    %% make crops
    im1=padarray(im,[100,100,0],0,'both');
    im2=im1./prctile(im,99.9,[1 2]);
    
    ctbcells=uint16(ctbcells);
    for n=1:size(ctbcells,1)
        cim1=im1(ctbcells(n,2):ctbcells(n,2)+200,ctbcells(n,1):ctbcells(n,1)+200,:);
        cim2=im2(ctbcells(n,2):ctbcells(n,2)+200,ctbcells(n,1):ctbcells(n,1)+200,:);
        cim2=reshape(cim2,size(cim2,1),[],1);
        cim2(:,rem(1:size(cim2,2),200)==0)=1;
        imwrite(uint16(cim1(:,:,1)),['origcropctb',num2str(n),'.tif']);
        for i=2:size(cim1,3)
            imwrite(uint16(cim1(:,:,i)),['origcropctb',num2str(n),'.tif'],'WriteMode','Append');
        end
        imwrite(cim2,['normcropctb',num2str(n),'.tif']);
    end
    
    nonctbcells=uint16(nonctbcells);
    for n=1:size(nonctbcells,1)
        cim1=im1(nonctbcells(n,2):nonctbcells(n,2)+200,nonctbcells(n,1):nonctbcells(n,1)+200,:);
        cim2=im2(nonctbcells(n,2):nonctbcells(n,2)+200,nonctbcells(n,1):nonctbcells(n,1)+200,:);
        cim2=reshape(cim2,size(cim2,1),[],1);
        cim2(:,rem(1:size(cim2,2),200)==0)=1;
        imwrite(uint16(cim1(:,:,1)),['origcropnonctb',num2str(n),'.tif']);
        for i=2:size(cim1,3)
            imwrite(uint16(cim1(:,:,i)),['origcropnonctb',num2str(n),'.tif'],'WriteMode','Append');
        end
        imwrite(cim2,['normcropnonctb',num2str(n),'.tif']);
    end
    
    %% draw cells
    files=dir('normcrop*.tif');
    files=sort_nat({files.name});
    
    
    cellshape={};
    figure;
    for i=1:numel(files)
        im=imread(files{i});
        imshow(im);
        if ctb(i)==1
            title('Draw CTB cell contours, then press Enter. If cell is bad, skip drawing and directly press Enter.');
        else
            title('Draw non-CTB cell contours, then press Enter. If cell is bad, skip drawing and directly press Enter.');
        end
        hold on;
        cellshape{i}=[];
        while 1
            [x,y]=myginput(1,'cross');
            if isempty(x)
                break
            else
                if isempty(cellshape{i})
                    scatter(x,y,3,[1 0 0],'filled');
                else
                    plot([cellshape{i}(end,1),x],[cellshape{i}(end,2),y],'-r','Linewidth',2);
                end
                cellshape{i}=[cellshape{i};x,y];
            end
        end
        hold off;
        if ~isempty(cellshape{i})
            cellshapeconvex{i}=convhull(cellshape{i});
        else
            cellshapeconvex{i}=[];
        end
    end
    
    close;
    %% tranform coordinates
    cellshape1=cellshape;
    for i=1:numel(cellshape)
        m=mean(cellshape{i});
        m=floor(m/200)*200;
        cellshape{i}=cellshape{i}-m;
    end
    
        
    save temp1
    %% find dots and count.
    origfiles=dir('origcrop*.tif');
    origfiles=sort_nat({origfiles.name});
    if ~exist('hybthresh1','var')
        hybthresh1=[30 30 30 30];
    end
    
    rols={};
    filtrols={};
    counts=[];
    for m=1:numel(files)
        if ~isempty(cellshapeconvex{m})
            r={};rsub=r;rsubfilt=rsub;
            for n=1:4
                a=imread(origfiles{m},n);
                CC = bwconncomp(imregionalmax(imreconstruct(max(a-hybthresh1(n),0),a)));
                r{n}=zeros(length(CC.PixelIdxList),1);
                for i=1:length(CC.PixelIdxList)
                    [~,I]=max(a(CC.PixelIdxList{i}));
                    r{n}(i)=CC.PixelIdxList{i}(I); %linear indexed peak positions
                end
                [y,x]=ind2sub(size(a),r{n});
                rsub{n}=[x,y];%rsub is in x, y
                %check if dots are within the convex hull
                rsubfilt{n}=rsub{n}(inpolygon(rsub{n}(:,1),rsub{n}(:,2),cellshape{m}(cellshapeconvex{m},1),cellshape{m}(cellshapeconvex{m},2)),:);
                counts(m,n)=size(rsubfilt{n},1);
            end
            rols{m}=rsub;
            filtrols{m}=rsubfilt;
        else
            rols{m}={};
            filtrols{m}={};
            counts(m,1:4)=0;
        end
    end
    
    save('results.mat','rols','filtrols','counts','cellshape','cellshapeconvex','cellpos','contours','depths','angles','ctb');
    %%
    cd ..
end