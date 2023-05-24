function makelocalimg2d(foldername,range)
%locally align maxproj images  to first cycle (and rolony coordinates in
%first cycle). basecall using pixel values



%% basecall first cycle: balance channels, read rolony signals, 
%read first cycle image.
olddir=cd(foldername);
files=dir('*.tif');
im=cell(length(files),1);
for i=1:length(files)
    im{i}=imread(files(i).name,1);
    for n=2:4
        im{i}(:,:,n)=imread(files(i).name,n);
    end
end
cd ..

%cut cycle 1 image into 400px by 400px tiles to prepare for alignment and
%basecalling of subsequent cycles. tiles overlap each other by 50%,so that
%every pixel is covered twice except at the edge. In each tile, only call
%rolonies at the center 300px by 300px range.
mkdir tiles;


tile1=1:300:floor(size(im{1},1)/300-1)*300+1;
tile2=1:300:floor(size(im{1},2)/300-1)*300+1;
tile1=repmat(tile1,length(tile2),1);
tile2=repmat(tile2',1,length(tile1));
tile=reshape(cat(3,tile1,tile2),[],2);%this has the starting index of each tile.

if range(1)==0
    range(1)=1;
    range(2)=size(tile,1);
end

for i=range(1):range(2)
    
    if tile(i,1)+798>size(im{1},1)
        e1=size(im{1},1);
    else
        e1=tile(i,1)+399;
    end
    
    if tile(i,2)+798>size(im{1},2)
        e2=size(im{1},2);
    else
        e2=tile(i,2)+399;
    end
    
    zmax1tile=im{1}(tile(i,1):e1,tile(i,2):e2,:);
    imwrite(zmax1tile(:,:,1),['tiles\',num2str(i),'.tif'])
    for m=2:4
        imwrite(zmax1tile(:,:,m),['tiles\',num2str(i),'.tif'],'WriteMode','Append');
    end
    imwrite(uint16(sum(zmax1tile,3)),['tiles\zmax',num2str(i),'.tif'])
    
end


%% for subsequent cycles:
for n=2:length(files)
    %cut images into 500px x 500 px images, and align each small stack with cycle 1.
    for i=range(1):range(2)
        
        if tile(i,1)+798>size(im{1},1)
            e1=size(im{1},1);
        else
            e1=tile(i,1)+399;
        end
        
        if tile(i,2)+798>size(im{1},2)
            e2=size(im{1},2);
        else
            e2=tile(i,2)+399;
        end
        
        im1tile=im{1}(tile(i,1):e1,tile(i,2):e2,:);
        im2tile=im{n}(tile(i,1):e1,tile(i,2):e2,:);
        
        
        %align tiles in xy plane allowing rotation and translation.
        %Then align a 3x3 area around each rolony along z using xcorr
        
        %make z-proj, calculate tform, and apply to the whole stack.
        %optimizer = registration.optimizer.OnePlusOneEvolutionary;
        %metric = registration.metric.MattesMutualInformation;
        %tform1 = imregtform(sum(im2tile,3),sum(im1tile,3),'rigid',optimizer,metric);%doesn't work for boundary slices due to the bounary being very bright. Need to detect boundary and exclude from alignment
        tform1 = imregcorr(sum(im2tile,3),sum(im1tile,3),'rigid','Window',0);%doesn't work for boundary slices due to the bounary being very bright. Need to detect boundary and exclude from alignment
        Rfixed=imref2d(size(im1tile(:,:,1)));
        warpedim2tile=imwarp(im2tile,tform1,'OutputView',Rfixed);%channel-summed stack warped in xy, used for z-alignment

    for m=1:4
        imwrite(warpedim2tile(:,:,m),['tiles\',num2str(i),'.tif'],'WriteMode','Append');
    end
        imwrite(uint16(sum(warpedim2tile,3)),['tiles\zmax',num2str(i),'.tif'],'WriteMode','Append');
    
       
    end
end


