function check_stitching(tformfname)
    if ~exist('tformfname','var')
        tformfname='40xto10x.mat';
    end
    
    load(tformfname,'tform40xto10x');
    [folders,pos]=get_folders();
    
    mkdir stitched40xto10x
    uniqpos=sort_nat(unique(pos));
    [~,posidx]=ismember(pos,uniqpos);
    %
    parfor i=1:numel(uniqpos)
        %get image size
        %%
        img10xfilename=['10x/',uniqpos{i},'.tif'];
        im10x=double(imread(img10xfilename,1))+double(imread(img10xfilename,2))+double(imread(img10xfilename,3))+double(imread(img10xfilename,4));
        subfolders=find(posidx==i);
        im=zeros(size(im10x));
        %%
        for n=1:numel(subfolders)
            imgfilename=[folders{subfolders(n)},'/original/n2vgeneseq01.tif'];
            im1=double(imread(imgfilename,1))+double(imread(imgfilename,2))+double(imread(imgfilename,3))+double(imread(imgfilename,4));
            im1=imwarp(im1,tform40xto10x{subfolders(n)},'OutputView',imref2d(size(im)));
            im=max(im,im1);
        end
        %%
        imwrite(imresize(uint16(im10x),0.2),['stitched40xto10x/small',uniqpos{i},'.tif']);
        imwrite(imresize(uint16(im),0.2),['stitched40xto10x/small',uniqpos{i},'.tif'],'WriteMode','Append');
        
        
    end
    fprintf('Registered 40x images in /stitched40xto10x. Visually check to make sure they are registered correctly')

end