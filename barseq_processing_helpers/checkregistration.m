
function checkregistration(fname,tformfname,scaling)
% make stitched images of RGB files to allow visual inspection of
% stitching/registration
[folders,pos]=get_folders();
if ~exist('scaling','var')
    scaling=1;
end

if ~exist('tformfname','var')
    tformfname='40xto10x.mat';
end

load(tformfname,'tform40xto10x');
reg_folder_name=['checkregistration'];
mkdir(reg_folder_name)
uniqpos=sort_nat(unique(pos));
[~,posidx]=ismember(pos,uniqpos);
%
f=dir([folders{1},'/RGB/*',fname,'*.tif']);
f={f.name};
parfor i=1:numel(uniqpos)
    %get image size
    %%
    img10xfilename=['10x/',uniqpos{i},'.tif'];
    img10xinfo=imfinfo(img10xfilename);
    subfolders=find(posidx==i);
    im=repmat({zeros(img10xinfo(1).Height,img10xinfo(1).Width,3,'uint8')},numel(f),1);
    %%
    for n=1:numel(subfolders)
        imgfiles=dir([folders{subfolders(n)},'/RGB/*',fname,'*.tif']);
        imgfiles={imgfiles.name};
        for m=1:numel(imgfiles)
            im1=imread([folders{subfolders(n)},'/RGB/',imgfiles{m}]);
            im1=imwarp(im1,tform40xto10x{subfolders(n)},'OutputView',imref2d([size(im{m},1),size(im{m},2)]));
            im1(1:10,1:10,:)=150;%add a border so that it's easy to count
            im{m}=max(im{m},im1);
        end
    end
    
    %%
    for m=1:numel(im)
        im{m}=im{m}*scaling;
        imwrite(imresize(uint8(im{m}),0.2),[reg_folder_name,'/',uniqpos{i},f{m}]);
    end
end
    fprintf('Stitched images in folder: checkregistration%s. Visually check to make sure they are registered correctly.\n',fname)

end

