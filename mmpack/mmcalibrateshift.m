function chshift=mmcalibrateshift(fname)
% given a series of stacked tif images, calculate the shift acorss
% channels.

%%
chshift=[];

files=dir(['*',fname,'*']);
files={files.name};
%%
for i=1:numel(files)
    %% read files
    fileinfo=imfinfo(files{i});
    im=zeros(fileinfo(1).Height,fileinfo(1).Width,numel(fileinfo));
    for n=1:size(im,3)
        im(:,:,n)=imread(files{i},n);
    end
    %% calculate shift
    tform={};
    tform1=[];
    for n=1:size(im,3)
        tform{n}=imregcorr(im(:,:,n),im(:,:,1),'translation');
        tform1(:,:,n)=tform{n}.T;
    end
    chshift(:,:,i)=reshape(tform1(3,1:2,:),2,[])';
end


