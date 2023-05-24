function zoffset=mmfinddicfocus(initz,endz,stepsize)
%Find focus point in stacks of dic images
% finding std dip is not reliable.
folders=dir('Pos*');
folders=sort_nat({folders.name});
for nn=1:numel(folders)
    %%
    cd(folders{nn})
    %% read files
    files=dir(['*.tif']);
    files=sort_nat({files.name});
    fileinfo=imfinfo(files{1});
    im=zeros(fileinfo(1).Height,fileinfo(1).Width,numel(files));
    parfor n=1:numel(files)
        im(:,:,n)=imread(files{n});
    end

    %% find focus
    im1=reshape(im,[],size(im,3));
    S=squeeze(std(im,0,[1 2]));
    figure;plot(S/mean(S))
    [~,I]=min(S);
        
    if endz>initz
        zoffset(nn)=initz+stepsize*(I-1);
    else
        zoffset(nn)=initz-stepsize*(I-1);
    end
    cd ..
end

save('zoffset.mat','zoffset');


end

