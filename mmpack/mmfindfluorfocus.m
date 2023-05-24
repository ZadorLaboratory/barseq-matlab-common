function offset=mmfindfluorfocus(fname,ffolder,initz,endz,stepz)
% gien a list of positions and a z-stack of fluorescent images, find the
% focus offset for each position.

%%
data=loadjson(fname);
data2=data;
offset=zeros(numel(data.POSITIONS),1);
%%
parfor n=1:numel(data.POSITIONS)
    %%
    pos=data.POSITIONS{n}.LABEL;
    cd([ffolder,'/',pos]);
    files=dir('*.tif');
    files={files.name};
    fileinfo=imfinfo(files{1});
    im=zeros(fileinfo(1).Height,fileinfo(1).Width,numel(files));
    for i=1:numel(files)
        im(:,:,i)=imread(files{i});
    end
    m=zeros(size(im,3),1);
    %%
    %calculate focus, using mean intensity excluding the top 0.05%
    for i=1:size(im,3)
        im1=reshape(im(:,:,i),1,[]);
        im1(im1>prctile(im1,99.95))=prctile(im1,99.95);
        m(i)=mean(im1);
    end
    [~,I]=max(medfilt1(m,3));
    if endz>=initz %the piezo and manual focus are reversed
        offset(n)=initz+(I-1)*stepz;
    else
        offset(n)=initz-(I-1)*stepz;
    end
    %%
    %find focusing device

    cd ../..
end

%%
for n=1:numel(data.POSITIONS)
    for i=1:numel(data.POSITIONS{n}.DEVICES)
        data2.POSITIONS{n}.PROPERTIES.OverlapUmX=0;
        data2.POSITIONS{n}.PROPERTIES.OverlapUmY=0;
        data2.POSITIONS{n}.PROPERTIES.OverlapPixelsX=0;
        data2.POSITIONS{n}.PROPERTIES.OverlapPixelsY=0;
        if ismember(data.POSITIONS{n}.DEVICES{i}.DEVICE,{'ManualFocus'})
            data2.POSITIONS{n}.DEVICES{i}.X=data2.POSITIONS{n}.DEVICES{i}.X-offset(n);%the piezo and manual focus are reversed
        end
    end
end

json=savejson('',data2,['offset',fname]);
end