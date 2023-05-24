function offset=niefindfluorfocus(posfname,focusfolder,initz,endz,stepz,ch)
% gien a list of positions and a z-stack of fluorescent images, find the
% focus offset for each position.

%% new version to allow defining which channel to use for sequencing
if ~exist('ch','var')
    ch=1;
end




%%
data=readmatrix(posfname);
data2=data;
offset=zeros(size(data,1),1);
%%
cd(focusfolder);
files=dir('*.tif');
files=sort_nat({files.name});

fileinfo=imfinfo(files{1});
pzc=cell(length(files),3);
for m=1:length(files)
    pzc(m,:)=textscan(files{m},'%*u %u %u %u %*s','Delimiter',{'xy','z','c','.'});
end
pzc=cell2mat(pzc);
parfor n=1:size(data,1)
    %%
    fidx=find(pzc(:,1)==n&pzc(:,3)==ch);
    im=zeros(fileinfo(1).Height,fileinfo(1).Width,numel(fidx));
    for i=1:numel(fidx)
        im(:,:,i)=imread(files{fidx(i)});
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

end
cd ..
%%
data2(:,3)=data2(:,3)+offset;%assuming same direction betwen the piezo and the focus
writematrix(data2,['offset',posfname],'Delimiter',';')
end