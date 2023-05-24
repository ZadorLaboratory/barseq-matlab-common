function mmcolorizebydepth(minz,maxz,r)
% Colorize rolonies by depths and generate max projection.
%%
if ~exist('r','var')
    r=10;
end


files=dir('*.tif');
files={files.name};

%parse file names
tcz=cell(length(files),3);
for m=1:length(files)
    tcz(m,:)=textscan(files{m},'%*s %u %s %u','Delimiter',{'_','.'});
    %find sequencing channels
    
end
tcz(:,2)=cellfun(@cell2mat,tcz(:,2),'UniformOutput',false);

ch=unique(tcz(:,2));
z=cell2mat(tcz(:,3));
%% for each channel, build a stack
mkdir colorizedepth
for i=1:numel(ch)
    im=[];
    idx=find(ismember(tcz(:,2),ch(i))&z>=minz&z<=maxz);
    for n=1:numel(idx)
        im(:,:,n)=double(imread(files{idx(n)}));
    end
    cmap1=jet(size(im,3));
    imc=reshape(reshape(im,[],size(im,3))*cmap1,size(im,1),size(im,2),[]);
    for n=1:size(imc,3)
        imc(:,:,n)=imtophat(imc(:,:,n),strel('ball',r,r));
        imc(:,:,n)=(imc(:,:,n))./max(max(imc(:,:,n)))*255;
    end
    imwrite(uint8(imc),['colorizedepth/',ch{i},'.tif']);
    
end