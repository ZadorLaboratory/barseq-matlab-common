%read file and median filter
info=imfinfo(filename);
for i=1:4
    im(:,:,i)=medfilt2(imread(filename,i));
end

end


%bgn subtraction
radius=100;
ball=strel('ball', radius, radius);
im=im-imopen(im,ball);

%% read 1st cycle and find rolonies
seqfiles=dir(fullfile('*stitched*.tif'));
[seqfiles,~]=sort_nat({seqfiles.name});
lroi1=[];
id1=[];
sig1=[];
%%
for m=1:length(seqfiles)
    %%
    parfor n=1:4
        lim(:,:,n)=deconvlucy(imgaussfilt(imread(seqfiles{m},n),1),psf{n},10); %current deconv settings are on the aggresive side.
    end
    
    parfor n=1:4
        im(:,:,n)=imgaussfilt(imread(seqfiles{m},n),1);
    end
    
    %%
    r={};rsub=r;
    
    parfor n=1:4
        a=lim(:,:,n);
        CC = bwconncomp(imregionalmax(imreconstruct(max(a-hybthresh1(n),0),a)));
        r{n}=zeros(length(CC.PixelIdxList),1);
        for i=1:length(CC.PixelIdxList)
            [~,I]=max(a(CC.PixelIdxList{i}));
            r{n}(i)=CC.PixelIdxList{i}(I); %linear indexed peak positions
        end
        [y,x]=ind2sub(size(lim(:,:,n)),r{n});
        rsub{n}=[x,y];%rsub is in x, y, consistent with lroi1
        
    end
    for n=1:4
        c=find(codes(:,m)==num2str(n));
        if ~isempty(c)
            a=im(:,:,n);
            sig1=[sig1;a(a(r{n})>bgn(n))];
            lroi1=[lroi1;rsub{n}(a(r{n})>bgn(n),:)];
            id1=[id1;repmat(c,sum(a(r{n})>bgn(n)),1)];
        end
    end
    
    %     %raw image
    %     r1={};r1sub=r1;
    %     parfor n=1:4
    %         a=lim(:,:,n);
    %         CC = bwconncomp(imregionalmax(imreconstruct(max(a-hybthresh1(n),0),a)));
    %         r1{n}=zeros(length(CC.PixelIdxList),1);
    %         for i=1:length(CC.PixelIdxList)
    %             [~,I]=max(a(CC.PixelIdxList{i}));
    %             r1{n}(i)=CC.PixelIdxList{i}(I); %linear indexed peak positions
    %         end
    %         [y,x]=ind2sub(size(lim(:,:,n)),r1{n});
    %         r1sub{n}=[y,x];
    %     end
%         %%
%         close all;
%        i=2;
%        a=im(:,:,i);
%     %    figure;imagesc(a);hold on;scatter(r1sub{i}(a(r1{i})>bgn(i),2),r1sub{i}(a(r1{i})>bgn(i),1),'+r');set(gca,'xlim',[1300 1500],'ylim',[1100 1300]);
%     figure;imagesc(a,[0 1000]);hold on;scatter(lroi1(id1==i,1),lroi1(id1==i,2),'+r');%set(gca,'xlim',[1300 1500],'ylim',[1100 1300]);
%     figure;imagesc(a,[0 300]);
end