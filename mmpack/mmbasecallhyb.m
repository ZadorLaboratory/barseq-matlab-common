function [id1,lroi1,sig1]=mmbasecallhyb(codebookhyb,hybthresh1,bgn,psf,chprofile)
%basecalling gene rolonies based on codebook using a soft basecalling based
%on rolonies identified in the first cycle based on intensity. works on
%single position images in the aligned folder. If psf is not provided,
%guess psf using deconvblind.chprofile provides channel bleedthrough
%profile if images have not been fixed completely.


%% parse inputs
codes=num2str(cell2mat(codebookhyb(:,2)));

if ~exist('hybthresh1')
    hybthresh1=[10 10 20 20];
end
if ~exist('chprofile')
    chprofile=eye(4);
end
%noise=uint16(rolthresh1/2);
noise=[1 1 1 1];

%% read 1st cycle and find rolonies
seqfiles=dir(fullfile('*hyb*.tif'));
[seqfiles,~]=sort_nat({seqfiles.name});
%%
if length(seqfiles)==1
    %%
    m=1;
    sig1=[];lroi1=[];id1=[];
    imoriginal=[];
    
    
    for n=1:4
        imoriginal(:,:,n)=imread(seqfiles{m},n);
    end
    imfixed=reshape(uint16(double(reshape(imoriginal,[],4))/chprofile),size(imoriginal,1),size(imoriginal,2),4);%subtract camera baseline and correct for bleeding
    %%
    parfor n=1:4
        lim(:,:,n)=deconvlucy(imgaussfilt(imfixed(:,:,n),1),psf{n},10); %current deconv settings are on the aggresive side.
        im(:,:,n)=imgaussfilt(imfixed(:,:,n),1);
    end
    
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
            sig1=[sig1;a(r{n}(a(r{n})>bgn(n)))];
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
else
    sig1={};lroi1={};id1={};
    for m=1:length(seqfiles)
        lim=[];im=[];
        imoriginal=[];
        
        
        for n=1:4
            imoriginal(:,:,n)=imread(seqfiles{m},n);
        end
        imfixed=reshape(uint16(double(reshape(imoriginal,[],4))/chprofile),size(imoriginal,1),size(imoriginal,2),4);%subtract camera baseline and correct for bleeding
        
        
        parfor n=1:4
            lim(:,:,n)=deconvlucy(imgaussfilt(imfixed(:,:,n),1),psf{n},10); %current deconv settings are on the aggresive side.
        end
        
        parfor n=1:4
            im(:,:,n)=imgaussfilt(imfixed(:,:,n),1);
        end
        
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
        sig1{m}=[];lroi1{m}=[];id1{m}=[];
        for n=1:4
            c=find(codes(:,m)==num2str(n));
            if ~isempty(c)
                a=im(:,:,n);
                sig1{m}=[sig1{m};a(a(r{n})>bgn(n))];
                lroi1{m}=[lroi1{m};rsub{n}(a(r{n})>bgn(n),:)];
                id1{m}=[id1{m};repmat(c,sum(a(r{n})>bgn(n)),1)];
            end
        end
    end
end



idhyb1=id1;
lroihyb1=lroi1;
sighyb1=sig1;
save('hybcalls.mat','idhyb1','lroihyb1','sighyb1');



%
% MIJ.exit;
end