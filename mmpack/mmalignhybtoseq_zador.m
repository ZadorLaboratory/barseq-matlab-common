function [templatesum,alignedim]=mmalignhybtoseq_zador(hybname,seqname,ch,chprofile,radius,nuclearch,cycle,chradius,tform_nuclear)
% align a single channel from hyb files to the sum of four sequencing
% channels from the first seq file. ch indicates the channel in hyb name to
% use for registration. chprofile is bleedthrough profile. radius is the
% radius for background subtraction. cycle is the sequencing cycle to align to. 
%This version also allows registration of the nuclear stain image after probes. 

%check if there is a nuclear channel.

if ~exist('nuclearch','var')
    nuclearch=0;
end

if ~exist('cycle','var')
    cycle=1;
end

if ~exist('chradius','var')
    chradius=30;
end

% check if hybfiles exist. If not, go into original folder and check again.
hybfiles=dir(['*',hybname,'*.tif']);
if ~isempty(hybfiles)
    hybfiles=sort_nat({hybfiles.name});
    o=0;
else
    cd original
    hybfiles=dir(['*',hybname,'*.tif']);
    hybfiles=sort_nat({hybfiles.name});
    cd ..
    o=1;
end

% check if seqfiles exist. If not, go into original folder and check again.
seqfiles=dir(['*',seqname,'*.tif']);
if ~isempty(seqfiles)
    seqfiles=sort_nat({seqfiles.name});
else
    cd original
    seqfiles=dir(['*',seqname,'*.tif']);
    seqfiles=sort_nat({seqfiles.name});
    seqfiles=cellfun(@(x) ['original/',x],seqfiles,'UniformOutput',0);
    cd ..
end

template1=imread(seqfiles{cycle},1);
template2=imread(seqfiles{cycle},2);
template3=imread(seqfiles{cycle},3);
template4=imread(seqfiles{cycle},4);
templatesum=double(template1+template2+template3+template4)./max(max(double(template1+template2+template3+template4)));

templatedic=imread(seqfiles{cycle},5);
%align hyb files ch channel to the first seqfiles.
if ~exist('tform_nuclear','var')
    tform_nuclear=repmat(affine2d,1,numel(hybfiles));
end

for i=1:length(hybfiles)
    if o==0
        info=imfinfo(hybfiles{i});
        im=zeros(info(1).Height,info(1).Width,length(info));
        for n=1:length(info)
            im(:,:,n)=medfilt2(imread(hybfiles{i},n));
        end
    else
        info=imfinfo(['original/',hybfiles{i}]);
        im=zeros(info(1).Height,info(1).Width,length(info));
        for n=1:length(info)
            im(:,:,n)=medfilt2(imread(['original/',hybfiles{i}],n));
        end
    end
    
    %subtract background
    ball=strel('ball', radius, radius);
    im(:,:,1:4)=im(:,:,1:4)-imopen(im(:,:,1:4),ball);
    
    %if there is a nuclear channel, also background subtract the nuclear
    %image
    if nuclearch>0
        im(:,:,nuclearch)=im(:,:,nuclearch)-imopen(im(:,:,nuclearch),ball);
    end
    
    %correct for bleeding
    im(:,:,1:4)=reshape(uint16(double(reshape(im(:,:,1:4),[],4))/chprofile),size(im,1),size(im,2),4);
    %subtract background of the rolony channel using smaller chradius
    im(:,:,ch)=im(:,:,ch)-imopen(im(:,:,ch),strel('ball',chradius,chradius));
    

%     %alignment using ECC
%     par.transform = 'euclidean';
%     par.levels = 4;
%     par.iterations = 100;
%     ransacWarp=iat_ecc(im(50:end-50,50:end-50,ch),templatesum(50:end-50,50:end-50),par);
%     %ransacWarp=iat_ecc(im(:,1300:end,ch),templatesum(:,1300:end),par);
%     
%     [M,N]=size(template4);
%     
%     %transform images
%     for n=1:length(info)
%         [alignedim(:,:,n),~]=iat_inverse_warping(im(:,:,n),ransacWarp,par.transform,1:N, 1:M);
%     end
    %align dic with imregcorr first
    imagedic=im(:,:,6);templatedic=imread(seqfiles{cycle},5);
    %tform0=imregcorr(imagedic(200:end-200,200:end-200),templatedic(200:end-200,200:end-200),'translation');
    tform0=tform_nuclear{i};
    Rfixed1=imref2d(size(imagedic));
    alignedim0=im;
    for n=1:length(info)

        alignedim0(:,:,n)=imwarp(im(:,:,n),tform0,'OutputView',Rfixed1);

    end

    
    %align using imregtform or imregcorr
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.InitialRadius = optimizer.InitialRadius/50;
    optimizer.GrowthFactor=1.01;
    optimizer.MaximumIterations=optimizer.MaximumIterations*10;
    tform = imregtform(alignedim0(100:end-100,100:end-100,ch), templatesum(100:end-100,100:end-100), 'translation', optimizer, metric,'PyramidLevels',4);
    % tform=imregcorr(imagesum(100:end-100,100:end-100), templatesum(100:end-100,100:end-100),'rigid','Window',0);
    Rfixed=imref2d(size(templatesum));
    alignedim=im;
    for n=1:length(info)
        alignedim(:,:,n)=imwarp(alignedim0(:,:,n),tform,'OutputView',Rfixed);
    end

    
    alignedfile=strcat('aligned',hybfiles{i});
    imwrite(uint16(alignedim(:,:,1)),alignedfile);
    for n=2:size(im,3)
        imwrite(uint16(alignedim(:,:,n)),alignedfile, 'WriteMode','append');
    end
    if o~=1
        movefile(hybfiles{i},['original/',hybfiles{i}]);
    end
end
