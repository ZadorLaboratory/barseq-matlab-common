function mmalignhybtoseqsinglech(hybname,seqname,hybch,seqch,radius,nuclearch,cycle)
% align a single channel from hyb files to the sum of four sequencing
% channels from the first seq file. ch indicates the channel in hyb name to
% use for registration. chprofile is bleedthrough profile. radius is the
% radius for background subtraction. cycle is the sequencing cycle to align to. 
%This version also allows registration of the nuclear stain image after probes. 
%This version performs registration on a single channel

%check if there is a nuclear channel.
if ~exist('nuclearch','var')
    nuclearch=0;
end

if ~exist('cycle','var')
    cycle=1;
end

% check if seqfiles exist. If not, go into original folder and check again.
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

% check if seqfiles exist. If not, go into aligned folder and check again.
seqfiles=dir(['*',seqname,'*.tif']);
if ~isempty(seqfiles)
    seqfiles=sort_nat({seqfiles.name});
else
    cd aligned
    seqfiles=dir(['*',seqname,'*.tif']);
    seqfiles=sort_nat({seqfiles.name});
    seqfiles=cellfun(@(x) ['aligned/',x],seqfiles,'UniformOutput',0);
    cd ..
end

templatesum=imread(seqfiles{cycle},seqch);

%align hyb files ch channel to the first seqfiles.
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
    
    %if there is a nuclear channel, also background subtract the nuclear
    %image
    ball=strel('ball', radius, radius);
    im(:,:,hybch)=im(:,:,hybch)-imopen(im(:,:,hybch),ball);
    if nuclearch>0
        
        im(:,:,nuclearch)=im(:,:,nuclearch)-imopen(im(:,:,nuclearch),ball);
        
    end
    
    %alignment
    par.transform = 'euclidean';
    par.levels = 4;
    par.iterations = 100;
    ransacWarp=iat_ecc(im(50:end-50,50:end-50,hybch),templatesum(50:end-50,50:end-50),par);
    
    [M,N]=size(templatesum);
    
    %transform images
    for n=1:length(info)
        [alignedim(:,:,n),~]=iat_inverse_warping(im(:,:,n),ransacWarp,par.transform,1:N, 1:M);
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
