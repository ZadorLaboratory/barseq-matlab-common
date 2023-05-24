function currrol=mmsegmentbc11(varargin)

%% parse inputs

p=inputParser;
addParameter(p,'bgn',50,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'lthresh',800,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'projthresh',0.2,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'rejthresh',0.8,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'cellmasksize',30,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'cellratio',0.5,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'maxcellrange',300,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'mincellarea',400,@(x)isnumeric(x)&&isscalar(x));

parse(p,varargin{:});
bgn  = p.Results.bgn;
lthresh = p.Results.lthresh;
projthresh  = p.Results.projthresh;
rejthresh  = p.Results.rejthresh;
cellmasksize  = p.Results.cellmasksize;
cellratio  = p.Results.cellratio;
maxcellrange  = p.Results.maxcellrange;
mincellarea  = p.Results.mincellarea;

% bgn=50;
% lthresh=800;
% projthresh=0.2;
% rejthresh=0.8;
% cellmasksize=30;
% cellratio=0.5;
% maxcellrange=300;
% mincellarea=400;



%segment cells using barcodes
files=dir('*aligned*BC*.tif');
files=sort_nat({files.name});

info=imfinfo(files{1});
im=zeros(info(1).Height,info(1).Width,4,length(files));
for i=1:length(files)
    for n=1:4
        im(:,:,n,i)=imread(files{i},n);
    end
end
%img=gpuarray(im);
img=im;
%find cells in each cycle/channel
imblurg=imgaussfilt(img,3);
%imblur=reshape(gather(imblurg),size(imblurg,1),size(imblurg,2),[]);
imblur=reshape(imblurg,size(imblurg,1),size(imblurg,2),[]);

r={};rsub={};
for n=1:4
    
    a=imblur(:,:,n);
    CC = bwconncomp(imregionalmax(imreconstruct(max(a-lthresh,0),a)));
    r{n}=zeros(length(CC.PixelIdxList),1);
    for i=1:length(CC.PixelIdxList)
        [~,I]=max(a(CC.PixelIdxList{i}));
        r{n}(i)=CC.PixelIdxList{i}(I); %linear indexed peak positions
    end
    [y,x]=ind2sub(size(imblur(:,:,n)),r{n});
    rsub{n}=[x,y];%rsub is in x, y, consistent with lroi1
    
    
end

currrol=[rsub{1};rsub{2};rsub{3};rsub{4}];







