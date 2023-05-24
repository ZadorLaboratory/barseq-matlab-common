function currrol=mmsegmentbc1(varargin)

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


%%
addpath('C:\Fiji.app\scripts');
javaaddpath 'C:\Program Files\MATLAB\R2018b\java\mij.jar'

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

currrol=[];
%find rolonies after gaussian blur with sigma=3
Miji(false);
for n=1:size(img,3)
    MIJ.createImage(imblur(:,:,n));
    MIJ.run('Find Maxima...', ['noise=',num2str(lthresh),' output=List']);
    currrol=[currrol;MIJ.getResultsTable];
    %cleanup
    MIJ.run('Clear Results');
    MIJ.closeAllWindows
end



MIJ.exit;











