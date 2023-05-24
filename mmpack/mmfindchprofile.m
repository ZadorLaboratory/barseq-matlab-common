function chprofile=mmfindchprofile(filename,baseline)
%find the channel bleedthrough profile and balance channels from the first
%sequencing image contianing filenmae
% 
% intprofile=[1 1 0.7 0.7];
% 
if ~exist('baseline')
    baseline=10;
end



files=dir(['*',filename,'*.tif']);
files=sort_nat({files.name});
im=imread(files{1});
for n=2:4
    im(:,:,n)=imread(files{1},n);
end


%bgn subtraction
radius=100;
ball=strel('ball', radius, radius);
im=im-imopen(im,ball);

%sparse-nmf
option.beta=0.3;
[A,Y]=sparsenmfnnls(reshape(double(im(200:end-200,200:end-200,:)),[],4)',4,option); %Y should be correspond to the four channels. Can tune the sparseness constraint if necessary (doesn't seem to affect much)
Y(Y<baseline)=nan;
Y=nanmedian(Y,2);


[~,I]=max(A,[],1);
[~,I]=sort(I);
chprofile=A(:,I)';%This assumes that for each color, the strongest channel is still the intended channel. (can't have two channels with max signal in the same imaging channel)

chprofile=chprofile.*repmat(Y(I),1,4);%normalize intensity
chprofile=chprofile./repmat(max(chprofile(:)),4,4);

