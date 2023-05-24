function tform=niecalibrateduocam(foldername,ch,transformationtype)
% calculate image transformation from dual-cam scopes
% ch gives a pair of channels to calibrate
%%
cd(foldername);
files=dir('*.tif');
files=sort_nat({files.name});

fileinfo=imfinfo(files{1});
pzc=cell(length(files),3);
for m=1:length(files)
    pzc(m,:)=textscan(files{m},'%*u %u %u %u %*s','Delimiter',{'xy','z','c','.'});
end
pzc=cell2mat(pzc);


fidx=find(pzc(:,3)==ch(1));
im1=zeros(fileinfo(1).Height,fileinfo(1).Width,numel(fidx));
for i=1:numel(fidx)
    im1(:,:,i)=imread(files{fidx(i)});
end
im1=max(im1,[],3);

fidx=find(pzc(:,3)==ch(2));
im2=zeros(fileinfo(1).Height,fileinfo(1).Width,numel(fidx));
for i=1:numel(fidx)
    im2(:,:,i)=imread(files{fidx(i)});
end
im2=max(im2,[],3);

%% manual registration by point pairs
[mp,fp]=cpselect(im2./max(im2,[],'all'),im1./max(im1,[],'all'),'Wait',true);
imref=imref2d(size(im2));
if ~isempty(mp)
    tform0=fitgeotrans(mp,fp,'affine');
    im2_0=imwarp(im2,tform0,'outputView',imref);
    figure;imshowpair(min(im2_0./max(im2_0,[],'all')*2,1),min(im1./max(im1,[],'all')*2,1))
else
    im2_0=im2;
    tform0=affine2d;
end

%%
%refine with correlation based registration
[optimizer,metric]=imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/5; %trial 1 conditions
optimizer.Epsilon=optimizer.Epsilon/10;
optimizer.MaximumIterations=optimizer.MaximumIterations*5;
tform1=imregtform(im2_0,im1,transformationtype,optimizer,metric);


cd ..
tform=tform1;
tform.T=tform0.T*tform1.T;

im2t=imwarp(im2,tform,'outputView',imref);
figure;imshowpair(min(im2t./max(im2t,[],'all')*2,1),min(im1./max(im1,[],'all')*2,1))

end

