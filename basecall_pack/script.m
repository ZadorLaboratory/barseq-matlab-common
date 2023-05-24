%script for basecalling from \seq
addpath(genpath(pwd));
old=cd('seq');
%do a rough alignment
Ultraviewseqalignment1('1seq');
old1=cd('aligned');
rgbout('aligned');
movefile('rgb*','..\rgb\');
%make a center crop
files=dir(fullfile('aligned*'));
files={files.name};
for i=1:length(files)
    im1=[];
    for n=1:4
        im1(:,:,n)=imread(files{i},n);
    end
    cropim1=im1(1100:1600,400:900,:);
    imwrite(cropim1(:,:,1),['crop',files{i}]);
    imwrite(cropim1(:,:,2),['crop',files{i}],'WriteMode','append');
    imwrite(cropim1(:,:,3),['crop',files{i}],'WriteMode','append');
    imwrite(cropim1(:,:,4),['crop',files{i}],'WriteMode','append');    
end

%find rolonies in all cycles
files=dir(fullfile('crop*'));
files={files.name};
rol=cell(4,1);

for i=1:length(files)
    im1=[];
    for n=1:4
        im1(:,:,n)=imread(files{i},n);
    end
    rol{i}=findrolonies(im1,20);
end
rol0=rol;
%do a precise alignment using rolonies
tt{1}=[0;0;0];
tr{1}=eye(3);

for i=2:length(files)
    [tr{i},tt{i}]=icp([rol{1}';ones(1,length(rol{1}))], ...
        [rol{i}';ones(1,length(rol{i}))],200);
    %align rolonies, currently only doing translation
    rol{i}=rol{i}+repmat(tt{i}(1:2)',length(rol{i}),1);
end

%% find matching rolonies
%use the 1st cycle as reference, find the rolonies in later cycles that are
%closest to the 1st cycle rolonies as matching rolonies

%check the distribution of min rolony distance between two cycles, 
d=pdist2(rol{1},rol{3});
rd=pdist2([rol{1}(:,1),rol{1}(randperm(length(rol{1})),2)],rol{3});

figure;histogram(min(d,[],2),0:0.4:10);
figure;histogram(min(rd,[],2),0:0.4:10);
% judging from the histograms, everything with less than 2-px distance
% should be real rolonies.

I=zeros(length(rol{1}),length(files));
mindist=I;
I(:,1)=(1:length(rol{1}))';
for i=2:length(files)
    d=pdist2(rol{1},rol{i});
    [mindist(:,i),I(:,i)]=min(d,[],2);%mindist is the distance between 1st 
    %cycle rolony and nth cycle rolony. I is the index for the
    %corresponding rolony in nth channel.
end
save temp2
%% basecall based on intensity in a 3x3 area.(could potentially switch to
%combine with the ID of rolonies.
sig=cell(4,1);
for i=1:length(files)
    im1=[];
    for n=1:4
        im1(:,:,n)=imread(files{i},n);
    end
    sig{i}=readintensities(rol0{i},im1,[1 1 1;1 1 1;1 1 1]);%read intensities using the original rol0
    
end
save temp3

%basecall
seq=zeros(length(rol{1}),length(files));
scores=seq;
for i=1:length(files)
    [~,seq(:,i)]=max(sig{i}(I(:,i),:),[],2);
end

%check the bleedthrough profile:
figure; scatter(sig{1}(:,1),sig{1}(:,2),'.')
figure; scatter(sig{2}(:,3),sig{1}(:,4),'.')


    
    



