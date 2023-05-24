%script for basecalling from \seq0 from prealigned images
addpath(genpath(pwd));
old=cd('seq2');
%do a rough alignment

%find rolonies in all cycles and get intensities
files=dir(fullfile('crop*'));
files={files.name};
rol=cell(4,1);
sig=rol;
parfor i=1:length(files)
    im1=[];
    for n=1:4
        im1(:,:,n)=imread(files{i},n);
    end
    rol{i}=findrolonies(im1,20);
    sig{i}=readintensities(rol{i},im1,[1 1 1;1 1 1;1 1 1]);%read intensities using the original rol0
    
end
rol0=rol;

save allrolonies;

%% rigid alignment
%do a precise alignment using rigid ICP
tt{1}=[0;0;0];
tr{1}=eye(3);

for i=2:length(files)
    [tr{i},tt{i}]=icp([rol{1}';ones(1,length(rol{1}))], ...
        [rol{i}';ones(1,length(rol{i}))],200);
    %align rolonies, currently only doing translation
    rol{i}=rol{i}+repmat(tt{i}(1:2)',length(rol{i}),1);
end
% find matching rolonies
%use the 1st cycle as reference, find the rolonies in later cycles that are
%closest to the 1st cycle rolonies as matching rolonies

I=zeros(length(rol{1}),length(files));
mindist=I;
I(:,1)=(1:length(rol{1}))';
for i=2:length(files)
    d=pdist2(rol{1},rol{i});
    [mindist(:,i),I(:,i)]=min(d,[],2);%mindist is the distance between 1st 
    %cycle rolony and nth cycle rolony. I is the index for the
    %corresponding rolony in nth channel.
end

%basecall
seq=zeros(length(rol{1}),length(files));
scores=seq;
for i=1:length(files)
    [~,seq(:,i)]=max(sig{i}(I(:,i),:),[],2);
    scores(:,i)=max(sig{i}(I(:,i),:),[],2)./sqrt(sum(sig{i}(I(:,i),:).^2,2));
end

seqC(seq==1)='G';
seqC(seq==2)='T';
seqC(seq==3)='A';
seqC(seq==4)='C';
seqC=reshape(seqC,[],length(files));
save temp2

%% non-rigid ICP

%try non-rigid ICP
rolnr=rol0;
target.vertices=[rol0{1},ones(length(rol0{1}),1)];
target.faces=delaunay(rol0{1}(:,1),rol0{1}(:,2));
Options=struct;
Options.normalWeighting = 0;
for i=2:length(files)
    source.vertices=[rol0{i},ones(length(rol0{i}),1)];
    source.faces=delaunay(rol0{i}(:,1),rol0{i}(:,2));
    [pointsTransformed, X] = nricp(source, target,Options);
    rolnr{i}=pointsTransformed(:,1:2);
end

% find matching rolonies
%use the 1st cycle as reference, find the rolonies in later cycles that are
%closest to the 1st cycle rolonies as matching rolonies

Inr=zeros(length(rolnr{1}),length(files));
mindistnr=Inr;
Inr(:,1)=(1:length(rolnr{1}))';
for i=2:length(files)
    d=pdist2(rolnr{1},rolnr{i});
    [mindistnr(:,i),Inr(:,i)]=min(d,[],2);%mindist is the distance between 1st 
    %cycle rolony and nth cycle rolony. I is the index for the
    %corresponding rolony in nth channel.
end

%basecall
seqnr=zeros(length(rolnr{1}),length(files));
scoresnr=seqnr;
for i=1:length(files)
    [~,seqnr(:,i)]=max(sig{i}(Inr(:,i),:),[],2);
    scoresnr(:,i)=max(sig{i}(Inr(:,i),:),[],2)./sqrt(sum(sig{i}(Inr(:,i),:).^2,2));
end

seqCnr(seqnr==1)='G';
seqCnr(seqnr==2)='T';
seqCnr(seqnr==3)='A';
seqCnr(seqnr==4)='C';
seqCnr=reshape(seqCnr,[],length(files));

save temp3


