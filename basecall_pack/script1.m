%script for basecalling from \seq1
addpath(genpath(pwd));
old=cd('seq11');
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
    cropim1=uint16(im1(1100:1600,400:900,:));
    imwrite(cropim1(:,:,1),['crop',files{i}]);
    imwrite(cropim1(:,:,2),['crop',files{i}],'WriteMode','append');
    imwrite(cropim1(:,:,3),['crop',files{i}],'WriteMode','append');
    imwrite(cropim1(:,:,4),['crop',files{i}],'WriteMode','append');    
end

%find rolonies in all cycles and get intensities
files=dir(fullfile('aligned*'));
files={files.name};
rol=cell(4,1);
sig=rol;
for i=1:length(files)
    im1=[];
    for n=1:4
        im1(:,:,n)=imread(files{i},n);
    end
    rol{i}=findrolonies(im1,20);
    sig{i}=readintensities(rol{i},im1,[1 1 1;1 1 1;1 1 1]);%read intensities using the original rol0
    
end
rol0=rol;

save allrolonies;

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
%% basecall based on intensity in a 3x3 area.
    %basecall
seq=zeros(length(rol{1}),length(files));
scores=seq;
for i=1:length(files)
    [~,seq(:,i)]=max(sig{i}(I(:,i),:),[],2);
    scores(:,i)=max(sig{i}(I(:,i),:),[],2)./sqrt(sum(sig{i}(I(:,i),:).^2,2));
end




%check the bleedthrough profile:
figure; scatter(sig{1}(:,1),sig{1}(:,2),'.')
figure; scatter(sig{1}(:,3),sig{1}(:,4),'.')

%convert seq

seqC(seq==1)='G';
seqC(seq==2)='T';
seqC(seq==3)='A';
seqC(seq==4)='C';
seqC=reshape(seqC,[],length(files));

%mark all non-TATT sequences (non SNAP25), check manually on images
snap25=ismember(seqC,'TATT','rows');
for i=1:length(files)
    im1=[];
    for n=1:4
        im1(:,:,n)=imread(files{i},n);
        figure;imshow(im1(:,:,n),[0 100]);
        hold on; scatter(rol{1}(snap25==0,2),rol{1}(snap25==0,1),'+r');
    end
    
end

%basecall using image, not rolonies
seq1=zeros(length(rol{1}),length(files));
scores1=seq1;
sig1=cell(4,1);

for i=1:length(files)
    im1=[];
    for n=1:4
        im1(:,:,n)=imread(files{i},n);
    end
    sig1{i}=readintensities(rol{1},im1,[1 1 1;1 1 1;1 1 1]);%read intensities using the original rol0
   
    [~,seq1(:,i)]=max(sig1{i},[],2);
    scores1(:,i)=max(sig1{i},[],2)./sqrt(sum(sig1{i}.^2,2));
    
end
seqC1(seq1==1)='G';
seqC1(seq1==2)='T';
seqC1(seq1==3)='A';
seqC1(seq1==4)='C';
seqC1=reshape(seqC1,[],length(files));

%mark all non-TATT sequences (non SNAP25), check manually on images
snap251=ismember(seqC1,'TATT','rows');
%compare the two seq qualities
figure;histogram(reshape(scores,1,[]));
set(gca,'xlim',[0.5 1]);
xlabel('Quality score');
ylabel('# of rolonies');
figure;histogram(reshape(scores1,1,[]));


figure;scatter(reshape(scores,1,[]),reshape(scores1,1,[]),'.');
set(gca,'xlim',[0.5 1],'ylim',[0.5 1]);

figure;scatter(reshape(scores1,1,[]),reshape(mindist,1,[]),'.');
set(gca,'xlim',[0.5 1]);
xlabel('Im basecall qual');
ylabel('matching rolony dist');

figure;scatter(reshape(scores,1,[]),reshape(mindist,1,[]),'.');
set(gca,'xlim',[0.5 1]);
xlabel('Im basecall qual');
ylabel('matching rolony dist');

%plot rolonies over images
files=dir(fullfile('RGB*'));
    files={files.name};
    figure;
for i=1:4
    im=imread(files{i});
    for n=1:length(rol0{i})
        im(rol0{i}(n,1),rol0{i}(n,2),:)=255;
    end
    imwrite(im,['rol',files{i}]);
end

%% check ICP quality
figure; scatter(rol0{1}(:,2),rol0{1}(:,1),'b+');
hold on;
scatter(rol0{2}(:,2),rol0{2}(:,1),'kx');
scatter(rol0{3}(:,2),rol0{3}(:,1),'ro');
scatter(rol0{4}(:,2),rol0{4}(:,1),'ks');
set(gca,'xlim',[250 350],'ylim',[250 350]);

figure; scatter(rol{1}(:,2),rol{1}(:,1),'b+');
hold on;
scatter(rol{2}(:,2),rol{2}(:,1),'kx');
scatter(rol{3}(:,2),rol{3}(:,1),'ro');
scatter(rol{4}(:,2),rol{4}(:,1),'ks');
set(gca,'xlim',[250 350],'ylim',[250 350]);

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

figure; scatter(rolnr{1}(:,2),rolnr{1}(:,1),'b+');
hold on;
scatter(rolnr{2}(:,2),rolnr{2}(:,1),'kx');
scatter(rolnr{3}(:,2),rolnr{3}(:,1),'ro');
scatter(rolnr{4}(:,2),rolnr{4}(:,1),'ks');
set(gca,'xlim',[250 350],'ylim',[250 350]);


save temp3
