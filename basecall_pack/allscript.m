%fix file names
for i=6:16
    rawfiles=dir(fullfile(['*seq',num2str(i,'%.2u'),'.tif']));
    [rawfiles,~]=sort_nat({rawfiles.name});
    for n=length(rawfiles):-1:1
        [~,I]=regexpi(rawfiles{n},'point ');
        idx=str2double(rawfiles{n}(I+1));
        if idx==8
            movefile(rawfiles{n},['XY point 1_seq',num2str(i,'%.2u'),'.tif']);
        else
            movefile(rawfiles{n},['XY point ',num2str(idx+1),'_seq',num2str(i,'%.2u'),'.tif']);
        end
    end
end
%normalize filenames
rawfiles=dir(fullfile(['*.tif']));
[rawfiles,~]=sort_nat({rawfiles.name});
for n=1:length(rawfiles)
    [~,idx1]=regexpi(rawfiles{n},'point ');
    idx2=regexpi(rawfiles{n},'seq');
    newfilename=rawfiles{n}([1:(idx1+1),idx2:end]);
    if ~strcmp(newfilename,rawfiles{n})
        movefile(rawfiles{n},newfilename);
    end
end

slicenames={}; %manually input slicenames
cyclenumber=16;
for n=1:length(slicenames)
    rawfiles=dir(fullfile(['*point ',num2str(n),'*.tif']));
    [rawfiles,~]=sort_nat({rawfiles.name});
    mkdir([slicenames{n},'\'])
    for m=1:length(rawfiles)
        movefile(rawfiles{m},[slicenames{n},'\',slicenames{n},'seq',num2str(m,'%02u'),'.tif']);
    end
end
for n=1:length(slicenames)
    old=cd(slicenames{n});
    rgbout('');
    mkdir('RGB/');
    movefile('RGB*.tif','RGB/');
    cd(old);
end


%process images and align
addpath(pwd);

for i=1:length(slicenames)
    old=cd(slicenames{i});
    Ultraviewseqalignment2('seq');
    cd(old);
end

for i=1:length(slicenames)
    old=cd([slicenames{i},'/aligned']);
    rgbout('');
    movefile('RGB*','../RGB');
end% alignment is pretty good for cell bodies, but not good for projection sites.

%% for injection sites, find rolonies and basecall based on aligned images.
for i=[1,3,5]
    old=cd([slicenames{i},'/aligned']);
    [rol,sig]=seqfindrolonies('seq01',[20]);%find rolonies from first cycle
    save('rolonies.mat','rol','sig');
    [seq,seqC,scores,sig1]=imgbasecall(rol{1},ones(3,3)); %basecall at 1st cycle rolonies
    %remove sequences with N and find unique sequences
    seqnonN=seqC(max(seqC=='N',[],2)==0,:);
    sig1nonN=sig1(max(seqC=='N',[],2)==0,:,:);
    scoresnonN=scores(max(seqC=='N',[],2)==0,:);
    rol1nonN=rol{1}(max(seqC=='N',[],2)==0,:);
    [uniqseq,~,J]=unique(seqnonN,'rows');
    occ=histc(J,1:size(uniqseq,1));
    [sortedocc,I]=sort(occ,'descend');
    % plot first 200 neurons
    files=dir(fullfile('*seq01.tif'));
    im=[];
    for n=1:4
        im(:,:,n)=imread(files.name,n);
    end
    im=mean(im,3);
    figure;imshow(im,[0 max((max(im)))]);hold on;
    for n=1:200
        scatter(rol{1}(ismember(seqC,uniqseq(I(n),:),'rows'),2),rol{1}(ismember(seqC,uniqseq(I(n),:),'rows'),1),'.');
        
    end
    title([slicenames{i}, 'first 200 neurons']);
    savefig('roloniesonimage.fig');
    
    save('seq.mat','seq','seqC','scores','seqnonN','sig1nonN','sig1','rol1nonN','occ','uniqseq');
    cd(old);
end


%% realign projection images
%process images
for i=[2,4,6 7 8]
    old=cd([slicenames{i},'/original']);
    Ultraviewfiximage('seq'); %fix bleedthrough
    cd(old);
end

%find rolonies
for i=[2,4,6,7,8]
    old=cd([slicenames{i},'/original']);
    [rol,sig]=seqfindrolonies('aligned',[30 30 30 20]);
    load rolonies.mat;
    [rolnr,tr,tt]=seqalignmenticp('aligned',rol); %align images using rigid icp
    save('rigidalignment.mat','rolnr','tr','tt');
    save('rolonies.mat','rol','sig');
    %[seq,seqC,scores]=nricpbasecall(rol,sig,2); %fine alignment (non-rigid ICP) and basecall
    [seq,seqC,scores,sig1]=imgbasecall(rol{1},ones(3,3));
    seqnonN=seqC(max(seqC=='N',[],2)==0,:);
    sig1nonN=sig1(max(seqC=='N',[],2)==0,:,:);
    scoresnonN=scores(max(seqC=='N',[],2)==0,:);
    rol1nonN=rol{1}(max(seqC=='N',[],2)==0,:);
    [uniqseq,~,J]=unique(seqnonN,'rows');
    occ=histc(J,1:size(uniqseq,1));
    [sortedocc,I]=sort(occ,'descend');
    figure;hold on;
    for i=1:200
        scatter(rol{1}(ismember(seqC,uniqseq(I(i),:),'rows'),2),rol{1}(ismember(seqC,uniqseq(I(i),:),'rows'),1),'.');
        
    end
    set(gca,'Ydir','reverse');
    set(gca,'xlim',[0 1406],'ylim',[0 1839]);
    
    save('seq.mat','seq','seqC','scores');
    cd(old);
end

%align images using rigid icp on the rolonies
for i=[2,4,6,7,8]
    old=cd([slicenames{i},'/original']);
    load rolonies.mat;
    
    mkdir('fixed/');
    movefile('fixed*.tif','fixed/');
    mkdir('aligned/');
    movefile('aligned*.tif','aligned/');
    movefile('RGB*.tif','RGB/');
    
    cd(old);
end


%script for aligning/basecalling rolonies for all files in a single folder
%find rolonies
    files=dir(fullfile('fixed*'));
    files=sort_nat({files.name});
    rol=cell(length(files),1);
    sig=rol;
    for i=1:length(files)
        im1=[];
        for n=1:4
            im1(:,:,n)=imread(files{i},n);
        end
        rol{i}=findrolonies(im1,10+10*(i==1)+10*(i<=3));
        sig{i}=readintensities(rol{i},im1,[1 1 1;1 1 1;1 1 1]);%read intensities using the original rol0

    end
    rol0=rol;
%visually check rolony positions across cycles.
    figure;
    for i=1:length(files)
        if i<=16
            subplot(4,4,i);
            scatter(rol{i}(:,2),rol{i}(:,1),'.');
            set(gca,'xlim',[0 size(im1,2)],'ylim',[0 size(im1,1)]);
            title(['cycle',num2str(i)]);
        end
    end

save allrolonies;
% non-rigid ICP to align images
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

%check transformed rolonies
figure;
for i=1:length(files)
    if i<=16
        subplot(4,4,i);
        scatter(rolnr{i}(:,2),rolnr{i}(:,1),'.');
        %set(gca,'xlim',[0 size(im1,2)],'ylim',[0 size(im1,1)]);
        title(['cycle',num2str(i)]);
    end
end


%find matching rolonies
I=zeros(length(rol{1}),length(files));
mindist=I;
I(:,1)=(1:length(rol{1}))';
for i=2:length(files)
    d=pdist2(rolnr{1},rolnr{i});
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

%convert seq
seqC(seq==1)='G';
seqC(seq==2)='T';
seqC(seq==3)='A';
seqC(seq==4)='C';
seqC=reshape(seqC,[],length(files));
%impose a mindist threshold of 5
seqC(mindist>=3)='N';

figure;scatter(reshape(scores,1,[]),reshape(mindist,1,[]),'.');
set(gca,'xlim',[0.5 1]);
xlabel('Im basecall qual');
ylabel('matching rolony dist');

%plot rolonies with >=5 identified rolonies in the first seven bases
figure; scatter(rolnr{1}(:,2),rolnr{1}(:,1),'.');
hold on; scatter(rolnr{1}(sum(mindist(:,1:3)<3,2)>=3,2),rolnr{1}(sum(mindist(:,1:3)<3,2)>=3,1),'o');
figure; scatter(rolnr{1}(:,2),rolnr{1}(:,1),'.');
hold on; scatter(rolnr{1}(sum(mindist(:,1:7)<3,2)>=5,2),rolnr{1}(sum(mindist(:,1:7)<3,2)>=5,1),'o');
