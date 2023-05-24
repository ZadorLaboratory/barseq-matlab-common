function [rolnr,seq,seqC,scores,I]=nricpbasecall(rol0,sig,thresh)
%non-rigid ICP to fine tune rolony alignment, then call bases.rol0 and sig
%are rolony positions and signal intensities generated by seqfindrolonies.
%thresh is the distance threshold allowed for rolonies between cycles 
% (typically set to 2 pixels for 20x images).
rolnr=rol0;
target.vertices=[rol0{1},ones(length(rol0{1}),1)];
target.faces=delaunay(rol0{1}(:,1),rol0{1}(:,2));
Options=struct;
Options.normalWeighting = 0;
for i=2:length(rol0)
    source.vertices=[rol0{i},ones(length(rol0{i}),1)];
    source.faces=delaunay(rol0{i}(:,1),rol0{i}(:,2));
    [pointsTransformed, ~] = nricp(source, target,Options);
    rolnr{i}=pointsTransformed(:,1:2);
end


%find matching rolonies
I=zeros(length(rol0{1}),length(rol0));
mindist=I;
I(:,1)=(1:length(rol0{1}))';
for i=2:length(rol0)
    for n=1:length(rolnr{1})
        [mindist(n,i),I(n,i)]=min(pdist2(rolnr{1}(n,:),rolnr{i}));
    end
    %d=pdist2(rolnr{1},rolnr{i});
    %[mindist(:,i),I(:,i)]=min(d,[],2);
    %mindist is the distance between 1st 
    %cycle rolony and nth cycle rolony. I is the index for the
    %corresponding rolony in nth channel.
end

%basecall
seq=zeros(length(rol0{1}),length(rol0));
scores=seq;
for i=1:length(rol0)
    [~,seq(:,i)]=max(sig{i}(I(:,i),:),[],2);
    scores(:,i)=max(sig{i}(I(:,i),:),[],2)./sqrt(sum(sig{i}(I(:,i),:).^2,2));
end

%convert seq
seqC(seq==1)='G';
seqC(seq==2)='T';
seqC(seq==3)='A';
seqC(seq==4)='C';
seqC=reshape(seqC,[],length(rol0));
%impose a mindist threshold 
seqC(mindist>thresh)='N';

figure;scatter(reshape(scores,1,[]),reshape(mindist,1,[]),'.');
set(gca,'xlim',[0.5 1]);
xlabel('Im basecall qual');
ylabel('matching rolony dist');
rolidx=I;
save('seq1.mat','seq','seqC','scores','rolnr','rolidx');
