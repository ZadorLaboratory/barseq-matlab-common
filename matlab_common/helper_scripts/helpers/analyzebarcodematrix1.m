function [B,x,filtaligneddepth,fisseqAllPN,LocAllPN]=analyzebarcodematrix1(barcodematrix,fisseqfile,spikes)
% given barcodematrix, individualdatawithdistancealigned, spikes, produce
% B, CoPNnorm, AsPNnorm,filtaligneddepth,fisseqCoPNnorm,fisseqAsPNnorm 
%% threshold matrix to contain only high confidence projecitons
load(barcodematrix);
load(spikes);

lowercutoff=10;
uppercutoff=10000;
sourcecutoff=20;

if ~strcmp(fisseqfile,'none')
    B=barcodematrixthresholded(max(barcodematrixthresholded(:,1:12),[],2)>lowercutoff & max(barcodematrixthresholded(:,1:12),[],2)<uppercutoff,1:12);
    seq=refbarcodesthresholded(max(barcodematrixthresholded(:,1:12),[],2)>lowercutoff & max(barcodematrixthresholded(:,1:12),[],2)<uppercutoff,:);
x=[];for i=1:12;x(i)=length(spikes(i).counts2u);end
else
    B=barcodematrixthresholded(max(barcodematrixthresholded(:,1:12),[],2)>lowercutoff & max(barcodematrixthresholded(:,1:12),[],2)<uppercutoff & barcodematrixthresholded(:,13)>sourcecutoff,1:12);
    seq=refbarcodesthresholded(max(barcodematrixthresholded(:,1:12),[],2)>lowercutoff & max(barcodematrixthresholded(:,1:12),[],2)<uppercutoff & barcodematrixthresholded(:,13)>sourcecutoff,:);
x=[];for i=1:12;x(i)=length(spikes(i).counts2u);end
end

%labels={'OB','OFC','Motor','Rstr','SSctx','Cstr','Amyg','IpVis','CVis','CAud','Thal','Tect'};



%Bnorm is normalized to spikes
Bnorm=B./repmat(x,length(B),1);
%C=ones(size(Bnorm,1),1);
Bnorm(:,13)=1;
Bnorm(:,14)=1;
%% hamming distance
seq15=seq(:,1:15);
[uniqueseq15,ia,ib]=unique(seq15,'rows');
% what is the hamming distance profile for these barcodes?
%for the full length reference
dfull=minhammingdist(seq,seq);
%for first15 only
d15=minhammingdist(uniqueseq15,uniqueseq15);

figure;
histogram(d15,-0.5:9.5,'FaceColor',[0 0.5 0.5]);xlabel('Min hamming distance');ylabel('Counts');
%% match in situ data to projections, with <2 mismatch.
if ~strcmp(fisseqfile,'none')
   load(fisseqfile);

    a=permute(squeeze(struct2cell(individualdatawithdistancealigned)),[2 1]);

        %[~,loc]=ismember(a(:,7),cellstr(char(seq15)));
        %[~,loccorr]=ismember(a(:,13),cellstr(char(seq15)));
        [dbetween,pos,counts]=minhammingdist_withhits(seq15,int8(cell2mat(a(:,13))));
        pos=pos.*(dbetween<2);
        pos(pos==0)=size(Bnorm,1)+1;
      Bnorm(end+1,:)=0;
     %   B(end+1,:)=0;

     %   quality=cell2mat(a(:,14));
     %   minquality=min(quality,[],2);
     %   qualityuncorr=cell2mat(a(:,8));
     %   minqualityuncorr=min(qualityuncorr,[],2);

    a(:,25)=num2cell(cellfun(@min,(a(:,14))));
    a(:,26:39)=num2cell(Bnorm(pos,:));
    [~,idx]=max(Bnorm(pos,1:12),[],2);
    a(:,40)=num2cell(idx.*(max(Bnorm(pos,1:12),[],2)));


    aligneddepth=cell2mat(a(:,[24 23 26:40 21 20 18]));
    depthslice=cell2mat(a(:,1));
    aligneddepth(:,end+1)=str2num(depthslice(:,end-1:end))*14;
    aligneddepth=aligneddepth(aligneddepth(:,18)>1.5,:);% filter out cells below cortex
    projectaligneddepth=aligneddepth(aligneddepth(:,16)~=0,:);% with matched MAPseq barcodes
figure;    histogram(((projectaligneddepth(:,19)+projectaligneddepth(:,20))));
a=input('\nInput lowerbound:');
b=input('\nInput upperbound:');

   filtaligneddepth=projectaligneddepth((projectaligneddepth(:,19)+projectaligneddepth(:,20))>a&(projectaligneddepth(:,19)+projectaligneddepth(:,20))<b,:);

fisseqAllPN=filtaligneddepth(:,3:14).*repmat(x,length(filtaligneddepth),1);
LocAllPN=[filtaligneddepth(:,[1,2,19,20,21])];



else
    filtaligneddepth=[];fisseqAllPN=[];LocAllPN=[];
end
end
