function [locall,locp,locnp]=findloc(barcodematrix,fisseqfile)
% given barcodematrix, individualdatawithdistancealigned, spikes, produce
% B, CoPNnorm, AsPNnorm,filtaligneddepth,fisseqCoPNnorm,fisseqAsPNnorm 
%% threshold matrix to contain only high confidence projecitons
load(barcodematrix);

lowercutoff=10;
uppercutoff=10000;
sourcecutoff=10;
    seq=refbarcodesthresholded(max(barcodematrixthresholded(:,1:12),[],2)>lowercutoff & max(barcodematrixthresholded(:,1:12),[],2)<uppercutoff,:);
seq15=seq(:,1:15);
%% match in situ data to projections, with <2 mismatch.
   load(fisseqfile);

    a=permute(squeeze(struct2cell(individualdatawithdistancealigned)),[2 1]);

        [dbetween,pos,~]=minhammingdist_withhits(seq15,int8(cell2mat(a(:,13))));
        pos=pos.*(dbetween<2);
 
     %   quality=cell2mat(a(:,14));
     %   minquality=min(quality,[],2);
     %   qualityuncorr=cell2mat(a(:,8));
     %   minqualityuncorr=min(qualityuncorr,[],2);
    aligneddepth=cell2mat(a(:,[18 20 22]));
    figure;    histogram(((aligneddepth(aligneddepth(:,3)>1.5&pos'>0,1)+aligneddepth(aligneddepth(:,3)>1.5&pos'>0,2))));
    a=input('\nInput lowerbound:');
    b=input('\nInput upperbound:');
locall=aligneddepth(((aligneddepth(:,1)+aligneddepth(:,2))>a&(aligneddepth(:,1)+aligneddepth(:,2))<b)&(aligneddepth(:,3)>1.5),1);
locp=aligneddepth((pos'>0)&((aligneddepth(:,1)+aligneddepth(:,2))>a&(aligneddepth(:,1)+aligneddepth(:,2))<b)&(aligneddepth(:,3)>1.5),1);
locnp=aligneddepth((pos'==0)&((aligneddepth(:,1)+aligneddepth(:,2))>a&(aligneddepth(:,1)+aligneddepth(:,2))<b)&(aligneddepth(:,3)>1.5),1);
end
