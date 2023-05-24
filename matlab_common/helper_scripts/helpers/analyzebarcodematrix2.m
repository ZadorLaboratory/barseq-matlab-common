function [B,x]=analyzebarcodematrix2(barcodematrix,spikes)
% given barcodematrix, individualdatawithdistancealigned, spikes, produce
% B, CoPNnorm, AsPNnorm,filtaligneddepth,fisseqCoPNnorm,fisseqAsPNnorm 
%% threshold matrix to contain only high confidence projecitons
load(barcodematrix);
load(spikes);

lowercutoff=10;
uppercutoff=10000;
sourcecutoff=20;

    B=barcodematrixthresholded(max(barcodematrixthresholded(:,1:12),[],2)>lowercutoff & max(barcodematrixthresholded(:,1:12),[],2)<uppercutoff & barcodematrixthresholded(:,13)>sourcecutoff,:);
    seq=refbarcodesthresholded(max(barcodematrixthresholded(:,1:12),[],2)>lowercutoff & max(barcodematrixthresholded(:,1:12),[],2)<uppercutoff & barcodematrixthresholded(:,13)>sourcecutoff,:);
x=[];for i=1:13;x(i)=length(spikes(i).counts2u);end

end
