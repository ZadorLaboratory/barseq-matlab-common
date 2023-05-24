function finalistn=reduceprobes(finalist,nprobes)
%pick nprobes from the finalist for each gene while trying to spread the
%probes along the whole sequence. If number of found probes is smaller than
%nprobes, then pick all.

ugii=unique(finalist(:,2));
finalistn=cell(0);
for i=1:length(ugii)
    idx=find(ismember(finalist(:,2),ugii(i)));
    if length(idx)<=nprobes
        finalistn=[finalistn;finalist(idx,:)];
    else 
        idx1=idx(floor((1:nprobes)/nprobes*length(idx)));
        finalistn=[finalistn;finalist(idx1,:)];
    end
end

