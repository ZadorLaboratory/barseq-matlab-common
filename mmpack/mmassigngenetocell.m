function [rolassignment,expmat,cellid]=mmassigngenetocell(lroi1,segim1,id221)
%assign rolonies to segmented image. rolony position is ins lroi as x,y.
%Produces rolony assignment and expression matrix

    rolassignment=segim1(sub2ind(size(segim1),lroi1(:,2),lroi1(:,1)));
    cellid=unique(segim1(:));
    cellid=cellid(cellid>0);
    uid=unique(id221);
    uid(uid==0)=[];
    
expmat=zeros(size(cellid,1),length(uid));
for n=1:size(cellid,1)
    expmat(n,:)=histcounts(id221(rolassignment==cellid(n),1),0.5:length(uid)+0.5);
end
end
