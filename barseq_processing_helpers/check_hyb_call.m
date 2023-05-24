
function check_hyb_call(slice)
% manually check the image of a few FOVs
load('../codebookhyb.mat','codebookhyb');
load('genehyb.mat','lroihyb','idhyb');
uniqhybid=unique(idhyb{slice});
figure;
hold on;
for i=1:numel(uniqhybid)
    scatter(lroihyb{slice}(idhyb{slice}==i,1),lroihyb{slice}(idhyb{slice}==i,2),3,'o','filled');
end
pbaspect([1 1 1]);
legend(codebookhyb(:,1),'Location','eastoutside');
set(gca,'ydir','reverse');
end