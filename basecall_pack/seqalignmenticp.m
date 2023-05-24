function [rolnr,tr,tt]=seqalignmenticp(filename,rol)
%Align tif images containing filename using rigid icp and write aligned 
%image files. rol is a list of rolony coordinates generated from seqfindrolonies
tt{1}=[0;0;0];
tr{1}=eye(3);
%rigid icp alignment
rolnr=rol;
for i=1:length(rol)
    [tr{i},tt{i}]=icp([rol{1}';ones(1,length(rol{1}))], ...
        [rol{i}';ones(1,length(rol{i}))],200);
    %align rolonies, currently only doing translation
    rolnr{i}=(tr{i}(1:2,1:2)*rol{i}'+repmat(tt{i}(1:2),1,length(rol{i})))';
end

%compare rolonies before alignment
figure;
for i=1:length(rol)
    if i<=16
        subplot(4,4,i);
        scatter(rol{i}(:,2),rol{i}(:,1),'.');
        set(gca,'xlim',[0 max(rol{1}(:,2)+5)],'ylim',[0 max(rol{1}(:,1)+5)]);
        title(['raw cycle',num2str(i)]);
    end
end

%compare rolonies after alignment
figure;
for i=1:length(rol)
    if i<=16
        subplot(4,4,i);
        scatter(rolnr{i}(:,2),rolnr{i}(:,1),'.');
        set(gca,'xlim',[0 max(rolnr{1}(:,2)+5)],'ylim',[0 max(rolnr{1}(:,1)+5)]);
        title(['aligned cycle',num2str(i)]);
    end
end


    





end
