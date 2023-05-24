function [depths,angle]=mmcalculatedepth(contours,ftlroi,pixelsize)
% calculate distance from ROIs to both surfaces of the cortex. the roi
% coordinates ftlroi 

depths=cell(length(ftlroi),1);
angle=depths;
for n=1:length(ftlroi)
    [p1,d1,~]=distance2curve(contours(n).top,ftlroi{n});
    [p2,d2,~]=distance2curve(contours(n).bottom,ftlroi{n});
    depths{n}=[d1 d2]*pixelsize;
    u=[p1-ftlroi{n},zeros(size(p1,1),1)];
    v=[p2-ftlroi{n},zeros(size(p2,1),1)];
    a=zeros(size(u,1),1);
    for i=1:size(u,1)
        a(i)=atan2d(norm(cross(u(i,:),v(i,:))),dot(u(i,:),v(i,:)));
    end
    angle{n}=a;
end

save('depths.mat','depths','angle');
