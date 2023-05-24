function cmap=hclrainbow(num,h,c,l)
%% generate hcl_rainbow from r in matlab
if ~exist('h','var')
    h=[0 360];
end
if ~exist('c','var')
    c=0.5;
end
if ~exist('l','var')
    l=0.7;
end


hcl=[linspace(h(1),h(2),num)',ones(num,1)*100*c,ones(num,1)*100*l];
lch=hcl(:,[3 2 1]);
lch=reshape(lch,1,size(lch,1),[]);
cmap=colorspace('lch->rgb',lch);
cmap=squeeze(cmap);

% %%
% figure;
% imagesc(0);
% colorbar;
% colormap(rgb)

