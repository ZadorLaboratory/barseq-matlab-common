
function plot_genes_on_fov(rolonies,gene_list,radius)
if ~exist('radius','var')
    radius=4;
end

imagesize=[3200 3200];

%make images of rolonies on each fov
output_dir=fullfile('..','figs',[gene_list{:}]);
mkdir(output_dir);
%gene_list={'RG(B19)'};
% make images with gene_list

uniq_FOV=unique(floor(rolonies.cellid/10000));
[~,gene_idx]=ismember(gene_list,rolonies.genes(:,1));
%colors=hclrainbow(numel(gene_idx)+1);
%colors=colors(1:numel(gene_idx),:);
colors=colormap('lines');
close
colors=colors(1:7,:);
colors=1-(1-colors)/2;%brighten the colors
if numel(gene_idx)>size(colors,1)
    colors=[colors;distinguishable_colors(numel(gene_idx)-size(colors,1),[colors;0 0 0;1 1 1])];
else
    colors=colors(1:numel(gene_idx),:);
end
folders=get_folders();
%%
for ii=1:numel(uniq_FOV)
    im=zeros([imagesize,numel(gene_idx)]);
    in_slice=floor(rolonies.cellid/10000)==uniq_FOV(ii);
    parfor n=1:numel(gene_idx)
        im(:,:,n)=full( ...
            sparse(rolonies.pos40x(in_slice&rolonies.id==gene_idx(n),2), ...
            rolonies.pos40x(in_slice&rolonies.id==gene_idx(n),1), ...
            1,imagesize(2),imagesize(1)) ...
            );
        im(:,:,n)=imdilate(im(:,:,n),offsetstrel('ball',radius,1));
        im(:,:,n)=im(:,:,n)-min(im(:,:,n),[],'all');
    end
    im=reshape(im,[],size(im,3));
    im_rgb=uint8(reshape(min(im*colors,1),imagesize(1),imagesize(2),[])*255);
    imwrite(im_rgb,fullfile(output_dir,[folders{ii},'.jpg']));
end
    
figure;%imshow(im_rgb);title(folders(ii),'Interpreter','none')
hold on;
for n=1:numel(gene_idx)
    scatter([],[],1,colors(n,:),'filled');
end
legend(gene_list, ...
    'FontSize',10, ...
    'Location','best')
axis off

exportgraphics(gcf,fullfile(output_dir,'legend.pdf'),'ContentType','vector')
close all;
end

