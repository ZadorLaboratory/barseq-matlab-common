% filter barcodes according to low and high abundance cutoffs and normalize
% by efficiency and across rows. then plot heatmap and clustergram

function [normalizedmatrix,filteredmatrix,filteredrefbarcodes,normalizedmatrix_tmp,normalizedmatrix_tmp_Aonly]=filterandnormalize(inputmatrix,actin,refbarcodes,areas,injsite,spikes,lowcutoff,highcutoff,plotting,plottingorder,clusters,plottitle);



%remove low abundance barcodes and high spikes

filteredmatrix=inputmatrix(max(inputmatrix(:,areas),[],2)>lowcutoff & max(inputmatrix(:,areas),[],2)<highcutoff,:);

%adjust refbarcodes accordingly
filteredrefbarcodes=refbarcodes(max(inputmatrix(:,areas),[],2)>lowcutoff & max(inputmatrix(:,areas),[],2)<highcutoff,:);


%plot max(targetsites) vs inj. site for all and filtered barcodes
if plotting==1
figure;
loglog(max(inputmatrix(:,areas),[],2),inputmatrix(:,injsite),'o');
hold on; loglog(max(filteredmatrix(:,areas),[],2),filteredmatrix(:,injsite),'o','Color','r');
xlabel('max(targetsites)');ylabel('inj. site');title([plottitle,'targetsite to inj. site correspondence']);
legend('removed','kept')
end


%normalize barcode matrix by efficiency
for i=1:size(filteredmatrix,2);sp(i)=size(spikes(i).counts2u,1);end
normalizedmatrix_tmp=filteredmatrix./repmat(sp,size(filteredmatrix,1),1);
% normalizedmatrix=normalizedmatrix_tmp./repmat(sum(normalizedmatrix_tmp(:,areas),2),1,size(normalizedmatrix_tmp,2));


%try normalizing by actin too
actinnormfactor=2.^(actin-actin(1));
normalizedmatrix_tmp_A=normalizedmatrix_tmp.*repmat(actinnormfactor,size(normalizedmatrix_tmp,1),1);

normalizedmatrix_tmp_Aonly=filteredmatrix.*repmat(actinnormfactor,size(filteredmatrix,1),1);


%now normalize across columns
normalizedmatrix=normalizedmatrix_tmp_A(:,areas)./repmat(sum(normalizedmatrix_tmp_A(:,areas),2),1,23);



%plot barcodes sorted by max
if plotting==1;
figure;
load ~/Documents/connectome/LCprojectome/cryocut_samples/ZL060/seq/nomoleculecountthreshold/greycolormap.mat
[m,loc]=max(normalizedmatrix(:,plottingorder),[],2);[s,ix]=sort(loc);colormap(cm);imagesc(normalizedmatrix(ix,plottingorder));
xlabel('areas')
ylabel('barcodes')
title(plottitle);
end



%plot clustergram
if clusters==1;
columnnames={'slice1','slice2','slice3','slice4','slice5','slice6','slice7','slice8','slice9','slice10','slice11','slice12','slice13','slice14','slice15','slice16',...
	'slice17','slice18','slice19','slice20','slice21','slice22','OBipsi'};
clustergram(normalizedmatrix,'ColumnLabels',columnnames,'Standardize',3,...
    'Symmetric','false','Colormap','Hot','LogTrans','false',...
    'OptimalLeafOrder','false','RowPDist','correlation','Cluster',1);
end