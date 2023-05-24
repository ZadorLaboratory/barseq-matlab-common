function [idxall,silh,p]=clustall1(x,dist,labels)
%Split data using kmean into two groups each time, calculate significance. Test significance using sigClust. Repeat until no further groupings found. 


idxall=ones(size(x,1),1);

p=zeros(1,3);
n=1;

while 1
    currgroups=unique(idxall);
    if max(rem(currgroups,10))==0
        break
    end
   
    for i=1:size(currgroups,1)
        currgroups(i)
        if rem(currgroups(i),10)~=0&&sum(idxall==currgroups(i))>size(x,1)/100 % only cluster if the current subcluster was a new cluster and if current cluster size is larger than 1% of whole sample
            subx=x(idxall==currgroups(i),:);
            %subnormx=normx(idxall==currgroups(i),:);
            [idx,~,pval,~]=kmean2(subx,dist);
            %[idx,~,pval,~]=RFclust2(subx,100,max(round(size(subx,1)/20),2));
            p(n,1)=currgroups(i);
            p(n,2)=pval;
            p(n,3)=size(subx,1);
            n=n+1;
            if pval<0.05/length(currgroups)
                idxall(idxall==currgroups(i))=idxall(idxall==currgroups(i))*10+idx;
            else
                idxall(idxall==currgroups(i))=idxall(idxall==currgroups(i))*10;
            end
        else
            idxall(idxall==currgroups(i))=idxall(idxall==currgroups(i))*10;
        end
    end
end
silh=mean(silhouette(x,idxall));
[~,I]=sort(idxall);
cmap=ones(64,3);
cmap(:,2)=1:-(1/63):0;
cmap(:,3)=1:-(1/63):0;


figure;imagesc(x(I,:));
colormap(cmap);
if exist('labels','var')
    set(gca,'xtick',1:numel(labels),'xticklabelrotation',90,'xticklabel',labels);
end
end
    