function [idxall,silh,p]=clustall(x,dist)
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
            if pval<0.05/100 % 100 being the max number of splits required to produce clusters each with >=1% of all cells %/(size(currgroups,1)+sum(rem(currgroups,10)>0))
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
cmap(:,1)=1:-(1/63):0;
cmap(:,2)=1:-(0.5/63):0.5;
cmap(:,3)=1:-(0.5/63):0.5;

figure;imagesc(x(I,:));
colormap(cmap);

end
    