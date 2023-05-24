function mi=calculatemi(x,y)
%for discrete data x and y,find the mutual information between them. x and
%y should be the same height.
%%
uniqx=unique(x,'rows');
uniqy=unique(y,'rows');

%%
Px=zeros(size(uniqx,1),1);
for i=1:size(uniqx,1)
    Px(i)=sum(ismember(x,uniqx(i,:),'rows'))/size(x,1);
end

Py=zeros(size(uniqy,1),1);
for i=1:size(uniqy,1)
    Py(i)=sum(ismember(y,uniqy(i,:),'rows'))/size(y,1);
end

%%
Hx=-Px'*log2(Px);



%%

Hxgiveny=[];
for i=1:size(uniqx,1)
    for m=1:size(uniqy,1)
        Pxy=sum(ismember(x,uniqx(i,:),'rows')&ismember(y,uniqy(m,:),'rows'))/size(x,1);
        Hxgiveny(i,m)=-Pxy*log2(Pxy/Py(m));
    end
end
Hxgiveny(isnan(Hxgiveny))=0;
Hxgiveny=sum(Hxgiveny(:));
%%
mi=Hx-Hxgiveny;
end