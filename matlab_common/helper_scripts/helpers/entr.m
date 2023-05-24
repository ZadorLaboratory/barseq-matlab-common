function S=entr(x,edges)
%calculate the entropy of x given bins defined by edges.
dx=discretize(x,edges);
px=zeros(length(edges)-1,1);
for i=1:length(edges)-1
    px(i)=sum(dx==i)/length(x);
end
%if probability=0, then adjust probability to 1 to avoid log2(0)
px(px==0)=1;
%calculate entropy
S=-(log2(px'))*px;
end
