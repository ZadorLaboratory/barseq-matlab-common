function [A,Y]=sparsenmf(B,k,b,labels)

option=struct;
option.iter=400;
option.dis=false;
option.beta=b;

W=[];H=[];spr=[];residual=[];

for m=1:10
[W(:,:,m),H(:,:,m),~,spr(m),residual(m)]=sparsenmfnnls(B,k,option);
end
[res,idx]=min(residual+sqrt(option.beta)*spr);
A=W(:,:,idx);
Y=H(:,:,idx);
[~,N]=max(A);
[~,i]=sort(N);
A=A(:,i);
Y=Y(i,:);

figure;imagesc(A);colormap('gray');
%set(gca,'ytick',1:11,'yticklabels',labels(2:end),'xtick',1:k,'xticklabels',1:k);
xlabel('Modules');
ylabel('Projections');
colorbar;
res=residual(idx)
%Ynorm=Y./repmat(mean(Y),k,1);
[R,P]=corr(Y');
figure;imagesc(R.*((P*(k-1)*k/2)<0.05),[-0.5 0.5]);
set(gca,'ytick',1:k,'xtick',1:k);
colorbar;

end