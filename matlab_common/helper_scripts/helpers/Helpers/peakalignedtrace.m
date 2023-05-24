function [averagetrace_s,traceerror_s,shuffledtrace_s,shuffledtraceerror_s,shuffledtrace_matrix]=peakalignedtrace(input,inputcounts,areas,plotting);

%% plot average peakaligned trace for cortex
input0_tmp=input(sum(input(:,areas),2)~=0,:);
inputcounts0_tmp=inputcounts(sum(input(:,areas),2)~=0,:);
input0=input0_tmp(max(inputcounts0_tmp(:,areas),[],2)>30,:);

[m,loc]=max(input0(:,areas),[],2);
%normalize to peak hight;
inputn=input0(:,areas)./repmat(m,1,size(areas,2));

matrix2=-1*ones(size(inputn,1),size(areas,2)*2-1);
for i=1:size(inputn,1);matrix2(i,(size(areas,2)+1)-loc(i):((size(areas,2)+1)-loc(i))+(size(areas,2)-1))=inputn(i,:);end

ms=[];
for j=1:100
    inputn_shuffle=[];
    for i=1:size(inputn,1);inputn_shuffle(i,:)=inputn(i,randperm(size(areas,2)));end
    [m,loc]=max(inputn_shuffle,[],2);
    matrix2s=-1*ones(size(inputn_shuffle,1),size(areas,2)*2-1);
    for i=1:size(inputn_shuffle,1);matrix2s(i,(size(areas,2)+1)-loc(i):((size(areas,2)+1)-loc(i))+(size(areas,2)-1))=inputn_shuffle(i,:);end
    for k=1:size(matrix2,2);ms(j,k)=mean(matrix2s(matrix2s(:,k)~=-1,k));end
end

averagetrace=[];
for i=1:size(matrix2,2);averagetrace(i)=mean(matrix2(matrix2(:,i)~=-1,i));end



%make plot

x=-6300:300:6600;
ind=isnan(m);

%plot decay only in one direction
matrix3=[matrix2(:,ceil(size(matrix2,2)/2):size(matrix2,2));matrix2(:,ceil(size(matrix2,2)/2):-1:1)];


averagetrace_s=[];for i=1:(size(areas,2));averagetrace_s(i)=mean(matrix3(matrix3(:,i)~=-1,i));end


matrix3_s=[];


for j=1:100
inputn_shuffle=[];
for i=1:size(inputn,1);inputn_shuffle(i,:)=inputn(i,randperm(size(areas,2)));end
[m,loc]=max(inputn_shuffle,[],2);
matrix2s=-1*ones(size(inputn_shuffle,1),size(areas,2)*2-1);
for i=1:size(inputn_shuffle,1);matrix2s(i,(size(areas,2)+1)-loc(i):((size(areas,2)+1)-loc(i))+(size(areas,2)-1))=inputn_shuffle(i,:);end
matrix3_s=[matrix2s(:,size(areas,2):2*size(areas,2)-1);matrix2s(:,size(areas,2):-1:1)];
for k=1:size(areas,2);shuffledtrace_matrix(j,k)=mean(matrix3_s(matrix3_s(:,k)~=-1,k));end
end

%make plot


m=[];
m=averagetrace_s;
traceerror_s=0.*averagetrace_s;
x=0:300:max(size(areas,2)-1)*300;
ind1=isnan(averagetrace_s);

shuffledtrace_s=mean(shuffledtrace_matrix);
shuffledtraceerror_s=std(shuffledtrace_matrix);
x=0:300:max(size(areas,2)-1)*300;
ind2=isnan(shuffledtrace_s);

if plotting==1
figure;boundedline(x(~ind1),averagetrace_s(~ind1),traceerror_s(~ind1),'b','alpha');
hold on;

boundedline(x(~ind2),shuffledtrace_s(~ind2),shuffledtraceerror_s(~ind2),'r','alpha');
axis([0 (size(areas,2)-2)*300 0 1])
xlabel('distance from peak (um)')
ylabel('relative barcode abundance')

% figure;plot(matrix3');





end
% matrix4=matrix2;
% matrix4(matrix4==-1)=0;
% clustergram(matrix4,'Symmetric','false','Standardize',3,'Cluster',1,'ColumnPDist','cosine')
% 
% 
