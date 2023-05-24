function [p1,p2]=proj2bias(mitralthreshn)
%calculate the double projection bias compared to independent projection.
%the outputs p1 and p2 counts 1 - the probability that the double proj 
%probability is smaller or larger than that predicted by random proj.  
mitralprob=sum(mitralthreshn>0)./size(mitralthreshn,1);
mitralprob2=zeros(size(mitralthreshn,2));
for i=1:size(mitralthreshn,2)
    for j=1:size(mitralthreshn,2)
        mitralprob2(i,j)=sum((mitralthreshn(:,i).*mitralthreshn(:,j))>0)./size(mitralthreshn,1);
    end
end
%shuffling data for expected prob2
mitralprob2rand=zeros(size(mitralthreshn,2),size(mitralthreshn,2),10000);
for k=1:10000
    mitralrand=mitralthreshn;
    for n=1:size(mitralthreshn,2)
        mitralrand(:,n)=mitralthreshn(randperm(length(mitralthreshn)),n);
    end
    for i=1:size(mitralrand,2)
        for j=1:size(mitralrand,2)
            mitralprob2rand(i,j,k)=sum((mitralrand(:,i).*mitralrand(:,j))>0)./size(mitralrand,1);
        end
    end
end
p1=sum(repmat(mitralprob2,1,1,10000)>=mitralprob2rand,3)/10000;
p2=sum(repmat(mitralprob2,1,1,10000)<=mitralprob2rand,3)/10000;
