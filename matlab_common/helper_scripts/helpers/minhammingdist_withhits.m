%calculate the minimum hamming distance between a set of reference
%sequences and test sequences (all in int8 format). perfect matches (i.e. same sequence in
%reference and test) will not be ignored. code returns 0 if there is only a
%perfect match.

function [D,pos,count]=minhammingdist_withhits(reference,test)
for i=1:size(test,1) %loop through test sequences
    for j=1:size(reference,2) %compare to all reference sequences
        A(:,j)=reference(:,j)-test(i,j); %do a columnwise substraction. substitutions are entries ~=0
    end
    substitutions=A~=0;
    totalsubs=sum(substitutions,2);
%    totalsubs(totalsubs==0)=size(substitutions,2)+1;
%     if sum(totalsubs==0)~=size(totalsubs,1)
    [D(i),pos(i)]=min(totalsubs);
    %count the number of min hamming distances
    count(i)=sum(totalsubs==D(i));
    %if ambivalent, then assign distance to 20.
    if sum(totalsubs==1)>1
    D(i)=20;
    end
        
    %     else
%     D(i)=0;
%     [tmp,pos(i)]=ismember(test(i,:),reference,'rows');
%     end
    
end
end

    