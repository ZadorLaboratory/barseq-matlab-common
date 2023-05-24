function [width]=findhalfmaximalwidth(peaks,input);

k=1;
for i=1:size(peaks,2)
    for j=1:size(peaks(i).position,2)
        if peaks(i).position(j)~=23; %ignore bulb peaks. they have no width;
        halfmax=input(i,peaks(i).position(j))/2;
        
        %find left half maximal width
        if peaks(i).position(j)==1;
            lefthalf=0;
        else
            lefthalf=peaks(i).position(j)-find(input(i,1:peaks(i).position(j))<halfmax,1,'last');
        end
        if isempty(lefthalf)==1
            lefthalf=0;
        end
        
        
        %find right half maximal width
        if peaks(i).position(j)==size(input,2);
            righthalf=0;
        else
            righthalf=find(input(i,peaks(i).position(j):end)<halfmax,1,'first')-1;
        end
        if isempty(righthalf)==1
            righthalf=0;
        end
        
        width(k)=lefthalf+righthalf;
        k=k+1;
        end
    end
   
end
