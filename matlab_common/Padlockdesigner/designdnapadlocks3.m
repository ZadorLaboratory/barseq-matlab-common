function [primer,padlock,Tm,loc]=designdnapadlocks3(sequence,barcode,overlap)
% design RT primers and padlock probes. This version has higher TM
% requirement and try to maximize the number of padlocks per gene, but allows
%a flexible number of overlap between rt primer and padlock defined by overlap.
% Rules: 
% skip the last 20nt.
% RT primers:26 bases, Tm>45, higher is better
% padlocks: <90 bases total, each arm 21-23 bases, TM~60 but >53. GC
% content for each arm: 40-60%, complexity >0.01
% Inter-probe distance: >=5 bases

backbone='GATCGTCGGACTGTAGAACTCTGAACCTGTCG';
backbone2='CCT';
primer={};
padlock={};
Tm=[];
loc=[];

if length(sequence)>60
    startidx=length(sequence)-20;
    lastpos=startidx+100;
    %find first primer/padlock.
    count=1;
    
    for idx=startidx:-1:60
        primerprop=oligoprop(sequence(idx-25:idx));
        if primerprop.Tm(5)>=45 && (lastpos-idx) >= 3  && ~contains(lower(sequence(idx-25:idx)),'gggg')&& ~contains(lower(sequence(idx-25:idx)),'cccc')% if primer Tm>62 and is atleast 5 bases from the closest sequence and no GGGG/CCCC
            for n=1:min(overlap+5,idx-53) %scan through all posible padlocks from having 21 overlaps with primer to having 5nt gap.
                %for each ligation site, check three padlock arms with
                %21-23 bases in arms
                arm5{1}=sequence(idx-25-n-overlap-22:idx-25-n-overlap); %idx-25-n-overlap is the position of the ligation site
                arm5{2}=sequence(idx-25-n-overlap-21:idx-25-n-overlap);
                arm5{3}=sequence(idx-25-n-overlap-20:idx-25-n-overlap);
                arm3{1}=sequence(idx-25-n-overlap+1:idx-25-n-overlap+23);
                arm3{2}=sequence(idx-25-n-overlap+1:idx-25-n-overlap+22);
                arm3{3}=sequence(idx-25-n-overlap+1:idx-25-n-overlap+21);
                qual=zeros(3,2);
                for i=1:3 %check GC content and check absence of GGGG/CCCC
                    arm5prop(i)=oligoprop(arm5{i});
                    arm3prop(i)=oligoprop(arm3{i});
                    qual(i,1)=(arm5prop(i).GC>=40&&arm5prop(i).GC<=60&&~contains(lower(arm5{i}),'gggg')&& ~contains(lower(arm5{i}),'cccc'))*arm5prop(i).Tm(5);
                    qual(i,2)=(arm3prop(i).GC>=40&&arm3prop(i).GC<=60&&~contains([lower(arm3{i}),'g'],'gggg')&& ~contains(lower(arm3{i}),'cccc'))*arm3prop(i).Tm(5);
                end
                if min(max(qual))>=53 %check if both arms have Tm >=58
                    I(1)=find(qual(:,1)>=53,1,'last');
                    I(2)=find(qual(:,2)>=53,1,'last');
                    %[~,I]=find(qual>=58); %Find indices of padlocks with highest Tm
                    
                    %check complexity of arms
                    if padlockcomplexity([lower(arm5{I(1)}),lower(arm3{I(2)})])>0.001 %filter out low-complexity padlock sequences
                        padlock{count}=[lower(arm3{I(2)}),backbone,barcode,backbone2,lower(arm5{I(1)})]; %find the padocks with tm closest to 60
                        Tm(1,count)=primerprop.Tm(5);
                        Tm(3,count)=qual(I(1),1);
                        Tm(2,count)=qual(I(2),2);
                        loc(count)=idx-25-n-overlap;
                        lastpos=idx-25-n-overlap-23+I(1); %most 5' position of padlock
                        primer{count}=seqrcomplement(sequence(idx-25:idx));
                        count=count+1;
                        break;
                    end
                end
            end
            
        end
    end
    
    if isempty(primer)
        primer={'NNNN'};
        padlock={'NNNN'};
        Tm=[0;0;0];
        loc=0;
    end
end
