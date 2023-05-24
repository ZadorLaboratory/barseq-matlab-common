function [primer,padlock,Tm,loc]=designdnapadlocksnoRT(sequence,barcode)
% design RT primers and padlock probes. This version has higher TM
% requirement and try to maximize the number of padlocks per gene, but allows
%a flexible number of overlap between rt primer and padlock defined by overlap.
% Rules:
% skip the last 20nt.
% padlocks: <90 bases total, each arm 21-23 bases, TM~60 but >50. GC content for each arm: 40-60%
% Inter-probe distance: >=5 bases
%%
backbone='GATCGTCGGACTGTAGAACTCTGAACCTGTCG';
backbone2='CCT';
primer={};
padlock={};
Tm=[];
loc=[];
overlap=0;

if length(sequence)>100
    startidx=length(sequence)-25;
    lastpos=startidx+100;
    %find first primer/padlock.
    count=1;
    
    for idx1=startidx:-1:60
        idx=idx1+25;
        if (lastpos-idx) >= 3
            for n=1
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
                if min(max(qual))>=50 %check if both arms have Tm >=50
                    I(1)=find(qual(:,1)>=50,1,'last');
                    I(2)=find(qual(:,2)>=50,1,'last');
                    %[~,I]=find(qual>=55); %Find indices of padlocks with highest Tm
                    
                    %check complexity of arms
                    if padlockcomplexity([lower(arm5{I(1)}),lower(arm3{I(2)})])>0.01 %filter out low-complexity padlock sequences
                        padlock{count}=[lower(arm3{I(2)}),backbone,barcode,backbone2,lower(arm5{I(1)})]; %find the padocks with tm closest to 60
                        Tm(1,count)=-1;
                        Tm(3,count)=qual(I(1),1);
                        Tm(2,count)=qual(I(2),2);
                        loc(count)=idx-25-n-overlap;
                        lastpos=idx-25-n-overlap-23+I(1); %most 5' position of padlock
                        primer{count}='NNNN';
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

