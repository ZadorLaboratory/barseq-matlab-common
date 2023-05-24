function [primer,padlock,Tm,loc]=designdnapadlocks1(sequence,barcode)
% design RT primers and padlock probes. This version has higher TM
% requirement and try to maximize the number of padlocks per gene.
% Rules: 
% skip the last 20nt.
% RT primers:26 bases, Tm>60, higher is better
% padlocks: <90 bases total, each arm 21-23 bases, TM~60 but >58. GC content for each arm: 40-60%
% Inter-probe distance: >=5 bases

backbone='GATCGTCGGACTGTAGAACTCTGAACCTGTCG';
backbone2='CCT';
primer={};
padlock={};
Tm=[];
loc=[];

if length(sequence)>273
    startidx=length(sequence)-20;
    lastpos=startidx+100;
    %find first primer/padlock.
    count=1;
    
    for idx=startidx:-1:60
        primerprop=oligoprop(sequence(idx-25:idx));
        if primerprop.Tm(5)>=62 && (lastpos-idx) >= 5  && ~contains(lower(sequence(idx-25:idx)),'gggg')&& ~contains(lower(sequence(idx-25:idx)),'cccc')% if primer Tm>62 and is atleast 5 bases from the closest sequence and no GGGG/CCCC
            for n=1:min(26,idx-53) %scan through all posible padlocks from having 21 overlaps with primer to having 5nt gap.
                %for each ligation site, check three padlock arms with
                %21-23 bases in arms
                arm5{1}=sequence(idx-22-29-n:idx-29-n); %idx-29-n is the position of the ligation site
                arm5{2}=sequence(idx-21-29-n:idx-29-n);
                arm5{3}=sequence(idx-20-29-n:idx-29-n);
                arm3{1}=sequence(idx+1-29-n:idx-29-n+23);
                arm3{2}=sequence(idx+1-29-n:idx-29-n+22);
                arm3{3}=sequence(idx+1-29-n:idx-29-n+21);
                qual=zeros(3,2);
                for i=1:3 %check GC content and check absence of GGGG/CCCC
                    arm5prop(i)=oligoprop(arm5{i});
                    arm3prop(i)=oligoprop(arm3{i});
                    qual(i,1)=(arm5prop(i).GC>=40&&arm5prop(i).GC<=60&&~contains(lower(arm5{i}),'gggg')&& ~contains(lower(arm5{i}),'cccc'))*arm5prop(i).Tm(5);
                    qual(i,2)=(arm3prop(i).GC>=40&&arm3prop(i).GC<=60&&~contains(lower(arm5{i}),'gggg')&& ~contains(lower(arm5{i}),'cccc'))*arm3prop(i).Tm(5);
                end
                if min(max(qual))>=58 %check if both arms have Tm >=58
                    [~,I]=max(qual); %Find indices of padlocks with highest Tm
                    padlock{count}=[lower(arm3{I(2)}),backbone,barcode,backbone2,lower(arm5{I(1)})]; %find the padocks with tm closest to 60
                    Tm(1,count)=primerprop.Tm(5);
                    Tm(3,count)=qual(I(1),1);
                    Tm(2,count)=qual(I(2),2);
                    loc(count)=idx-29-n;
                    lastpos=idx-23-29-n+I(1); %most 5' position of padlock
                    primer{count}=seqrcomplement(sequence(idx-25:idx));
                    count=count+1;
                    break;
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
