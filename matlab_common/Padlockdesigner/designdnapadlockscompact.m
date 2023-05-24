function [primer,padlock,Tm,loc]=designdnapadlockscompact(sequence,barcode)
% design RT primers and padlock probes
% Rules: 
% avoid first/last 100bp
% RT primers:26 bases, Tm>60, higher is better
% padlocks: <90 bases total, each arm 19-21 bases, TM~58 but >51. GC content for each arm: 40-60%
% Inter-probe distance: >=5 bases
%this version overlaps the 5' arm of padlocks with the last 12 bases of the RT primer

backbone='GATCGTCGGACTGTAGAACTCTGAACCTGTCG';
backbone2='CCT';
prime={};
padlock={};
Tm=[];
loc=[];

if length(sequence)>150
    startidx=length(sequence)-20;
    lastpos=startidx+100;
    %find first primer/padlock.
    count=1;
    
    for idx=startidx:-1:60
            primerprop=oligoprop(sequence(idx-25:idx));
             if primerprop.Tm(5)>62 && lastpos-idx > 75 ...
                    && isempty(strfind(lower(sequence(idx-71:idx)),'gggg')) ...
                    && isempty(strfind(lower(sequence(idx-71:idx)),'cccc'))
                arm5{1}=sequence(idx-20-32:idx-32);
                arm5{2}=sequence(idx-19-32:idx-32);
                arm5{3}=sequence(idx-18-32:idx-32);
                arm3{1}=sequence(idx+1-32:idx-32+21);
                arm3{2}=sequence(idx+1-32:idx-32+20);
                arm3{3}=sequence(idx+1-32:idx-32+19);
                qual=zeros(3,2);
                for n=1:3
                    arm5prop(n)=oligoprop(arm5{n});
                    arm3prop(n)=oligoprop(arm3{n});
                    qual(n,1)=(arm5prop(n).GC>=40&&arm5prop(n).GC<=60)*arm5prop(n).Tm(5);
                    qual(n,2)=(arm3prop(n).GC>=40&&arm3prop(n).GC<=60)*arm3prop(n).Tm(5);
                end
                if min(abs(qual(:,1)-58))<7 && min(abs(qual(:,2)-58))<7
                    padlock{count}=[lower(arm3{abs(qual(:,2)-58)==min(abs(qual(:,2)-58))}),backbone, ...
                        barcode,backbone2,lower(arm5{abs(qual(:,1)-58)==min(abs(qual(:,1)-58))})];
                    lastpos=idx;
                    Tm(1,count)=primerprop.Tm(5);
                    Tm(3,count)=qual(abs(qual(:,1)-58)==min(abs(qual(:,1)-58)),1);
                    Tm(2,count)=qual(abs(qual(:,2)-58)==min(abs(qual(:,2)-58)),2);
                    loc(count)=idx-47;
                    primer{count}=seqrcomplement(sequence(idx-25:idx));
                    count=count+1;
                end
            end

        end
    end
    
end
