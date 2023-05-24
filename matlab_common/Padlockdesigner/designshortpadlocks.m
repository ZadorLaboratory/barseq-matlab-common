function [primer,padlock,Tm,loc]=designshortpadlocks(sequence,barcode)
% design RT primers and padlock probes
% Rules: 
% avoid first/last 100bp
% RT primers:26 bases, Tm>60, higher is better
% padlocks: <90 bases total, each arm 21-23 bases, TM~58 but >51. GC content for each arm: 40-60%
% Inter-probe distance: >=5 bases

backbone='GGACGGACTCGACCGTCG';
backbone2='';
prime={};
padlock={};
Tm=[];
loc=[];

    if length(sequence)>273
        startidx=length(sequence)-100;
        lastpos=startidx+100;
        %find first primer/padlock.
        count=1;

        for idx=startidx:-1:173
            primerprop=oligoprop(sequence(idx-25:idx));
             if primerprop.Tm(5)>62 && lastpos-idx > 75 ...
                    && isempty(strfind(lower(sequence(idx-71:idx)),'gggg')) ...
                    && isempty(strfind(lower(sequence(idx-71:idx)),'cccc'))
                arm5{1}=sequence(idx-15-43:idx-43);
                arm5{2}=sequence(idx-14-43:idx-43);
                arm5{3}=sequence(idx-13-43:idx-43);
                arm3{1}=sequence(idx+1-43:idx-43+16);
                arm3{2}=sequence(idx+1-43:idx-43+15);
                arm3{3}=sequence(idx+1-43:idx-43+14);
                qual=zeros(3,2);
                for n=1:3
                    arm5prop(n)=oligoprop(arm5{n});
                    arm3prop(n)=oligoprop(arm3{n});
                    qual(n,1)=(arm5prop(n).GC>=40&&arm5prop(n).GC<=60)*arm5prop(n).Tm(5);
                    qual(n,2)=(arm3prop(n).GC>=40&&arm3prop(n).GC<=60)*arm3prop(n).Tm(5);
                end
                if max(qual(:,1))>=50 &&  max(qual(:,2))>=50
                    padlock{count}=[lower(arm3{qual(:,2)==max(qual(:,2))}),backbone, ...
                        barcode,backbone2,lower(arm5{qual(:,1)==max(qual(:,1))})];
                    lastpos=idx;
                    Tm(1,count)=primerprop.Tm(5);
                    Tm(3,count)=qual(qual(:,1)==max(qual(:,1)),1);
                    Tm(2,count)=qual(qual(:,2)==max(qual(:,2)),2);
                    loc(count)=idx-47;
                    primer{count}=seqrcomplement(sequence(idx-25:idx));
                    count=count+1;
                end
            end

        end
    end
    
end
