function [primer,padlock,Tm,loc]=designdnapadlocks2_custom(sequence,barcode,overlaprange,noprimer,backbone,backbone2,padlockTmtarget)
% design RT primers and padlock probes. This version has higher TM
% requirement and try to maximize the number of padlocks per gene, but allows
%a flexible number of overlap between rt primer and padlock defined by overlap.
% Rules: 
% skip the last 20nt.
% RT primers:26 bases, Tm>60, higher is better
% padlocks: <90 bases total, each arm 21-23 bases, TM~60 but >58. GC content for each arm: 40-60%
% Inter-probe distance: >=5 bases
% has the option of not using RT primers.
%%
if ~exist('backbone','var')
    backbone='GATCGTCGGACTGTAGAACTCTGAACCTGTCG';
end
if ~exist('backbone2','var')
    backbone2='CCT';
end
if ~exist('padlockTmtarget','var')
    padlockTmtarget=58;
end
if ~exist('overlaprange','var')
    overlaprange=[-5,5];
end
if numel(overlaprange)==1
    overlaprange=[-5,overlaprange];
end
if overlaprange(1)>overlaprange(2)
    overlaprange=overlaprange([2,1]);
end
%%

primer={};
padlock={};
Tm=[];
loc=[];

offset=23;%distance from RT primer end to the padlock ligation site

if ~exist('noprimer','var')
    useprimer=1;
elseif noprimer~=0
    useprimer=0;
elseif noprimer == 0
    useprimer=1;
end


if length(sequence)>150
    startidx=length(sequence)-20-(1-useprimer)*5;
    lastpos=startidx+100;
    %find first primer/padlock.
    count=1;
    
    for idx=startidx:-1:76
        primerprop=oligoprop(sequence(idx-25:idx));
        if useprimer==0||(primerprop.Tm(5)>=62 && (lastpos-idx) >= 3  && ~contains(lower(sequence(idx-25:idx)),'gggg')&& ~contains(lower(sequence(idx-25:idx)),'cccc')) % if primer Tm>62 and is atleast 5 bases from the closest sequence and no GGGG/CCCC
            for n=overlaprange(2):-1:overlaprange(1) %scan through all posible padlocks within the overlap range.
                %for each ligation site, check three padlock arms with
                %21-23 bases in arms
                arm5{1}=sequence(idx-25+n-22-offset:idx-25+n-offset); %idx-25-n-overlap is the position of the ligation site
                arm5{2}=sequence(idx-25+n-21-offset:idx-25+n-offset);
                arm5{3}=sequence(idx-25+n-20-offset:idx-25+n-offset);
                arm3{1}=sequence(idx-25+n+1-offset:idx-25+n+23-offset);
                arm3{2}=sequence(idx-25+n+1-offset:idx-25+n+22-offset);
                arm3{3}=sequence(idx-25+n+1-offset:idx-25+n+21-offset);
                qual=zeros(3,2);
                for i=1:3 %check GC content and check absence of GGGG/CCCC
                    arm5prop(i)=oligoprop(arm5{i});
                    arm3prop(i)=oligoprop(arm3{i});
                    qual(i,1)=(arm5prop(i).GC>=40&&arm5prop(i).GC<=60&&~contains(lower(arm5{i}),'gggg')&& ~contains(lower(arm5{i}),'cccc'))*arm5prop(i).Tm(5);
                    qual(i,2)=(arm3prop(i).GC>=40&&arm3prop(i).GC<=60&&~contains([lower(arm3{i}),'g'],'gggg')&& ~contains(lower(arm3{i}),'cccc'))*arm3prop(i).Tm(5);
                end
                if min(max(qual))>=padlockTmtarget %check if both arms have Tm >=58
                    I(1)=find(qual(:,1)>=padlockTmtarget,1,'last');
                    I(2)=find(qual(:,2)>=padlockTmtarget,1,'last');
                    %[~,I]=find(qual>=58); %Find indices of padlocks with highest Tm
                    
                    %check complexity of arms
                    if padlockcomplexity([lower(arm5{I(1)}),lower(arm3{I(2)})])>0.01 %filter out low-complexity padlock sequences
                        padlock{count}=[lower(arm3{I(2)}),backbone,barcode,backbone2,lower(arm5{I(1)})]; %find the padocks with tm closest to 60
                        Tm(1,count)=primerprop.Tm(5);
                        Tm(3,count)=qual(I(1),1);
                        Tm(2,count)=qual(I(2),2);
                        loc(count)=idx-25+n-offset;
                        lastpos=idx-25+n-offset-23+I(1); %most 5' position of padlock
                        if useprimer>0
                            primer{count}=seqrcomplement(sequence(idx-25:idx));
                        else
                            primer{count}='not_used';
                        end
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
