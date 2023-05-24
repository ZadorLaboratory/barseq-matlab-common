function idx = maximizeentropy(proj1,depths1,binwidth,mingroupsize)
%binary split dataset using binarized proj1 to maximize information gain on
%depths1
%%
mind=min(depths1);
maxd=max(depths1);
idx=zeros(length(depths1),1);

m=1; %index for which layer
projidx=0;
typeidx=0;
if ~exist('mingroupsize','var')
    mingroupsize=20;
end

if range(depths1)>200
    c=histcounts(depths1,mind:binwidth:maxd);
    logc=log2(c/sum(c));
    logc(isnan(logc))=0;
    stotal=sum(-c.*logc)/length(depths1); %entropy
    scurr=stotal;
    
    while numel(depths1)>mingroupsize&&m<=4
        
        %split using each projection,calculate
        mins=zeros(size(proj1,2),1);
        I=mins;
        idxnew=repmat(idx,1,size(proj1,2));
        for i=1:size(proj1,2)
            
            idxnew01=idx;
            idxnew00=idxnew01;
            idxnew10=idxnew00;
            idxnew11=idxnew00;
            
            idxnew01(idx==0&proj1(:,i)>0)=1;
            idxnew00(idx==0&proj1(:,i)==0)=1;
            idxnew10(idx==1&proj1(:,i)>0)=0;
            idxnew11(idx==1&proj1(:,i)==0)=0;
            
            c011=histcounts(depths1(idxnew01==1),mind:binwidth:maxd);
            c010=histcounts(depths1(idxnew01==0),mind:binwidth:maxd);
            c001=histcounts(depths1(idxnew00==1),mind:binwidth:maxd);
            c000=histcounts(depths1(idxnew00==0),mind:binwidth:maxd);
            c111=histcounts(depths1(idxnew11==1),mind:binwidth:maxd);
            c110=histcounts(depths1(idxnew11==0),mind:binwidth:maxd);
            c101=histcounts(depths1(idxnew10==1),mind:binwidth:maxd);
            c100=histcounts(depths1(idxnew10==0),mind:binwidth:maxd);
            
            
            s01=(nansum(-c011.*log2(c011/sum(c011)))+nansum(-c010.*log2(c010/sum(c010))))/length(depths1);
            s00=(nansum(-c001.*log2(c001/sum(c001)))+nansum(-c000.*log2(c000/sum(c000))))/length(depths1);
            s11=(nansum(-c111.*log2(c111/sum(c111)))+nansum(-c110.*log2(c110/sum(c110))))/length(depths1);
            s10=(nansum(-c101.*log2(c101/sum(c101)))+nansum(-c100.*log2(c100/sum(c100))))/length(depths1);
            
            s01rand=zeros(100,1);
            s00rand=s01rand;
            s11rand=s01rand;
            s10rand=s01rand;
            %statistical test
            parfor n=1:100
                idxrand=randperm(length(depths1));
                c011=histcounts(depths1(idxnew01(idxrand)==1),mind:binwidth:maxd);
                c010=histcounts(depths1(idxnew01(idxrand)==0),mind:binwidth:maxd);
                c001=histcounts(depths1(idxnew00(idxrand)==1),mind:binwidth:maxd);
                c000=histcounts(depths1(idxnew00(idxrand)==0),mind:binwidth:maxd);
                c111=histcounts(depths1(idxnew11(idxrand)==1),mind:binwidth:maxd);
                c110=histcounts(depths1(idxnew11(idxrand)==0),mind:binwidth:maxd);
                c101=histcounts(depths1(idxnew10(idxrand)==1),mind:binwidth:maxd);
                c100=histcounts(depths1(idxnew10(idxrand)==0),mind:binwidth:maxd);
                
                s01rand(n)=(nansum(-c011.*log2(c011/sum(c011)))+nansum(-c010.*log2(c010/sum(c010))))/length(depths1);
                s00rand(n)=(nansum(-c001.*log2(c001/sum(c001)))+nansum(-c000.*log2(c000/sum(c000))))/length(depths1);
                s11rand(n)=(nansum(-c111.*log2(c111/sum(c111)))+nansum(-c110.*log2(c110/sum(c110))))/length(depths1);
                s10rand(n)=(nansum(-c101.*log2(c101/sum(c101)))+nansum(-c100.*log2(c100/sum(c100))))/length(depths1);
            end
            
            p01=sum(s01rand>s01);
            p00=sum(s00rand>s00);
            p11=sum(s11rand>s11);
            p10=sum(s10rand>s10);
            if p01<=95
                s01=stotal;
            end
            if p00<=95
                s00=stotal;
            end
            if p11<=95
                s11=stotal;
            end
            if p10<=95
                s10=stotal;
            end
            [mins(i),I(i)]=min([s01,s00,s11,s10]);
            switch I(i)
                case 1
                    idxnew(:,i)=idxnew01;
                case 2
                    idxnew(:,i)=idxnew00;
                case 3
                    idxnew(:,i)=idxnew10;
                case 4
                    idxnew(:,i)=idxnew11;
            end
            
        end
        
        %check if any of the splits result in lower entropy. If yes, update
        %current entropy and current idx. If not, stop.
        if min(mins)<scurr
            [scurr,projidx(m)]=min(mins);
            typeidx(m)=I(projidx(m));
            %update idx
            if min(sum(idx==0),sum(idx==1))<min(sum(idxnew(:,projidx(m))==0),sum(idxnew(:,projidx(m))==1))
                idx=idxnew(:,projidx(m));
            elseif min(sum(idxnew(:,projidx(m))==0),sum(idxnew(:,projidx(m))==1))>mingroupsize
                idx=idxnew(:,projidx(m));
            end
            
            m=m+1;
        else
            break
        end
        
        
        
    end
    if numel(unique(idx))>1
        idx=idx+1;
    end
end





end

