function finalist=filteredvalidtargets(validtargets,padlocks,primers,allnames)
%filter out bad padlock/primers using blast results in validtargets
%get the gene abbreviation for each query frmo the blast results. If
%bypassing the blast search, then use {} for validtargets.
%%
names=allnames(:,3);
if ~isempty(validtargets)
    padlocklist=validtargets(:,1);
    for n=1:size(padlocklist,1)
        padlocklist{n}=padlocklist{n}(1:(regexpi(padlocklist{n},'-','once')-1));
    end
else
    padlocklist={};
end

%%
genelist=unique(padlocklist);
genename={};
if ~isempty(validtargets)
    for n=1:size(genelist,1)
        tic
        m=1;
        while m==1
            m=0;
            try
                idx=find(genelist{n}=='.');
                if idx>0
                    data=getgenbank(genelist{n}(1:idx-1));
                else
                    data=getgenbank(genelist{n});
                end
                genename{n}=data.Definition;
                idx1=regexpi(genename{n},'[()]');
                genename{n}=genename{n}(idx1(end-1)+1:idx1(end)-1);
            catch
                m=1;
            end
        end
        if rem(n,5)==0
            fprintf('%u of %u genes checked \n',n,size(genelist,1));
        end
        t(n)=toc;
    end
end
%%
%get the subject gene abbreviation from the blast results
if ~isempty(validtargets)
    sublist=unique(validtargets(:,3));

    for n=1:size(sublist,1)
        subgenelist{n}=sublist{n};
        idx1=regexpi(subgenelist{n},'[()]');
        if ~isempty(idx1)
            subgenelist{n}=subgenelist{n}(idx1(end-1)+1:idx1(end)-1);
        else
            subgenelist{n}=[];
        end
    end
end
%%
%match subjects to query. good padlocks are marked 1 in names(:,2)
for n=1:size(padlocks,1)
    names{n,2}=1;
    if ~isempty(validtargets)
        hits=validtargets(ismember(validtargets(:,1),names{n,1}),:);
        % for each padlock name, search for all hits, check if the hits are the
        % expected gene.
        for m=1:size(hits,1)
            hitquery=hits{m,1}(1:(regexpi(hits{m,1},'-','once')-1));
            hitqueryname=genename(ismember(genelist,hitquery));
            hitsubname=subgenelist(ismember(sublist,hits{m,3}));
            if ~ismember(hitqueryname,hitsubname)&&~isempty(hitsubname{1})
                names{n,2}=0;
                break
            end
        end
    end
end
%%
%generate gene list from padlock names. pick all padlocks for each.
for n=1:size(names,1)
    Accnum{n}=names{n,1}(1:(regexpi(names{n,1},'-','once')-1));
end
uniqAccnum=unique(Accnum);

finalist={};
for i=1:length(uniqAccnum)
    idx=find(ismember(Accnum,uniqAccnum{i}));
    
    idx(cell2mat(names(idx,2))==0)=[]; %remove idx corresponding to promiscuous padlocks
    %if length(idx)>=3
    %    finalist=[finalist;[allnames(idx(1:3),:) primers(idx(1:3)) padlocks(idx(1:3))]];
    %elseif ~isempty(idx) && length(idx)<3
    if ~isempty(idx)
        finalist=[finalist;[allnames(idx,:) primers(idx) padlocks(idx)]];
    elseif isempty(idx)
        finalist=[finalist;[{} {} {}]];
    end
end
[~,I]=sort((finalist(:,2)));
finalist=finalist(I,:);





