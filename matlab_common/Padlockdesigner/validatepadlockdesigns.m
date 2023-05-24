function validtargetsnm32=validatepadlockdesigns(padlocks,names,bclength,blastoptfile)
%blast the primers/padlocks on NCBI against mouse genome + transcriptome,
%download the blast results and read them. bclength is the length of the
%RPI used
%This version allows choosing different blast configuration files

if ~exist('blastoptfile','var')
    blastoptfile='NMonly.asn';
end


%write arms into a fasta file
names=cellfun(@(x)x(x~=' '),names,'UniformOutput',false);

armidx=cellfun(@(x)regexpi(x,['GATCGTCGGACTGTAGAACTCTGAACCTGTCG',repmat('.',1,bclength),'CCT']),padlocks);
arms=cell(length(padlocks),1);
for i=1:length(padlocks)
    arms{i}=[padlocks{i}(armidx(i)+42:end),padlocks{i}(1:armidx(i)-1)];
end
armlengths=cellfun(@length,arms);
ligationpos=armlengths-armidx+1;%length of the 3' arm.

if exist('arms1.fasta','file')
    delete arms1.fasta;
end
fastawrite('arms1.fasta',names,arms);

%blast arm sequences, request btop and starting base of alignment
eval(['!blastn -query arms1.fasta -out results.txt -import_search_strategy ',blastoptfile,' -remote -outfmt "6 qseqid sseqid qstart btop evalue"']);

%!blastn -query test.fasta -out testresults.txt -import_search_strategy template.asn -remote -outfmt "6 qseqid sseqid qstart btop evalue"


%parse result file into alignment string. In the alignment string, 0 is
%match, 1 is mismatch, 2 is a gap (query present, but no match in subject),
%10,20,30.... is an insertion of 1, 2 ,3 following the base.


while 1
    
    resultfound=1;
    try
        fileid=fopen('results.txt');
    catch
        resultfound=0;
        fprintf('Results not ready yet, retry in 30 sec ...\n');
        pause(30);
    end
    if resultfound==1
        break
    end
end


a=textscan(fileid,'%s','delimiter','\r\n');
fclose(fileid);
b=cellfun(@strsplit,a{1},'UniformOutput',false);
b=vertcat(b{:});
    %b(:,2) has subject acc numbers, b(:,3) is alignment start idx, b(:,5) is
    %btop.

    %make alignment string
alignstring=cell(size(b,1),1);
for n=1:size(b,1)
    [matchstart,matchend]=regexpi(b{n,4},'\d+');
    [misstart,misend]=regexpi(b{n,4},'[-atcg]+');
    if str2double(b{n,3})>1
        alignstring{n}=repmat(2,1,str2double(b{n,3})-1);
    else
        alignstring{n}=double([]);
    end
    for i=1:size(matchstart,2)
        %fill the matching part
        alignstring{n}=[alignstring{n} zeros(1,str2double(b{n,4}(matchstart(i):matchend(i))))];
        %fill the mismatch part
        if i<size(matchstart,2)
            for m=misstart(i):2:misend(i)-1
                if b{n,4}(m)=='-' %insertion
                    alignstring{n}(end)=alignstring{n}(end)+10;
                elseif b{n,4}(m+1)=='-' %gap
                    alignstring{n}=[alignstring{n},2];
                else %mismatch
                    alignstring{n}=[alignstring{n},1];
                end
            end
        end
    end
    if length(alignstring{n})<armlengths(ismember(names,b{n,1}))
        alignstring{n}=[alignstring{n},2*ones(1,armlengths(ismember(names,b{n,1}))-length(alignstring{n}))];
    end
    
end

%evaluate alignment strings. The rule for detection is: 
% 1. perfect match of 3nt on each side of the ligation junction
% 2. no gap/insertion within 7nt on each side of the ligation junction
% 3. The matching bases of each arm has a Tm of at least 50 (changed to 37)
validtargets=[b(:,1:2),num2cell(zeros(size(b,1),1))];


for n=1:size(b,1)

    if sum(alignstring{n}(ligationpos(ismember(names,b{n,1}))-2:ligationpos(ismember(names,b{n,1}))+3))==0 ... %#1
            && sum(alignstring{n}(ligationpos(ismember(names,b{n,1}))-6:ligationpos(ismember(names,b{n,1}))+7)>1)==0%#2
        arm1=arms{ismember(names,b{n,1})}(1:ligationpos(ismember(names,b{n,1})));
        align1=alignstring{n}(1:ligationpos(ismember(names,b{n,1})));
        arm2=arms{ismember(names,b{n,1})}(ligationpos(ismember(names,b{n,1}))+1:end);
        align2=alignstring{n}(ligationpos(ismember(names,b{n,1}))+1:end);
        if sum(align1==0) && sum(align2==0)
            prop1=oligoprop(arm1(align1==0));
            prop2=oligoprop(arm2(align2==0));
            if prop1.Tm(5)>37 && prop2.Tm(5)>37
                validtargets{n,3}=1;
            end
            
        end
    end
end

validtargets=validtargets(cellfun(@(x)x,validtargets(:,3))==1,1:3);

subaccnum=cellfun(@(x)strsplit(x,'|'),validtargets(:,2),'UniformOutput',0);
subaccnum=vertcat(subaccnum{:});

validtargets(:,2)=subaccnum(:,4);
save('validtargets.mat','validtargets');
%get the gene name for the NM subjects, output gene name and accession
%number
save('blastresults.mat','b');
%nm=cellfun(@(x)~isempty(regexp(x,'^NM', 'once')),validtargets(:,2));
%validtargetsnm32=validtargets(nm,:);
%%
validtargetsnm32=validtargets;
for i=1:size(validtargetsnm32,1)
    idx=find(validtargetsnm32{i,2}=='.');
    if ~isempty(idx)
        validtargetsnm32{i,2}=validtargetsnm32{i,2}(1:idx-1);
    end
end

uniqtargets=unique(validtargetsnm32(:,2));
uniqnames=uniqtargets;

t=zeros(length(uniqnames),1);
for n=1:size(uniqtargets,1)
    tic
    m=1;
    while m>0&&m<10 % give up if more than 10 tries
        m=0;
        try
            data=getgenbank(uniqtargets{n});
            uniqnames{n}=data.Definition;
        catch e
            fprintf('%s \n', e.identifier)
            pause(10)
            m=m+1;
        end
    end
    if rem(n,1)==0
        fprintf('%u of %u genes checked \n',n,size(uniqtargets,1));
    end
    t(n)=toc;
end

%%
[~,I]=ismember(validtargetsnm32(:,2),uniqtargets);
validtargetsnm32(:,3)=uniqnames(I);


end

        
