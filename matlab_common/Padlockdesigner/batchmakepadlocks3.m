function [allprimers,allpadlocks,allTm,allloc,allgenes]=batchmakepadlocks3(accnum, barcodes,giiidx,idx,overlap)
% given NCBI accession numbers accnum and the barcodes. make padlocks for
% all genes. accnum should be n x 1 cell arrays, with each cell containing
% an accession number. barcodes should be n x m char arrays for m digit
% barcodes.
%giiidx is an array of indices indicating which gii to use. idx is a logical array indicating the genes in accnum to design.
%This version uses whole mRNA sequence and try to maximize number of
%padlock probes, but allows flexible number of overlap between RT and padlock to 5nts.
%This version checks complexity

allprimers={};
allpadlocks={};
allTm=[];
allloc=[];
allgenes={};

if iscell(barcodes)
    barcodes=char(barcodes);
end



for n=1:length(accnum)
    if idx(n)
        err=0;
        try
            data=getgenbank(accnum{n});
            sequence=data.Sequence;
        catch
            warning(['Error reading sequence ',accnum{n},', skipping to next sequence\n']);
            err=1;
        end
        if err~=1
            
            [primer,padlock,Tm,loc]=designdnapadlocks2(sequence,barcodes(giiidx(n),:),overlap);
            allprimers=[allprimers;primer'];
            allpadlocks=[allpadlocks;padlock'];
            allTm=[allTm;Tm'];
            allloc=[allloc;loc'];
            allgenes=[allgenes;(repmat({data.Definition,['GII',num2str(giiidx(n))]},length(loc),1)),cellstr([repmat([data.LocusName,'-'],length(loc),1),num2str((1:length(loc))')])];
        end
    end
end
allgenes(:,3)=cellfun(@(x)x(x~=' '),allgenes(:,3),'UniformOutput',false);

end



