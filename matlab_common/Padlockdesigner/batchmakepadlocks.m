function [allprimers,allpadlocks,allTm,allloc,allgenes]=batchmakepadlocks(accnum, barcodes,startidx,usecds,idx)
% given NCBI accession numbers accnum and the barcodes. make padlocks for
% all genes. accnum should be n x 1 cell arrays, with each cell containing
% an accession number. barcodes should be n x m char arrays for m digit
% barcodes.
%startidx is the first barcode to be used. usecds indicate whether to use
%cds. idx is a logical array indicating the genes in accnum to design.

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
        try
            cds=sequence(data.CDS.indices(1):data.CDS.indices(2));
        catch
            warning(['No cds found in sequence ',accnum{n},', using whole sequence\n']);
            err=2;
        end
        if err~=1
            if usecds==1
                if err==2
                    [primer,padlock,Tm,loc]=designdnapadlocks(sequence,barcodes(n+startidx-1,:));
                else
                    [primer,padlock,Tm,loc]=designdnapadlocks(cds,barcodes(n+startidx-1,:));
                    %if length(loc)<5
                    %    [primer,padlock,Tm,loc]=designdnapadlocks(sequence,barcodes(n+startidx-1,:));
                    %end
                end
            else 
                [primer,padlock,Tm,loc]=designdnapadlocks(sequence,barcodes(n+startidx-1,:));
            end


            allprimers=[allprimers;primer'];
            allpadlocks=[allpadlocks;padlock'];
            allTm=[allTm;Tm'];
            allloc=[allloc;loc'];
            allgenes=[allgenes;(repmat({data.Definition,['GII',num2str(n+startidx-1)]},length(loc),1)),cellstr([repmat([data.LocusName,'-'],length(loc),1),num2str((1:length(loc))')])];
        end
    end
end
allgenes(:,3)=cellfun(@(x)x(x~=' '),allgenes(:,3),'UniformOutput',false);

end



