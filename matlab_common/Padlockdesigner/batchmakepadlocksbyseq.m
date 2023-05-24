function [allprimers,allpadlocks,allTm,allloc,allgenes]=batchmakepadlocksbyseq(genenames,sequences, barcodes,giiidx)
% given NCBI accession numbers accnum and the barcodes. make padlocks for
% all genes. accnum should be n x 1 cell arrays, with each cell containing
% an accession number. barcodes should be n x m char arrays for m digit
% barcodes.
%giiidx is an array of indices indicating which gii to use. idx is a logical array indicating the genes in accnum to design.
%This version uses whole mRNA sequence and try to maximize number of
%padlock probes, but allows flexible number of overlap between RT and padlock to 5nts.

%no RT primers

allprimers={};
allpadlocks={};
allTm=[];
allloc=[];
allgenes={};

if iscell(barcodes)
    barcodes=char(barcodes);
end



for n=1:size(sequences,1)
    [primer,padlock,Tm,loc]=designdnapadlocksnoRT(sequences{n},barcodes(giiidx(n),:));
    allprimers=[allprimers;primer'];
    allpadlocks=[allpadlocks;padlock'];
    allTm=[allTm;Tm'];
    allloc=[allloc;loc'];
   allgenes=[allgenes;(repmat({genenames{n},['GII',num2str(giiidx(n))]},length(loc),1)),cellstr([repmat([genenames{n},'-'],length(loc),1),num2str((1:length(loc))')])];
end

end



