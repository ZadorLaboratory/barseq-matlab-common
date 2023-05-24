function seqC=convert_seq(sequence)
% convert a sequence from 1,2,3,4 to G,T,A,C and vice versa
if ischar(sequence)
    seqC=zeros(size(sequence),'int8');
    seqC(sequence=='G')=1;
    seqC(sequence=='T')=2;
    seqC(sequence=='A')=3;
    seqC(sequence=='C')=4;
else
    seqC=char(sequence);
    seqC(sequence==1)='G';
    seqC(sequence==2)='T';
    seqC(sequence==3)='A';
    seqC(sequence==4)='C';
end
end