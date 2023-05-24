function [seq,seqC,scores,sig]=imgbasecall(rol0,convmat)
%Basecall using aligned images. rol0 is the first cycle rolony positions (Nx2). 
files=dir(fullfile('aligned*'));
files=sort_nat({files.name});

%read out signals from all images
sig=zeros(size(rol0,1),length(files),4);
for i=1:length(files)
    im=[];
    imc=[];
    for n=1:4 %convolute images
        im(:,:,n)=imread(files{i},n);
        imc(:,:,n)=conv2(im(:,:,n),convmat,'same');
    end
    for m=1:size(rol0,1)
        sig(m,i,:)=imc(round(rol0(m,1)),round(rol0(m,2)),:);
    end
end

%basecall
[~,seq]=max(sig,[],3);
seq(sum(sig,3)==0)=0;
scores=max(sig,[],3)./sqrt(sum(sig.^2,3));

%convert seq
seqC(seq==1)='G';
seqC(seq==2)='T';
seqC(seq==3)='A';
seqC(seq==4)='C';
seqC(seq==0)='N';
seqC=reshape(seqC,[],length(files));
end
