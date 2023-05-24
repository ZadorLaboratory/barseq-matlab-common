function c=padlockcomplexity(seq2)
% test padlock sequence complexity using a 15-mer sliding window with 10-nt
% step size. Usually a value < 1e-3 should be filtered out, but < 1e-2 can
% also be used resulting in more stringent filter and few extra padlocks
% being filtered out.

c1=ones(floor(length(seq2)/10),1);
for n=1:floor(length(seq2)/10)
    seq3=seq2(((n-1)*10+1):min((n-1)*10+16,length(seq2)));
    u=lingseqcomplexity(seq3);
    c1(n)=prod(u);
end
c=prod(c1);