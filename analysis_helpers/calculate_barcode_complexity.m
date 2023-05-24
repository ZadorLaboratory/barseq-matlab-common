function filt_neurons=calculate_barcode_complexity(filt_neurons)
% calculate complexity of all_bc from filt_neuron, add field bc_complexity

for i=1:length(filt_neurons.id)

    r=size(filt_neurons.all_bc{i,1},1);
    if r==0  
       filt_neurons.bc_complexity{i}=0;
    else
        complexity=zeros(r,1);
        for j=1:r
            bc=char(48+filt_neurons.all_bc{i, 1}(j,:));
            comp=log10(prod(lingseqcomplexity(char(bc))));
            complexity(j)=comp;
        end
        filt_neurons.bc_complexity{i}=complexity;
    end
end
