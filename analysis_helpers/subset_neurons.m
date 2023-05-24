
function filt_neurons_subset=subset_neurons(filt_neurons,cellid)
%given a list of cellid, subset filt_neurons
filt_neurons_subset=filt_neurons;
[~,idx]=ismember(cellid,filt_neurons.id);
if sum(idx==0)>0
    error('Some cellids are not found in filt_neuron. Are they the same dataset?\n')
end



fnames=fieldnames(filt_neurons);
cellnum=numel(filt_neurons.id);
for n=1:numel(fnames)
    v=filt_neurons.(fnames{n});
    s=size(v,1);
    if s==cellnum
        filt_neurons_subset.(fnames{n})=v(idx,:);
    end
end


end
