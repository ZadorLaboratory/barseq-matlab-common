
function filt_neurons_comb=combine_neurons(filt_neurons1,filt_neurons2)
%given a list of cellid, subset filt_neurons
filt_neurons_comb=struct;

fnames=fieldnames(filt_neurons1);
cellnum=numel(filt_neurons1.id);
for n=1:numel(fnames)
    v=filt_neurons1.(fnames{n});
    s=size(v,1);
    if s==cellnum
        filt_neurons_comb.(fnames{n})=[filt_neurons1.(fnames{n});filt_neurons2.(fnames{n})];
    else
        filt_neurons_comb.(fnames{n})=filt_neurons1.(fnames{n});
    end
end
%check if neuron ids are unique. If not, throw a warning, but continue
if max(ismember(filt_neurons1.id,filt_neurons2.id))>0
    warning('Some neurons are repeated in the two. Values from the second dataset will be replaced if in conflict.\n');
    cellid=unique(filt_neurons_comb.id);
    filt_neurons_comb=subset_neurons(filt_neurons_comb,cellid);
end




end