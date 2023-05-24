
function [cellid,bc,bc_count]=expand_bc(filt_neurons,cellid)
%given filt_neurons, concatenate all_bc
if ~exist('cellid','var')
    cellid=filt_neurons.id;
end
filt_neurons.total_bc_count=cellfun(@sum,filt_neurons.all_bc_count);

has_bc=filt_neurons.total_bc_count>0;

idx=ismember(filt_neurons.id,cellid);
bc=cell2mat(filt_neurons.all_bc(idx&has_bc));
try
    bc_count=cell2mat(filt_neurons.all_bc_count(idx&has_bc)')';
catch
    bc_count=cell2mat(filt_neurons.all_bc_count(idx&has_bc));
end
    
cellid=cellfun(@(x,y) repmat(x,numel(y),1), ...
     num2cell(filt_neurons.id(idx&has_bc)), ...
     filt_neurons.all_bc_count(idx&has_bc), ...
     'UniformOutput',0);
cellid=cell2mat(cellid);

end