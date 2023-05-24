function match_axonal_bc(filt_neurons_fname,varargin)
%match axonal bc to soma_bc

%is the filt_neurons filename provided?
if ~exist('filt_neurons_fname','var')
    %load filt_neurons, need to have at least soma data
    try
        T=load('filt_neurons-bc-somas-proj.mat');
    catch
        T=load('filt_neurons-bc-somas.mat');
    end

else
    T=load(filt_neurons_fname);
end
filt_neurons=T.filt_neurons;

% parse thresholds
if ~isempty(varargin)
    varargin=reshape(varargin,2,[]);
    for i=1:size(varargin,2)
        switch lower(varargin{1,i})
            case 'mismatch'
                mismatch=varargin{2,i};
            case 'score'
                qual_thresh=varargin{2,i};
            case 'signal'
                sig_thresh=varargin{2,i};
            case 'bc_fname'
                bc_fname=varargin{2,i};
        end
    end
end


if ~exist('bc_fname','var')
    %load bcs
    f=dir('alldata*bc.mat');
    f={f.name};
    S=load(f{end});
else
    S=load(bc_fname);
end
bc=S.bc;

%assign default thresholds
if ~exist('mismatch','var')
    mismatch=1;
end
if ~exist('qual_thresh','var')
    qual_thresh=0.8;
end
if ~exist('sig_thresh','var')
    sig_thresh=0;% not using a sig thresh, but keeping this option open
end
%% Assign axonal/dendritic barcodes to neurons, match to barcode library, and assign subclass-level labels.
% make mock bc CCF positions. 

%initialize axonalbc fields
filt_neurons.axonalbc_id=cell(numel(filt_neurons.id),1);
filt_neurons.axonalbc_pos=cell(numel(filt_neurons.id),1);
filt_neurons.axonalbc_slice=cell(numel(filt_neurons.id),1);
filt_neurons.axonalbc_count=cell(numel(filt_neurons.id),1);



%match barcodes
barcoded_idx=find(filt_neurons.is_barcoded);
d=pdist2(filt_neurons.soma_bc(barcoded_idx,:),bc.seq,'hamming')*size(filt_neurons.soma_bc,2);
matched_bc=d<=mismatch;

% median (across cycles) of max (across channels) signal
% median of qual (across cycles)
%%
axonalbc_id={};
axonalbc_pos={};
axonalbc_slice={};
axonalbc_count={};
parfor i=1:numel(barcoded_idx) 
    %tic
    filt_matched_bc=find(matched_bc(i,:)'& ...
        median(max(bc.sig,[],3),2)>=sig_thresh & ...
        median(bc.qual,2)>=qual_thresh)';
    axonalbc_id{i}=filt_matched_bc;
    axonalbc_pos{i}=bc.pos(filt_matched_bc,:);
    axonalbc_slice{i}=bc.slice(filt_matched_bc);
    axonalbc_count{i}=numel(filt_matched_bc);
    %toc
end
filt_neurons.axonalbc_id(barcoded_idx)=axonalbc_id;
filt_neurons.axonalbc_pos(barcoded_idx,:)=axonalbc_pos;
filt_neurons.axonalbc_slice(barcoded_idx)=axonalbc_slice;
filt_neurons.axonalbc_count(barcoded_idx)=axonalbc_count;

%%
save('filt_neurons-bc-somas-proj-axonal.mat','filt_neurons','-v7.3')
end


