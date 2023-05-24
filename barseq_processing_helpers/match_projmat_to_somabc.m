

function match_projmat_to_somabc(mismatch,projmat_fname,bc_info_fname)
curr_dir=cd('..');
if ~exist('projmat_fname','var')
    f=dir('*BarcodeMatrix.mat');
    f={f.name};
    projmat_fname=f{1}; % if bc_fname isn't provided, use the last alldata file that contains bc
end
if ~exist('bc_info_fname','var')
    f=dir('*sampleinfo.csv');
    f={f.name};
    bc_info_fname=f{1}; % if bc_fname isn't provided, use the last alldata file that contains bc
end
if ~exist('mismatch','var')
    mismatch=1; %allowing 1 mismatch by default
end
cd(curr_dir);
S=load('filt_neurons-bc-somas.mat','filt_neurons');
filt_neurons=S.filt_neurons;

T=load(fullfile('..',projmat_fname));
U=readtable(fullfile('..',bc_info_fname));

filt_neurons.proj_area=U.sample_name;
filt_neurons.proj_brain=U.brain;

%make projection matrix.
Bseq=convert_seq(char(T.Bseq));
Bseq=Bseq(:,1:size(filt_neurons.soma_bc,2));

%check if mismatch number is fine for the projection matrix
d_projmat=squareform(pdist(Bseq,'hamming'));
d_projmat=(d_projmat+eye(size(d_projmat)))*size(Bseq,2);
d_projmat=min(d_projmat);
figure('Position',[50 50 600 300]);
subplot(1,2,1)
histogram(d_projmat,-0.5:size(Bseq,2)+0.5);
title('min hamming dist of projmat')
if sum(d_projmat==mismatch)>0
    warning('%u barcodes in MAPseq cannot be distinguished with the allowed mismatch.\n',sum(d_projmat<=mismatch));
end

%check if mismatch number is fine for somas
d_somas=squareform(pdist(filt_neurons.soma_bc(filt_neurons.is_barcoded,:),'hamming'));
d_somas=(d_somas+eye(size(d_somas)))*size(filt_neurons.soma_bc,2);
d_somas=min(d_somas);
subplot(1,2,2)
histogram(d_somas,-0.5:size(filt_neurons.soma_bc,2)+0.5);
title('min hamming dist of somas')
if sum(d_somas==mismatch)>0
    warning('%u barcodes in somas cannot be distinguished with the allowed mismatch.\n',sum(d_somas<=mismatch));
end
exportgraphics(gcf,'min_hamming_dist.pdf','ContentType','vector');

%match barcodes
filt_neurons.projmat=zeros(numel(filt_neurons.id),size(T.B,2));
filt_neurons.projmat_norm=filt_neurons.projmat;
has_match=zeros(size(filt_neurons.id));
filt_neurons.duplicate_bc=has_match;

is_barcoded_idx=find(filt_neurons.is_barcoded);

d=pdist2(filt_neurons.soma_bc(is_barcoded_idx,:),Bseq,'hamming')*size(filt_neurons.soma_bc,2);
[min_d,idx]=min(d,[],2); %index in projmat to assign to each filt_neuron.
has_match=min_d<=mismatch&sum(d==min_d,2)==1;
[uniq_idx,~,ic]=unique(idx(has_match));
c=histcounts(ic,0.5:max(ic)+0.5);

filt_neurons.has_proj=zeros(numel(filt_neurons.id),1);
filt_neurons.has_proj(is_barcoded_idx)=has_match;

filt_neurons.duplicate_bc(is_barcoded_idx(has_match))=ismember(idx(has_match),uniq_idx(c>1));

filt_neurons.projmat(is_barcoded_idx(has_match),:)=T.B(idx(has_match),:);
filt_neurons.projmat_norm(is_barcoded_idx(has_match),:)=T.Bnorm(idx(has_match),:);


save('filt_neurons-bc-somas-proj.mat','filt_neurons');


end
