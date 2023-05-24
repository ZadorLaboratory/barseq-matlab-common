
function basecall_hyb(hybthresh,hybbgn)
[folders,~]=get_folders();
    
load('../codebookhyb.mat','codebookhyb');
load('psf.mat','psf')
% read sequential signals.
idhyb={};
lroihyb={};
sighyb={};
parfor(i=1:length(folders),20)
    cd(folders{i});
    cd aligned
    [idhyb{i},lroihyb{i},sighyb{i}]=mmbasecallhyb(codebookhyb,hybthresh,hybbgn,psf);
    cd ../..
end
save('genehyb.mat','idhyb','lroihyb','sighyb','-v7.3');
end