function psf=mmfindpsf(filename,l)
%estimate four psf from rolony images of four channels
if ~exist('l')
    l=11;
end
psfi=ones(l);

files=dir(['*',filename,'*.tif']);
files=sort_nat({files.name});
files=files{1};

psf=cell(4,1);
parfor i=1:4
    [~,psf{i}]=deconvblind(imread(files,i),psfi,15);
end
save('psf.mat','psf');
end

