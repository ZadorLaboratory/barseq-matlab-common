function [cellid1,seq1,seqC1,sig1,score1]=mmbasecallbccells_xc1(filename,nucprofile,c)
%using segmented images, call barcodes from cells based on INTENSITY. If a nucprofile is
%provided, then reduce nuclear background. nucprofile should be a 1x4
%vector that gives the relative strengths of nuc background signals in the
%four channels.

folders=get_folders();
for i=1:numel(folders)
    %load cell mask
    original_dir=cd(fullfile(folders{i},'aligned'));
    S=load('cellmask.mat');
    maski=S.maski;
    [cellid1,seq1,seqC1,sig1,score1]=mmbasecallcells('bcseq',maski);
    save('bcsoma.mat','cellid1','seq1','seqC1','sig1','score1');
    cellid{i}=cellid1;
    seq{i}=seq1;
    seqC{i}=seqC;
    sig{i}=sig1;
    score{i}=score1;
    cd(original_dir)
end

save('all_bccells_intensity.mat','cellid','seq','seqC','sig','score');
