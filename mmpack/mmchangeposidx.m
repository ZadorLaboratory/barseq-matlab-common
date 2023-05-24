function mmchangeposidx(c)
%Add c to the Position number of all tiles generated micromanager. Runs in
%the processed folder in standard pipeline.
if ~exist('c','var')
    c=30;%the number to add to Pos.
end

folders10x=dir('MAX*');
folders10x(~[folders10x.isdir])=[];
folders10x=sort_nat({folders10x.name});
folders10xnew=folders10x;
for i=1:numel(folders10x)
    idx=str2num(folders10x{i}(8:end));
    folders10xnew{i}=[folders10x{i}(1:7),num2str(idx+c)];
    movefile(folders10x{i},folders10xnew{i});
end

cd 10x
filestiles=dir('stitched*.tif');
filestiles=sort_nat({filestiles.name});
filestilesnew=filestiles;
for i=1:numel(filestiles)
    idx=str2num(filestiles{i}(12:end-4));
    filestilesnew{i}=[filestiles{i}(1:11),num2str(idx+c),filestiles{i}(end-3:end)];
    movefile(filestiles{i},filestilesnew{i});
end
