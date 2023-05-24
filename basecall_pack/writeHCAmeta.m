cd('../point 2');
imfiles=dir('original/smalltiles/*.tif');
lowresfiles=dir('XY point 210xnissl/*.tif');


metadata=cell(length(imfiles)+length(lowresfiles),6);
for i=1:length(imfiles)
    metadata{i,1}=['point 2/point2/',imfiles(i).name];
    metadata{i,2}='420508-10-2';
    metadata{i,3}='tiff';
    info=imfinfo(['original/smalltiles/',imfiles(i).name]);
    [~,md5]=system(['CertUtil -hashfile "',info.Filename,'" MD5']);
    idx=regexp(md5,newline);
    metadata{i,4}=md5(idx(1)+1:idx(2)-1);
    metadata{i,5}='zador_singlerolseq_1';
    metadata{i,6}='zador_baristaseq_1';
end
for i=1:length(lowresfiles)
    metadata{i+length(imfiles),1}=['point 2/point2nissl10x/',lowresfiles(i).name];
    metadata{i+length(imfiles),2}='420508-10-2';
    metadata{i+length(imfiles),3}='tiff';
    info=imfinfo(['XY point 210xnissl/',lowresfiles(i).name]);
    [~,md5]=system(['CertUtil -hash file "',info.Filename,'" MD5']);
    idx=regexp(md5,newline);
    metadata{i+length(imfiles),4}=md5(idx(1)+1:idx(2)-1);
    metadata{i+length(imfiles),5}='zador_singlerolseq_1';
    metadata{i+length(imfiles),6}='zador_baristaseq_1';
end
    


%%
cd ('../point 3');  
imfiles=dir('original/smalltiles/*.tif');
lowresfiles=dir('XY point 310xnissl/*.tif');


metadata3=cell(length(imfiles)+length(lowresfiles),6);
for i=1:length(imfiles)
    metadata3{i,1}=['point 3/point3/',imfiles(i).name];
    metadata3{i,2}='420508-17-1';
    metadata3{i,3}='tiff';
    info=imfinfo(['original/smalltiles/',imfiles(i).name]);
    [~,md5]=system(['CertUtil -hashfile "',info.Filename,'" MD5']);
    idx=regexp(md5,newline);
    metadata3{i,4}=md5(idx(1)+1:idx(2)-1);
    metadata3{i,5}='zador_singlerolseq_1';
    metadata3{i,6}='zador_baristaseq_1';
end
for i=1:length(lowresfiles)
    metadata3{i+length(imfiles),1}=['point 3/point3nissl10x/',lowresfiles(i).name];
    metadata3{i+length(imfiles),2}='420508-17-1';
    metadata3{i+length(imfiles),3}='tiff';
    info=imfinfo(['XY point 310xnissl/',lowresfiles(i).name]);
    [~,md5]=system(['CertUtil -hashfile "',info.Filename,'" MD5']);
    idx=regexp(md5,newline);
    metadata3{i+length(imfiles),4}=md5(idx(1)+1:idx(2)-1);
    metadata3{i+length(imfiles),5}='zador_singlerolseq_1';
    metadata3{i+length(imfiles),6}='zador_baristaseq_1';
end

%%
cd ('../point 1 copy');  
imfiles=dir('original/smalltiles/*.tif');
lowresfiles=dir('XY point 110xnissl/*.tif');


metadata1=cell(length(imfiles)+length(lowresfiles),6);
for i=1:length(imfiles)
    metadata1{i,1}=['point 1/point1/',imfiles(i).name];
    metadata1{i,2}='420508-10-1';
    metadata1{i,3}='tiff';
    info=imfinfo(['original/smalltiles/',imfiles(i).name]);
    [~,md5]=system(['CertUtil -hashfile "',info.Filename,'" MD5']);
    idx=regexp(md5,newline);
    metadata1{i,4}=md5(idx(1)+1:idx(2)-1);
    metadata1{i,5}='zador_singlerolseq_1';
    metadata1{i,6}='zador_baristaseq_1';
end
for i=1:length(lowresfiles)
    metadata1{i+length(imfiles),1}=['point 1/point1nissl10x/',lowresfiles(i).name];
    metadata1{i+length(imfiles),2}='420508-17-1';
    metadata1{i+length(imfiles),3}='tiff';
    info=imfinfo(['XY point 110xnissl/',lowresfiles(i).name]);
    [~,md5]=system(['CertUtil -hashfile "',info.Filename,'" MD5']);
    idx=regexp(md5,newline);
    metadata1{i+length(imfiles),4}=md5(idx(1)+1:idx(2)-1);
    metadata1{i+length(imfiles),5}='zador_singlerolseq_1';
    metadata1{i+length(imfiles),6}='zador_baristaseq_1';
end
metadataall=[metadata1;metadata;metadata3];


T = cell2table(metadata(2:end,:),'VariableNames',metadata(1,:));
% Write the table to a CSV file
writetable(T,'../Imagefile.csv')