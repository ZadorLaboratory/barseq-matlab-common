%Writes all starfish required metadata files. nissl images should be in the
%5th channel in the first round.


mkdir('starfish');
%% write experiment.json
fid=fopen('starfish/experiment.json','wt');
fprintf(fid, ...
'{\n   "version": "0.0.0",\n');

%manifest files
fprintf(fid,'   "images":{\n');
fprintf(fid, ...
    '        "primary": "primary_images.json",\n'); %primary images
fprintf(fid, ...
    '        "nissl": "nissl.json"\n'); %nissl
fprintf(fid,'   },\n');

%codebook
fprintf(fid,'    "codebook": "codebook.json",\n');

%finish
fprintf(fid,'}');
fclose(fid);

%% write maifest files
fid=fopen('starfish/primary_images.json','wt');
fprintf(fid, ...
'{\n   "version": "0.0.0",\n');

%FOV files
fprintf(fid,'   "contents":{\n');
fprintf(fid, ...
    '        "primaryfov_000": "primaryFOV_000.json"\n'); %primary images
fprintf(fid,'   },\n');

%finish
fprintf(fid,'}');
fclose(fid);

fid=fopen('starfish/nissl.json','wt');
fprintf(fid, ...
'{\n   "version": "0.0.0",\n');

%FOV files
fprintf(fid,'   "contents":{\n');
fprintf(fid, ...
    '        "nisslfov_000": "nisslFOV_000.json"\n'); %primary images
fprintf(fid,'   },\n');

%finish
fprintf(fid,'}');
fclose(fid);
%% write primary fov file
fid=fopen('starfish/primaryFOV_000.json','wt');


pixelsize=0.3225;
zspacing=3;

fprintf(fid, ...
'{\n   "version": "0.0.0",\n');

%dimensions
fprintf(fid, ...
'   "dimensions": [\n        "x",\n        "y",\n        "z",\n        "r",\n        "c"\n    ],\n');

%shapes
imfiles=dir('smalltiles/*.tif');
for i=1:length(imfiles)
    rtczxy(i,:)=cell2mat(textscan(imfiles(i).name,'seq%ualignedfixedT%uC%uZ%uX%uY%u'));
end
fprintf(fid, ...
'   "shape": [\n        "x":%u,\n        "y":%u,\n        "z":%u,\n        "r":%u,\n        "c":4\n    ],\n',max(rtczxy(:,[5 6 4 1])));


%file format
fprintf(fid,...
    '    "default_tile_format": "TIFF",\n');

%default shape
info=imfinfo(['smalltiles/',imfiles(1).name]);

fprintf(fid,...
    '    "default_tile_shape": [\n        %u,\n        %u\n    ],\n', info.Width,info.Height);

%tiles
fprintf(fid, '    "tiles": [\n');

firstentry=1;
for i=1:length(imfiles)
    info=imfinfo(['smalltiles/',imfiles(i).name]);
    if rtczxy(i,3)<=4 %c<=4 are sequencing channels
        %coordinates
        if firstentry==0
            fprintf(fid,',\n');
        end
        fprintf(fid, ...
            '        {\n           "coordinates":{\n                "x":[%.2f, %.2f],\n                "y":[%.2f, %.2f],\n                "z":[%.2f, %.2f],\n           },\n', ...
            (rtczxy(i,5)-1)*info.Width*pixelsize,(rtczxy(i,5))*info.Width*pixelsize, ...
            (rtczxy(i,6)-1)*info.Height*pixelsize,(rtczxy(i,6))*info.Height*pixelsize, ...
            (rtczxy(i,4)-1)*zspacing,(rtczxy(i,5))*zspacing);
        %indices
        fprintf(fid, ...
            '           "indices":{\n                "x":%u,\n                "y":%u,\n                "z":%u,\n                "r":%u,\n                "c":%u,\n           },\n', ...
            rtczxy(i,[5 6 4 1 3])-1);
        %filename
        fprintf(fid,'           "file": "%s",\n',imfiles(i).name);
        fprintf(fid,'           "tile_format": "TIFF",\n');
        %sha256
        [~,sha]=system(['CertUtil -hashfile "',info.Filename,'" SHA256']);
        idx=regexp(sha,newline);
        fprintf(fid,'           "sha256": "%s",\n',sha(idx(1)+1:idx(2)-1));
        %extras: tile size
        fprintf(fid,...
             '           "extras":{\n                "tile_shape": [\n                    %u,\n                    %u\n                ]\n            }\n', info.Width,info.Height);

        fprintf(fid,'        }');
        firstentry=0;
    end
end
fprintf(fid,'\n    ]\n');

%finish
fprintf(fid,'}');
fclose(fid);

%% write nissl fov file
fid=fopen('starfish/nisslFOV_000.json','wt');


pixelsize=0.3225;
zspacing=3;

fprintf(fid, ...
'{\n   "version": "0.0.0",\n');

%dimensions
fprintf(fid, ...
'   "dimensions": [\n        "x",\n        "y",\n        "z",\n        "r",\n        "c"\n    ],\n');

%shapes
imfiles=dir('smalltiles/*.tif');
for i=1:length(imfiles)
    rtczxy(i,:)=cell2mat(textscan(imfiles(i).name,'seq%ualignedfixedT%uC%uZ%uX%uY%u'));
end
fprintf(fid, ...
'   "shape": [\n        "x":%u,\n        "y":%u,\n        "z":%u,\n        "r":1,\n        "c":41\n    ],\n',max(rtczxy(:,[5 6 4])));


%file format
fprintf(fid,...
    '    "default_tile_format": "TIFF",\n');

%default shape
info=imfinfo(['smalltiles/',imfiles(1).name]);

fprintf(fid,...
    '    "default_tile_shape": [\n        %u,\n        %u\n    ],\n', info.Width,info.Height);

%tiles
fprintf(fid, '    "tiles": [\n');
firstentry=1;
for i=1:length(imfiles)
    info=imfinfo(['smalltiles/',imfiles(i).name]);
    if rtczxy(i,3)==5&&rtczxy(i,1)==1 %round 1 channel 5 is nissl
        %coordinates
        if firstentry==0
            fprintf(fid,',\n');
        end
        fprintf(fid, ...
            '        {\n           "coordinates":{\n                "x":[%.2f, %.2f],\n                "y":[%.2f, %.2f],\n                "z":[%.2f, %.2f],\n           },\n', ...
            (rtczxy(i,5)-1)*info.Width*pixelsize,(rtczxy(i,5))*info.Width*pixelsize, ...
            (rtczxy(i,6)-1)*info.Height*pixelsize,(rtczxy(i,6))*info.Height*pixelsize, ...
            (rtczxy(i,4)-1)*zspacing,(rtczxy(i,5))*zspacing);
        %indices
        fprintf(fid, ...
            '           "indices":{\n                "x":%u,\n                "y":%u,\n                "z":%u,\n                "r":%u,\n                "c":%u,\n           },\n', ...
            rtczxy(i,[5 6 4 1 3])-1);
        %filename
        fprintf(fid,'           "file": "%s",\n',imfiles(i).name);
        fprintf(fid,'           "tile_format": "TIFF",\n');
        %sha256
        [~,sha]=system(['CertUtil -hashfile "',info.Filename,'" SHA256']);
        idx=regexp(sha,newline);
        fprintf(fid,'           "sha256": "%s",\n',sha(idx(1)+1:idx(2)-1));
        %extras: tile size
        fprintf(fid,...
             '           "extras":{\n                "tile_shape": [\n                    %u,\n                    %u\n                ]\n            }\n', info.Width,info.Height);

        fprintf(fid,'        }');
        firstentry=0;
    end
end
fprintf(fid,'\n    ]\n');

%finish
fprintf(fid,'}');
fclose(fid);

%% write codebook
fid=fopen('starfish/codebook.json','wt');

codebook=cell(1);%copy codebook from excel file.

[genes,I,~]=unique(codebook(2:end,1));
codes=upper(char(codebook(I+1,4)));
codes1=uint8(codes);
codes1(codes=='G')=0;
codes1(codes=='T')=1;
codes1(codes=='A')=2;
codes1(codes=='C')=3;




fprintf(fid, ...
'{\n  "version": "0.0.0",\n');
%mapping
fprintf(fid,...
    '  "mappings": [\n');
for i=1:length(genes)
    %codeword
    fprintf(fid,'    {\n        "codeword": [\n');
    for n=1:size(codes1,2)
        if n==size(codes1,2)
            fprintf(fid,'            {"r": %u, "c": %u, "v": 1}\n',n-1,codes1(i,n));
        else
            fprintf(fid,'            {"r": %u, "c": %u, "v": 1},\n',n-1,codes1(i,n));
        end
    end
    fprintf(fid,'        ],\n        "target":"%s"\n',genes{i});
    if i==length(genes)
        fprintf(fid,'    }\n');
    else
        fprintf(fid,'    },\n');
    end
end
fprintf(fid, '  ]\n');
%finish
fprintf(fid,'}');

fclose(fid);    

    


