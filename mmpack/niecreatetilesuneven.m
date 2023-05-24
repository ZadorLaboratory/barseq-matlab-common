function niecreatetilesuneven(fname,imwidth,pixelsize,overlap)
% given four edge locations, create tile with interpolated z positions.
% Must provide positions in sets of 4s.
%require jsonlab
%create tiles from saved pos files.

%%
data=readmatrix(fname); %this is the x, y, z positions


tileconfig={};
poslist=zeros(0,3);
labellist={};
for i=1:floor(size(data,1)/4)
    %% find coordinates for each edge
    x=data((i-1)*4+1:i*4,1)*1000; %in mm?
    y=data((i-1)*4+1:i*4,2)*1000; %in mm?
    z=data((i-1)*4+1:i*4,3);

    %% calculate midpoint xy
    midpointx=range(x)/2+min(x);
    midpointy=range(y)/2+min(y);
    %% regress z slope and midpoint on xy
    z1=[(x-midpointx) (y-midpointy) ones(numel(x),1)]\z;
    zslopex=z1(1);
    zslopey=z1(2);
    midpointz=z1(3);
    %% calculate tile config
    tileconfig{i}=[ceil(range(x)/(imwidth*(1-overlap/100)*pixelsize)) ceil(range(y)/(imwidth*(1-overlap/100)*pixelsize))]+1; %make it slightly bigger just in case
    midpoint=tileconfig{i}/2-0.5;

    %%
    for n=1:tileconfig{i}(1)*tileconfig{i}(2)

        %% add grid positions
        [GRID_COL,GRID_ROW]=ind2sub(tileconfig{i},n);
        GRID_COL=GRID_COL-1;
        GRID_ROW=GRID_ROW-1;
        %% change LABEL
        LABEL=['Pos',num2str(i), ...
            '_',num2str(GRID_COL,'%.3u'), ...
            '_',num2str(GRID_ROW,'%.3u')];
        %% change XY positions
        %find the device of XYstage

        Yoffset=(GRID_ROW-midpoint(2))*imwidth*(1-overlap/100)*pixelsize;
        Xoffset=(GRID_COL-midpoint(1))*imwidth*(1-overlap/100)*pixelsize;
        Y=round(midpointy+Yoffset);
        X=round(midpointx+Xoffset);

        %% change Z positions
        %find the device of XYstage

        Zoffset=Yoffset*zslopey + Xoffset*zslopex;
        Z=midpointz+Zoffset;
        poslist=[poslist;X/1000,Y/1000,Z];
        labellist=[labellist;{LABEL}];
    end
end
writematrix(poslist,['tiled',fname],'Delimiter',';');

save(['Posinfo_',fname(1:end-4),'.mat'],'labellist','poslist','tileconfig');
