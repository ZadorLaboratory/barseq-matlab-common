function data2=mmcreatetileuneven(fname,imwidth,pixelsize,overlap)
% given four edge locations, create tile with interpolated z positions.
% Must provide positions in sets of 4s.
%require jsonlab
%create tiles from saved pos files.

%%
data=loadjson(fname);

data2=data;


%%
lastpos=0;
for i=1:floor(length(data.POSITIONS)/4)
    %% find coordinates for each edge
    x=[];
    y=[];
    z=[];
    for n=1:4
        for m=1:length(data.POSITIONS{(i-1)*4+n}.DEVICES)
            if contains(data.POSITIONS{(i-1)*4+n}.DEVICES{m}.DEVICE,'XYStage')
                y(n)=data.POSITIONS{(i-1)*4+n}.DEVICES{m}.Y;
                x(n)=data.POSITIONS{(i-1)*4+n}.DEVICES{m}.X;
            end
            if contains(data.POSITIONS{(i-1)*4+n}.DEVICES{m}.DEVICE,'ManualFocus')
                z(n)=data.POSITIONS{(i-1)*4+n}.DEVICES{m}.X;
            end
            
        end
    end
    %% calculate midpoint xy
    midpointx=range(x)/2+min(x);
    midpointy=range(y)/2+min(y);
    %% regress z slope and midpoint on xy
    z1=[(x-midpointx)' (y-midpointy)' ones(numel(x),1)]\z'; 
    zslopex=z1(1);
    zslopey=z1(2);
    midpointz=z1(3);
    %% calculate tile config
    tileconfig=[ceil(range(x)/(imwidth*(1-overlap/100)*pixelsize)) ceil(range(y)/(imwidth*(1-overlap/100)*pixelsize))]+1; %make it slightly bigger just in case
    midpoint=tileconfig/2-0.5;

    %% new position name
    posname=['Pos',num2str(i-1)];
    
    
    
    
    
    %% make new positions
    data2.POSITIONS(lastpos+(1:tileconfig(1)*tileconfig(2)))=data.POSITIONS(i*4-3);
    %%
    for n=1:tileconfig(1)*tileconfig(2)
        %% create data.Positions{}.properties
        data2.POSITIONS{lastpos+n}.PROPERTIES.OverlapUmX=num2str(imwidth*overlap*pixelsize/100);
        data2.POSITIONS{lastpos+n}.PROPERTIES.OverlapUmY=num2str(imwidth*overlap*pixelsize/100);
        data2.POSITIONS{lastpos+n}.PROPERTIES.OverlapPixelsX=num2str(imwidth*overlap/100);
        data2.POSITIONS{lastpos+n}.PROPERTIES.OverlapPixelsY=num2str(imwidth*overlap/100);
        %% add grid positions
        [data2.POSITIONS{lastpos+n}.GRID_COL,data2.POSITIONS{lastpos+n}.GRID_ROW]=ind2sub(tileconfig,n);
        data2.POSITIONS{lastpos+n}.GRID_COL=data2.POSITIONS{lastpos+n}.GRID_COL-1;
        data2.POSITIONS{lastpos+n}.GRID_ROW=data2.POSITIONS{lastpos+n}.GRID_ROW-1;
        %% change LABEL
        data2.POSITIONS{lastpos+n}.LABEL=[posname, ...
            '_',num2str(data2.POSITIONS{lastpos+n}.GRID_COL,'%.3u'), ...
             '_',num2str(data2.POSITIONS{lastpos+n}.GRID_ROW,'%.3u')];
        %% change XY positions
        %find the device of XYstage
        for m=1:length(data.POSITIONS{i*4-3}.DEVICES)
            if contains(data.POSITIONS{i*4-3}.DEVICES{m}.DEVICE,'XYStage')
                data2.POSITIONS{lastpos+n}.DEVICES{m}.Y=round(midpointy+((data2.POSITIONS{lastpos+n}.GRID_ROW-midpoint(2))*imwidth*(1-overlap/100)*pixelsize));
                data2.POSITIONS{lastpos+n}.DEVICES{m}.X=round(midpointx+(data2.POSITIONS{lastpos+n}.GRID_COL-midpoint(1))*imwidth*(1-overlap/100)*pixelsize);
                Yoffset=((data2.POSITIONS{lastpos+n}.GRID_ROW-midpoint(2))*imwidth*(1-overlap/100)*pixelsize);
                Xoffset=(data2.POSITIONS{lastpos+n}.GRID_COL-midpoint(1))*imwidth*(1-overlap/100)*pixelsize;
            end
        end
        %% change Z positions
        %find the device of XYstage
        for m=1:length(data.POSITIONS{i*4-3}.DEVICES)
            if contains(data.POSITIONS{i*4-3}.DEVICES{m}.DEVICE,'ManualFocus')
               Zoffset=Yoffset*zslopey + Xoffset*zslopex;
               data2.POSITIONS{lastpos+n}.DEVICES{m}.X=round(midpointz+Zoffset);
                
            end
        end
    end
    lastpos=lastpos+tileconfig(1)*tileconfig(2);
end
json=savejson('',data2,['tile',fname]);

                