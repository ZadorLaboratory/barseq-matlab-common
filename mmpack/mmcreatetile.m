function data2=mmcreatetile(fname,tileconfig, imwidth,pixelsize,overlap)
%create tiles from saved pos files.
data=loadjson(fname);

data2=data;
%FOR EACH POSITIONS, create additional tiles specified by tileconfig
midpoint=tileconfig/2-0.5;

for i=1:length(data.POSITIONS)
    data2.POSITIONS((i-1)*tileconfig(1)*tileconfig(2)+(1:tileconfig(1)*tileconfig(2)))=data.POSITIONS(i);
    for n=1:tileconfig(1)*tileconfig(2)
        %create data.Positions{}.properties
        data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.PROPERTIES.OverlapUmX=num2str(imwidth*overlap*pixelsize/100);
        data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.PROPERTIES.OverlapUmY=num2str(imwidth*overlap*pixelsize/100);
        data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.PROPERTIES.OverlapPixelsX=num2str(imwidth*overlap/100);
        data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.PROPERTIES.OverlapPixelsY=num2str(imwidth*overlap/100);
        %add grid positions
        [data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.GRID_COL,data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.GRID_ROW]=ind2sub(tileconfig,n);
        data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.GRID_COL=data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.GRID_COL-1;
        data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.GRID_ROW=data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.GRID_ROW-1;
        %change LABEL
        data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.LABEL=[data.POSITIONS{i}.LABEL, ...
            '_',num2str(data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.GRID_COL,'%.3u'), ...
             '_',num2str(data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.GRID_ROW,'%.3u')];
        %change XY positions
        %find the device of XYstage
        for m=1:length(data.POSITIONS{i}.DEVICES)
            if contains(data.POSITIONS{i}.DEVICES{m}.DEVICE,'XYStage')
                data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.DEVICES{m}.Y=data.POSITIONS{i}.DEVICES{m}.Y+((data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.GRID_ROW-midpoint(2))*imwidth*(1-overlap/100)*pixelsize);
                data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.DEVICES{m}.X=data.POSITIONS{i}.DEVICES{m}.X+(data2.POSITIONS{(i-1)*tileconfig(1)*tileconfig(2)+n}.GRID_COL-midpoint(1))*imwidth*(1-overlap/100)*pixelsize;
                
            end
        end
    end
end
json=savejson('',data2,['tile',fname]);

                