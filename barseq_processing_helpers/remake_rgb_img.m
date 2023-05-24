
function remake_rgb_img(fname,rgb_intensity)
    folders=dir('MAX*');
    folders=sort_nat({folders.name});
    %worker pools
    p=gcp('nocreate');
    if isempty(p)
        parpool("local",24);
    elseif p.NumWorkers<12
        parpool("local",24);
    end
    
    % remake rgb images
    parfor i=1:length(folders)
        %for i=13
        cd(folders{i});
        cd aligned
        rgbout1(['alignedfixed',fname],rgb_intensity);
        mkdir('../RGB');
        movefile('RGB*.tif','../RGB/');
        cd ../..
    end
end
