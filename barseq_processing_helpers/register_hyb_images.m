
function register_hyb_images(ch,radius,nuclearch,reg_cycle,chradius,chprofilehyb,chshifthyb, rgbintensity)
    %register hyb to first seq cycle, copy nuclear images
%%  
    clc
    
    p=gcp('nocreate');
    if isempty(p)
        parpool("local",24);
    elseif p.NumWorkers<12
        parpool("local",24);
    end
    
    
    folders=dir('MAX*');
    folders=sort_nat({folders.name});
    fprintf(['Registering hyb to seq cycle ', num2str(reg_cycle),' on ', num2str(p.NumWorkers) ' workers...\n'])
    parfor i=1:length(folders)
    %for i=37
        cd(folders{i});
        mmalignhybtoseq('nuclear','geneseq',ch,chprofilehyb,radius,nuclearch,reg_cycle,chradius,chshifthyb);
        rgbout1('aligned',rgbintensity);
        movefile('aligned*.tif','aligned/');
        movefile('RGB*.tif','RGB/');
        cd ..
    
    end
    fprintf('Registration finished, copying files to make a hyb copy ...\n')
    % Make sequential rounds images
    for i=1:length(folders)
        copyfile([folders{i},'/aligned/alignednuclear.tif'],[folders{i},'/aligned/alignedhyb01.tif']);
    end
    fprintf('All done!\n')
end
