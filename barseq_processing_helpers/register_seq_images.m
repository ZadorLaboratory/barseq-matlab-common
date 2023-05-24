
function register_seq_images(fname,chprofile20x,ball_radius,chshift20x,rgb_intensity,local)
    if ~exist('local','var')
        local=1;
    end
    folders=dir('MAX*');
    folders=sort_nat({folders.name});
    %worker pools
    p=gcp('nocreate');
    if isempty(p)
        parpool("local",32);
    elseif p.NumWorkers<12
        parpool("local",32);
    end
    p=gcp('nocreate');
    fprintf('Starting image registration on %u workers...\n',p.NumWorkers)
    tic
    parfor i=1:length(folders)
    %for i=1
    %    for i=13
        cd(folders{i});
        if isfolder('original')
            cd original
            try
                movefile([fname,'*.tif'],'../')
            catch ME
                %rethrow(ME)
            end
            cd ..
        end
        if local==1
            [~,~,warnmsg]=mmseqalignmentsinglethread_local(fname,chprofile20x,ball_radius,chshift20x);
        else
            [~,~,warnmsg]=mmseqalignmentsinglethread(fname,chprofile20x,ball_radius,chshift20x);
        end

        if ~isempty(warnmsg)
            warning('%s: %s\n',folders{i},warnmsg);
        end
        cd aligned
        rgbout1(['alignedfixed',fname],rgb_intensity);
        mkdir('../RGB');
        movefile('RGB*.tif','../RGB/');
        cd ../..
    end
    regtime=toc;
    fprintf('Registration finished, total elapsed time is %u hours %u minutes %u seconds',round(regtime/3600),round(rem(regtime,3600)/60),round(rem(regtime,60)))
end