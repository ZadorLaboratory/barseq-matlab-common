function mm2maxprojimaging(fname,chseq,medorder,pausetime,remotelogin,remotepath)
%make max proj of single-plane tiff files in each folder generated by micromanager. The
%folders contain "name", and chseq indicate the sequence of channels in the
%final file. This script should be run in the folder for each
%acquisition. Unlike mmmaxproj, mmmaxprojfast does not require the
%bio-format package and processes files faster, but requires manually
%determining the sequence of channels.
%This version can do median filtering on the z-axis to reduce noise without
%affecting x-y resolution.
%remote path should be username@domain:path. make sure that path exists
% This version works on micro-manager2 produced files.
%%
workernum=[2 8];
if ~exist('medorder','var')
    medorder=0;
end


if ~exist('pausetime','var')
    pausetime=30;
end
passedtime=0;

if ~exist('remotepath','var')
    remotepath=[];
else
    if remotepath(end)=='/'
        remotepath=remotepath(1:end-1);
    end
end

if ~exist('remotelogin','var')
    remotelogin=[];
end

mkdir maxproj
finishedfolders={};
nn=1;
%%
%start parpool
delete(gcp('nocreate'));
parpool(workernum);%limit the number of workers to limit IO.
%%
%make remote folder

if ~isempty(remotelogin)&&~isempty(remotepath)
    eval(['!ssh ',remotelogin,' mkdir -p ',remotepath])
    fprintf('Created remote folder.\n');
end
%%
firsttime=1;
while nn<5
    %%
    % wait then look for folders
    fprintf('Waiting %u sec.\n',max(pausetime-passedtime,1));
    pause(max(pausetime-passedtime,1));
    tic
    folders=dir(['*',fname,'*']);
    folders=sort_nat({folders.name});
    % find the finished folders by counting files
    newfolders=folders(~ismember(folders,finishedfolders));
    filecounts=cellfun(@(x) numel(dir([x,'/*.tif'])),newfolders);
    if firsttime %if first time, count the largest filecounts
        maxfilecounts=max(filecounts);
        firsttime=0;
    end
    %all newfolders with the max file counts are completed.
    
    if isempty(filecounts)%if no new folder, increase nn. Otherwise nn==1;
        nn=nn+1;
        newfolders={};
        fprintf('No new folder completed,counter = %u.\n',nn);
    else
        nn=1;
        newfolders=newfolders(filecounts==maxfilecounts);
    end
    %%
    if nn==5
        newfolders=folders(~ismember(folders,finishedfolders));
        %if there are still folders that doesn't have max file counts after
        %5 tries, then try to maxproj anyway before quiting the loop. Throw
        %a warning message. This should not happen under normal imaging
        %conditions.
        warning('Incomplete folders after 5 tries, maxproj anyway.\n');
    end
    %%
    
    newfiles={};
    parfor n=1:length(newfolders)
        cd(newfolders{n});
        files=dir('*.tif');
        files={files.name};
        info=imfinfo(files{1});
        %parse file names
        ctz=cell(length(files),3);
        for m=1:length(files)
            ctz(m,:)=textscan(files{m},'%*s %s %*s %s %u','Delimiter',{'_','.','_z'});
            %find sequencing channels

        end
        ctz(:,1)=cellfun(@cell2mat,ctz(:,1),'UniformOutput',false);

        ch=unique(ctz(:,1));
        %sort channels using the sequence given.
        ch=ch(chseq);


        %make max projectios
        for m=1:length(ch)
            idx=find(strcmp(ctz(:,1),ch(m)));
            im=zeros(info.Height,info.Width,length(idx));
            for q=1:length(idx)
                im(:,:,q)=imread(files{idx(q)});
            end
            if medorder==0||size(im,3)==1
                if m==1
                    imwrite(uint16(max(im,[],3)),['../maxproj/MAX_',newfolders{n},'.tif']);
                else
                    imwrite(uint16(max(im,[],3)),['../maxproj/MAX_',newfolders{n},'.tif'],'WriteMode','Append');
                end
            else
                if m==1
                    imwrite(uint16(max(medfilt1(im,medorder,[],3),[],3)),['../maxproj/MAX_',newfolders{n},'.tif']);
                else
                    imwrite(uint16(max(medfilt1(im,medorder,[],3),[],3)),['../maxproj/MAX_',newfolders{n},'.tif'],'WriteMode','Append');
                end
            end

        end
        newfiles{n}=['maxproj/MAX_',newfolders{n},'.tif'];

        cd ..

    end
    finishedfolders=[finishedfolders,newfolders];
    
    passedtime=toc;
    for n=1:length(newfolders)
        if ~isempty(remotelogin)&&~isempty(remotelogin)
            eval(['!start /B scp ',newfiles{n},' ',remotelogin,':',remotepath,'/']);
        end
    end

end


% if ~isempty(remotelogin)&&~isempty(remotelogin)
%     eval(['!scp -r maxproj/ ',remotelogin,':',remotepath,'/']);
% end
end
