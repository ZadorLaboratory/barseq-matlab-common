function mmmaxproj(name)
% For files acquired from mm, read ome tiff data, make maxproj. this script
%runs in each mda acquisition folder.
%Requires bio-format package for matlab and sort_nat.

files=dir(['*',name,'*.ome.tif']);
folders=dir(['*',name,'*']);

if isempty(files) %each position is stored as single-plane images in folders
    sp=1;
else
    sp=0;
end

mkdir('maxproj');
switch sp
    case 1 %positions are stored as single-plane images
        foldersname={folders.name};
        foldersname=foldersname(cell2mat({folders.isdir})&cellfun(@(x) x(1)~='.'& ~ismember(x,'maxproj','rows'),foldersname));
        filedone=zeros(length(foldersname),1);
        for n=1:length(foldersname)
            %check if the file has been processed.
            if filedone(n)==0
                cd(foldersname{n});
                imfiles=dir('*.tif');
                data=bfopen(imfiles(1).name);
                
                %match series name to file names
                idx=zeros(size(data,1),1);
                for i=1:size(data,1)
                    h=textscan(data{i,1}{1,2},'%s','delimiter',';');
                    h=h{1}{2};
                    [~,idx(i)]=ismember(h,foldersname);
                end
                %for each series matching a filename, make maxproj
                for m=1:length(idx)
                    if idx(m)>0
                        %parse description of files to find image channel and z
                        %information
                        
                        for i=1:size(data{m,1},1)
                            h=textscan(data{m,1}{i,2},'%s','delimiter',';');
                            h=h{1};
                            tf=contains(h,'Z')&contains(h,'=');
                            z=textscan(h{find(tf,1)},'%*s %u %u','delimiter',{'=','/'});
                            zc(i,1)=z{1};
                            tf=contains(h,'C')&contains(h,'=');
                            c=textscan(h{find(tf,1)},'%*s %u %u','delimiter',{'=','/'});
                            zc(i,2)=c{1};
                        end
                        %read channel names from metadata
                        ch=cell(c{2},1);
                        for i=1:c{2}
                            ch{i}=data{m,4}.getChannelName(0,i-1).toCharArray';
                            if i==1
                                chall=ch{i};
                            else
                                chall=[chall,'_',ch{i}];
                            end
                            
                        end
                        
                        %make maxproj
                        maxproj=zeros(size(data{m,1}{1,1},1),size(data{m,1}{1,1},2),c{2});
                        
                        for i=1:c{2}
                            maxproj(:,:,i)=max(cell2mat(reshape(data{m,1}(zc(:,2)==i,1),1,1,z{2})),[],3);
                        end
                        maxproj=uint16(maxproj);
                        %writemaxproj
                        imwrite(maxproj(:,:,1),['../maxproj/',foldersname{idx(m)},'_',chall,'.tif']);
                        if size(maxproj,3)>1
                            for i=2:size(maxproj,3)
                                imwrite(maxproj(:,:,i),['../maxproj/',foldersname{idx(m)},'_',chall,'.tif'],'WriteMode','Append');
                            end
                        end
                        %end
                        
                        filedone(idx(m))=1;
                    end
                end
                cd ..
            end
        end
    case 0 %positions are stored directly as stacks in this folder
        files=sort_nat({files.name});
        filedone=zeros(length(files),1);
        for n=1:length(files)
            %check if the file has been processed.
            if filedone(n)==0
                data=bfopen(files{n});
                %match series name to file names
                idx=zeros(size(data,1),1);
                for i=1:size(data,1)
                    h=textscan(data{i,1}{1,2},'%s','delimiter',';');
                    h=h{1}{2};
                    [~,idx(i)]=ismember([h,'.ome.tif'],files);
                end
                %for each series matching a filename, make maxproj
                for m=1:length(idx)
                    if idx(m)>0
                        %parse description of files to find image channel and z
                        %information
                        for i=1:size(data{m,1},1)
                            h=textscan(data{m,1}{i,2},'%s','delimiter',';');
                            h=h{1};
                            tf=contains(h,'Z')&contains(h,'=');
                            z=textscan(h{find(tf,1)},'%*s %u %u','delimiter',{'=','/'});
                            zc(i,1)=z{1};
                            tf=contains(h,'C')&contains(h,'=');
                            c=textscan(h{find(tf,1)},'%*s %u %u','delimiter',{'=','/'});
                            zc(i,2)=c{1};
                        end
                        %read channel names from metadata
                        ch=cell(c{2},1);
                        for i=1:c{2}
                            ch{i}=data{m,4}.getChannelName(0,i-1).toCharArray';
                            if i==1
                                chall=ch{i};
                            else
                                chall=[chall,'_',ch{i}];
                            end
                            
                        end
                        
                        %make maxproj
                        maxproj=zeros(size(data{m,1}{1,1},1),size(data{m,1}{1,1},2),c{2});
                        
                        for i=1:c{2}
                            maxproj(:,:,i)=max(cell2mat(reshape(data{m,1}(zc(:,2)==i,1),1,1,z{2})),[],3);
                        end
                        %writemaxproj
                        maxproj=uint16(maxproj);
                        imwrite(maxproj(:,:,1),['maxproj/',files{idx(m)},'_',chall,'.tif']);
                        if size(maxproj,3)>1
                            for i=2:size(maxproj,3)
                                imwrite(maxproj(:,:,i),['maxproj/',files{idx(m)},'_',chall,'.tif'],'WriteMode','Append');
                            end
                        end
                        filedone(idx(m))=1;
                    end
                end
                
            end
        end
end








