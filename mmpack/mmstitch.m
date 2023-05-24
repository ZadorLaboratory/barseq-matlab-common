function mmstitch(ch,refch,overlap)

%make max projections of stacks given channel numbers, stitch by reading positions from filename. Output stitched maxprojection images.


%% initialize

%ch=3; %number of channels
%overlap=0.15; %fraction of image overlap
%refch=2; %channel used for stitching

%% Make max projections and flip both x and y, also collects file names.
files=dir('*.tif');
files=sort_nat({files.name});
mkdir('max');
sipxy=cell(length(files),5);
parfor i=1:length(files)
    sipxy(i,:)=textscan(files{i},'%s %s %*s %s %u %u %*s %*s','Delimiter',{'_','.'});
    info=imfinfo(files{i});
    im=zeros(info(1).Height,info(1).Width,ch,length(info)/ch);
    for m=1:ch
        for n=1:floor(length(info)/ch)
            im(:,:,m,n)=imread(files{i},(m-1)*(length(info)/ch)+n);
        end
    end
    %max projection
    im=max(im,[],4);
    % fix camera rotation
    im=imrotate(im,180);
    %write files
    for n=1:ch
        imwrite(uint16(im(:,:,n)),['max/MAX_ch',num2str(n),'_',files{i}]);
    end
end

sipxy(:,1:3)=cellfun(@cell2mat,sipxy(:,1:3),'UniformOutput',false);
        


%%	Stitch reference channel
%find all samples
a=pwd;
files=dir('*.tif');
files=sort_nat({files.name});
sipxy=cell(length(files),5);
for i=1:length(files)
    sipxy(i,:)=textscan(files{i},'%s %s %*s %s %u %u %*s %*s','Delimiter',{'_','.'});
end
sipxy(:,1:3)=cellfun(@cell2mat,sipxy(:,1:3),'UniformOutput',false);
        

addpath('C:\Fiji.app\scripts');
javaaddpath 'C:\Program Files\MATLAB\R2018b\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2021b\java\mij.jar'
%javaaddpath 'C:\Program Files\MATLAB\R2018b\java'
Miji
cd(a)
cd max
uniqpos=sort_nat(unique(sipxy(:,3)));
mkdir stitch
for i=1:length(uniqpos)
	%for each position, count the number of columns and rols
    colnum=max(cell2mat(sipxy((ismember(sipxy(:,3),uniqpos(i))),4)))+1;
    rownum=max(cell2mat(sipxy((ismember(sipxy(:,3),uniqpos(i))),4)))+1;
    %first stitch using ImageJ
    MIJ.run('Grid/Collection stitching', ['type=[Filename defined position] ',...
        'order=[Defined by filename         ] ',...
        'grid_size_x=',num2str(colnum),' grid_size_y=',num2str(rownum),' '....
        'tile_overlap=',num2str(floor(overlap*100)),' ',...
        'first_file_index_x=0 first_file_index_y=0 ',...
        'directory=[',pwd,'] ',...
        'file_names=',['MAX_ch',num2str(refch),'_',sipxy{1,1},'_',sipxy{1,2},'_MMStack_',uniqpos{i},'_{xxx}_{yyy}','.ome.tif '],...
        'output_textfile_name=',['stitch\',uniqpos{i},'TileConfiguration.txt '],...
        'fusion_method=[Linear Blending] ',...
        'regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 ',...
        'compute_overlap ',...
        'computation_parameters=[Save computation time (but use more RAM)] ',...
        'image_output=[Fuse and display]']);
    MIJ.run('Save',['Tiff..., path=[',pwd,'\stitch\','MAX_stitch_ch',num2str(refch),'_',uniqpos{i},'.tif]']);
    MIJ.run("Close All")
end
cd ..
 MIJ.exit; 
	


%%





