function mmstitchmax(ch,refch,overlap,singleslice)

%make max projections of stacks given channel numbers, stitch by reading
%positions from filename. Output stitched maxprojection images.singleslice
%controls whether output each channel as a separate image or all channels
%as a stack. Useful for large tiles that 


%% initialize

%ch=3; %number of channels
%overlap=0.15; %fraction of image overlap
%refch=2; %channel used for stitching

if ~exist('singleslice','var')
    singleslice=0;
end


%% Make max projections and flip both x and y, also collects file names.
files=dir('*.tif');
files=sort_nat({files.name});
%mkdir('stitch');
sipxy=cell(length(files),3);




for i=1:length(files)
    sipxy(i,:)=textscan(files{i},'%*s %s %u %u',1,'Delimiter',{'_','.'});
    info=imfinfo(files{i});
end

xy=cell2mat(sipxy(:,2:3));      
pos=cellfun(@cell2mat,sipxy(:,1),'UniformOutput',false);
uniqpos=sort_nat(unique(pos));

%build common filename. this may needs to be standardized in the future.
c=textscan(files{i},'%s','Delimiter','_');
c=c{1}';
c(2,:)={'_'};
c=cell2mat(c(6:end-1));

%% write channel to be used for stitching in stitch folder
mkdir stitchchannel
for n=1:ch
    mkdir(['stitchchannel/ch',num2str(n)]);
    for i=1:length(files)
        im=imrotate(imread(files{i},n),180); %read and invert image.
        imwrite(im,['stitchchannel/ch',num2str(n),'/',files{i}]);
    end
end

%%	Stitch reference channel using ImageJ
addpath('C:\Fiji.app\scripts');
%javaaddpath 'C:\Program Files\MATLAB\R2018b\java\mij.jar'
javaaddpath('C:\Program Files\MATLAB\R2021b\java\mij.jar')
Miji(false);
 cd stitchchannel
mkdir stitchedref
for i=1:length(uniqpos)
    %for each position, count the number of columns and rols
    colnum=max(xy(ismember(pos,uniqpos(i)),1))+1;
    rownum=max(xy(ismember(pos,uniqpos(i)),2))+1;
    %first stitch using ImageJ
    MIJ.run('Grid/Collection stitching', ['type=[Filename defined position] ',...
        'order=[Defined by filename         ] ',...
        'grid_size_x=',num2str(colnum),' grid_size_y=',num2str(rownum),' '....
        'tile_overlap=',num2str(floor(overlap*100)),' ',...
        'first_file_index_x=0 first_file_index_y=0 ',...
        'directory=[',[pwd,'\ch',num2str(refch)],'] ',...
        'file_names=',['MAX_',uniqpos{i},'_{xxx}_{yyy}','.tif '],...
        'output_textfile_name=',[uniqpos{i},'TileConfiguration.txt '],...
        'fusion_method=[Linear Blending] ',...
        'regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 ',...
        'compute_overlap ',...
        'computation_parameters=[Save computation time (but use more RAM)] ',...
        'image_output=[Fuse and display]']);
    MIJ.run('Save',['Tiff..., path=[',pwd,'\stitchedref\','stitchedref_',uniqpos{i},'.tif]']);
    MIJ.run("Close All")
end


%% transform other channels.
for n=1:ch
    mkdir(['ch',num2str(n),'stitched']);
    for i=1:length(uniqpos)
        MIJ.run('Grid/Collection stitching', ['type=[Positions from file] ',...
            'order=[Defined by TileConfiguration] ',...
            'directory=[',[pwd,'\ch',num2str(n)],'] ',...
            'layout_file=',['..\ch',num2str(refch),'\',uniqpos{i},'TileConfiguration.registered.txt '],...
            'fusion_method=[Linear Blending] ',...
            'regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 ',...
            'computation_parameters=[Save computation time (but use more RAM)] ',...
            'image_output=[Fuse and display]']);
        MIJ.run('Save',['Tiff..., path=[',pwd,'\ch',num2str(n),'stitched\',uniqpos{i},'.tif]']);
         MIJ.run("Close All")
    end
end
MIJ.exit; 
%% piece channels together.
for i=1:length(uniqpos)
    im=imread(['ch1stitched/',uniqpos{i},'.tif']);
    if singleslice==0
        imwrite(im,['stitched',uniqpos{i},'.tif']);
        for n=2:ch
            im=imread(['ch',num2str(n),'stitched/',uniqpos{i},'.tif']);
            imwrite(im,['stitched',uniqpos{i},'.tif'],'WriteMode','Append');
        end
    else
        for n=1:ch
            im=imread(['ch',num2str(n),'stitched/',uniqpos{i},'.tif']);
            imwrite(im,['stitched',uniqpos{i},'_ch',num2str(n),'.tif']);
        end
    end
end

cd ..



%%


	





