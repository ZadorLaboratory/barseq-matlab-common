
function pos=stitch_10x_images(ch_count,reg_ch,overlap,x,trim,rescale_factor)
% stitch images using ImageJ and save tforms. Requires imageJ installation 
    cam_x=3200;
    cam_y=3200;
    fprintf('Stitching whole-slice image.\n');
    niestitchgenessingleplane(ch_count,reg_ch,overlap,x,trim); % This doesn't work reliably with DIC channel, adjust imaging settings.
    
    
    fprintf('Making 10x equivalent stitcjed images.\n')
    % make fake 10x images. 
    mkdir 10x
    cd(fullfile(['stitchchannel',num2str(x)],'ch1stitched'));
    f=dir('*.tif');
    f=sort_nat({f.name});
    for i=1:length(f)
        im=imread(f{i});
        for n=2:5
            im(:,:,n)=imread(['../ch',num2str(n),'stitched/',f{i}]);
        end
        %im=imresize(imrotate(im,180),0.5);
        im=imresize(im,rescale_factor);
        imwrite(im(:,:,1),['../../10x/',f{i}]);
        for n=2:5
            imwrite(im(:,:,n),['../../10x/',f{i}],'WriteMode','Append');
        end
    end
    
    cd ../..
    fprintf('Parse tileconfig for each position.\n');
    % parse tileconfig for each pos, sensitive to file name
%    [folders,pos,xxx,yyy]=get_folders();
            

%     %  Make initial transofrmation
%     itform={};
%     tileconfig={};
%     uniqpos=sort_nat(unique(pos));
%     for m=1:numel(uniqpos)
%         tileconfig{m}=[max(xxx(ismember(pos,uniqpos(m)))),max(yyy(ismember(pos,uniqpos(m))))];
%         itform{m}={};
%         for i=1:tileconfig{m}(2) %y
%             for n=1:tileconfig{m}(1) %x
%                 %itform{m}{i,n}=[1 0 0; 0 1 0; (tileconfig{m}(1)-n)*2048/2*0.85 (tileconfig{m}(2)-i)*2048/2*0.85 1];%these are theoretical values
%                 itform{m}{i,n}=[1 0 0; 0 1 0; (tileconfig{m}(1)-n)*(cam_x/2*(1-overlap)) (i-1)*(cam_y/2*(1-overlap)) 1];%these are empirical measures
%             end
%         end
%     end
%     
%     fprintf('Register 20x sequencing images to 10x sequencing images.\n')
%     % Register 40x sequencing images to 10x sequencing images.
%     tic
%     cd 10x
%     allfiles10x=dir('*.tif');
%     allfiles10x=sort_nat({allfiles10x.name});
%     cd ..
%     
%     tform40xto10x={};
%     parfor(i=1:length(folders),4) %not enough memory when slices are large
%         [~,posidx]=ismember(pos{i},uniqpos);
%         file10x=['../../10x/',allfiles10x{posidx}];
%         cd(folders{i})
%         cd original
%         file40x=dir('*seq*.tif');
%         file40x=sort_nat({file40x.name});
%         file40x=file40x{1};
%         itform1=itform{posidx}{yyy(i),xxx(i)};segim1=[];scalediff=rescale_factor;
%         [tform40xto10x{i},~]=mmregister(file40x,file10x,itform1,segim1,scalediff);
%         %copyfile('40xto10xreg.tif',['../../40xto10x/',folders{i},'.tif']);
%         cd ../..
%         
%     end 
%     toc

    %Read 40xto10x tform directly from the ImageJ stitching outputs.
    %reg_ch=5; %for testing
    %x=1;% for testing
    [folders,pos,xxx,yyy]=get_folders();

    fname_all=dir(fullfile(['stitchchannel',num2str(x)],['ch',num2str(reg_ch)],'Pos*.registered.txt'));
    fname_all=sort_nat({fname_all.name})';
    fname_all=cellfun(@(y) fullfile(['stitchchannel',num2str(x)],['ch',num2str(reg_ch)],y),fname_all,'UniformOutput',0);
    
    
    
    %stitchedfname=dir(fullfile('10x','*.tif'));
    %stitchedfname=sort_nat({stitchedfname.name})';
    %stitchedfname=cellfun(@(y) fullfile('10x',y),stitchedfname,'UniformOutput',0);
    %folders=get_folders();
    %
    tforms_converted={};
    
    for i=1:numel(fname_all)
        fname=fname_all{i};
        fid=fopen(fname);
        c=textscan(fid,'%s %*s %s','Delimiter',';');
        fclose(fid);
        tifnames=c{1};
        tforms=c{2};
        I=contains(tifnames,'tif');
        tifnames=tifnames(I);
        tifnames=cellfun(@(x) x(1:end-4),tifnames,'UniformOutput',0);%truncate filenames to get folder names
        tforms=tforms(I);
        tform_xy=cellfun(@(x)textscan(x,'%*f %f %f %*f','Delimiter',',()'),tforms,'UniformOutput',0);
        tform_xy=cell2mat(cellfun(@cell2mat,tform_xy,'UniformOutput',0))/2;
    
        
        %stitchedfileinfo=imfinfo(stitchedfname{i});
        %stitched_width=stitchedfileinfo(1).Width;
        tform_xy(:,1)=max(tform_xy(:,1))-tform_xy(:,1);
        tform_xy(:,1)=tform_xy(:,1)-cam_x*trim/2; %add back the trimed portion
        tform_xy(:,2)=tform_xy(:,2)-cam_y*trim/2; %add back the trimed portion
        
        [~,I1]=ismember(tifnames,folders); %match up files
    
    
        for n=1:numel(I1)
            T=[0.5,0,0;0,0.5,0;tform_xy(n,:),1];
            tforms_converted{I1(n)}=affine2d(T);
        end
    
    end
    tform40xto10x=tforms_converted;

    save('40xto10x.mat','tform40xto10x');
    fprintf('All done. Saved transformations to 40xto10x.mat.\n')
end