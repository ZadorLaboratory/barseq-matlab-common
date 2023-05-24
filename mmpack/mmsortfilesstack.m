%sort sequencing images generated by seq-o-matic. Sequence files(folders)
%should be named by sample name followed by 'seq' followed by an index
%(usually automaticaly generated by mm.
folders=dir('*seq*');
folders=sort_nat({folders.name});
for i=1:length(folders)
    rootdir=cd(folders{i});
    poslist=dir('*Pos*.tif');
    poslist=sort_nat({poslist.name});
    
    for n=1:length(poslist)
        info=imfinfo(poslist{n});
        im=zeros(info(1).Height,info(1).Width,4,floor(length(info)/4));
        
        for m=1:floor(length(info)/4)
            for p=1:4
                im(:,:,p,m)=imread(poslist{n},m+(p-1)*7);
            end
        end
        im=max(im,[],4);
        %make max projectios
        for m=1:4
            if m==1
                imwrite(uint16(im(:,:,m)),['../MAX_',folders{i},'_',poslist{n},'_seq',num2str(i),'.tif']);
            else
                imwrite(uint16(im(:,:,m)),['../MAX_',folders{i},'_',poslist{n},'_seq',num2str(i),'.tif'],'WriteMode','Append');
            end
        end
        
    end
    cd ..
end
%% separate by Pos
for i=1:length(poslist)
    poslist{i}=['Pos',num2str(i-1)];
end

for i=1:length(poslist)
    mkdir(poslist{i});
    files=dir(['*',poslist{i},'*.tif']);
    files={files.name};
    for n=1:length(files)
        movefile(files{n},[poslist{i},'/',files{n}]);
    end
end
    




%% make process max proj files to make rgb files
for i=1:length(poslist)
    cd(poslist{i});
    %The next two lines are for when image is taken during cleavage
    %xlightsubtractseqbgn('MAX',0.65)
    %xlightseqalignment1('bgnsub');
    %The next line is for when no image is taken during cleavage
    xlightseqalignment('MAX');
    cd aligned
    rgbout1('aligned',3);
    mkdir('../RGB');
    movefile('RGB*.tif','../RGB/');
    cd ..
    cd ..
    
end


%% redo RGB
for i=1:length(poslist)
    cd(poslist{i});
    cd aligned
    rgbout1('aligned',0.4);
    mkdir('../RGB');
    movefile('RGB*.tif','../RGB/');
    cd ..
    cd ..
    
end

%% make movie
for i=1:length(poslist)
    cd(poslist{i});
    cd RGB;
    file=dir('*.tif');
    file=sort_nat({file.name});
    myVideo = VideoWriter('sequence.avi');
    myVideo.FrameRate=2;
    
    open(myVideo);
    for n=1:length(file)
        writeVideo(myVideo,imread(file{n}));
    end
    close(myVideo);
    cd ..
    cd ..
    
end






%% manually select cells and basecall
for i=1:length(poslist)
    rootdir=cd(poslist{i});
    cd RGB;
    files=dir('RGB*.tif');
    files=sort_nat({files.name});
    I=imread(files{1});
    for i=2:length(files)
        I(:,:,:,i)=imread(files{i});
    end
    implay(I,2);
    im=I(:,:,:,1);
    
    cd ..
    
    
    xall=[];yall=[];
    figure;imshow(im);hold on;title('Press "Enter" when done.');
    while 1
        [x,y]=myginput(1,'crosshair');
        if numel(x)==1
            xall=[xall;x];
            yall=[yall;y];
            scatter(x,y,10,'r','filled');
        else
            break
        end
    end
    close all;
    
    
    cd aligned;
    files=dir('aligned*.tif');
    files=sort_nat({files.name});
    im=conv2(imread(files{1}),eye(3));
    for n=2:4
        im(:,:,n)=conv2(imread(files{1},n),eye(3));
    end
    for i=2:length(files)
        for n=1:4
            im(:,:,n,i)=conv2(imread(files{i},n),eye(3));
        end
    end
    
    sig=zeros(length(xall),1,size(im,3),size(im,4));
    %basecall all selected points
    for i=1:length(xall)
        sig(i,1,:,:)=im(uint16(yall(i)),uint16(xall(i)),:,:);
    end
    sig=permute(squeeze(sig),[1 3 2]);
    [intensity,seq]=max(sig,[],3);
    qual=intensity./sqrt(sum(sig.^2,3));
    seqC=char(seq);
    seqC(seq==1)='G';
    seqC(seq==2)='T';
    seqC(seq==3)='A';
    seqC(seq==4)='C';
    
    save('seq.mat','seq','seqC','qual','intensity','sig');
    cd(rootdir);
end

%% plot quality and intensity over cycles
for i=1:length(poslist)
    rootdir=cd(poslist{i});
    cd aligned;
    load seq.mat
    figure
    subplot(1,2,1);
    errorbar(log2(1+mean(intensity)),((log2(1+mean(intensity))-log2(1+std(intensity)))/sqrt(size(intensity,1))),'LineWidth',2,'Color','k');
    set(gca,'ylim',[0 14]);
    %plot(log2(1+mean(intensity)),'LineWidth',2,'Color','k');
    xlabel('Cycles');
    ylabel('Intensity(bits)');
    title(poslist{i});
    subplot(1,2,2);
    errorbar(mean(qual),std(qual)/sqrt(size(qual,1)),'Linewidth',2,'Color','k');
    set(gca,'ylim',[0.5 1.05],'ytick',0.5:0.1:1);
    xlabel('Cycles');
    ylabel('Quality core');
    cd(rootdir);
    title(poslist{i});
end

%% remake RGB
for i=1:length(poslist)
    cd(poslist{i});
    cd aligned
    rgbout1('aligned',1);
    mkdir('../RGB');
    movefile('RGB*.tif','../RGB/');
    cd ..
    cd ..
    
end

