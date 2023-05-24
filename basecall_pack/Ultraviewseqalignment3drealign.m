function Ultraviewseqalignment3drealign(filename)
%Redo stack alignment in fixed images based on Ultraviewseqalignemnt3d. 

%This version uses FFT phase correlation for image alignment and
%cuts off the border for alignment. This version takes 3d stacks as inputs 
%but only aligns in the xy plane and uses. The z-alignment is done later
%locally during basecalling. This version also uses sparse nmf to find 
%channel bleedthrough.


%fix and preprocess tif files:
files=dir(['fixed*',filename,'*']);
for i=1:length(files)
    fixedfolders{i}=files(i).name;
end


%copy the first cycle fixed images directly as aligned images
%copyfile(fixedfolders{1},['aligned',fixedfolders{1}]);
%cd(['aligned',fixedfolders{1}]);
%files=dir('*.tif');
%for i=1:length(files)
%    movefile(files(i).name,['aligned',files(i).name]);
%end
%cd ..
%chprofiles=cell(1);
%fixedfolders=dir('fixed*');
%fixedfolders={fixedfolders.name};

%align subsequent cycles to cycle 1

parfor(k=2:length(fixedfolders),2);
    alignseq(fixedfolders{k},fixedfolders{1});
end


%movefile('aligned*.tif','aligned/');
%delete('fixed*tif');

end


function alignseq(imagefolder,templatefolder)
%read template sum across seq channels
cd(templatefolder);
filestemplate=dir('*.tif');
tcztemplate=zeros(length(filestemplate),3);
for i=1:length(filestemplate)
    tcztemplate(i,:)=cell2mat(textscan(filestemplate(i).name,'fixedT%uC%uZ%u'));
end

filestemplate(tcztemplate(:,2)>4)=[];%remove non-sequencing channels
tcztemplate(tcztemplate(:,2)>4,:)=[];%remove non-sequencing channels

imtemplate=cell(1,1,length(filestemplate));
for i=1:length(filestemplate)
    imtemplate{i}=imread(filestemplate(i).name);
end
imtemplate=cell2mat(imtemplate);
imsumtemplate=zeros(size(imtemplate,1),size(imtemplate,2),4);
for i=1:4
    imsumtemplate(:,:,i)=max(imtemplate(:,:,tcztemplate(:,2)==i),[],3); %max proj
end
imsumtemplate=sum(imsumtemplate,3);%z-proj of channel-summed template
clearvars imtemplate filestemplate tcztemplate %the non-summed template images are not used after this point

%read moving image and sum
cd(['../',imagefolder]);
files=dir('*.tif');
tcz=zeros(length(files),3);
for i=1:length(files)
    tcz(i,:)=cell2mat(textscan(files(i).name,'fixedT%uC%uZ%u'));
end

im=cell(1,1,length(files));
for i=1:length(files)
    im{i}=imread(files(i).name);
end
im=cell2mat(im);

imsum=zeros(size(im,1),size(im,2),4);
for i=1:4
    imsum(:,:,i)=max(im(:,:,tcz(:,2)==i),[],3); 
end
imsum=sum(imsum,3); %z-proj of channel-summed moving

%align moving sum to template sum
%tform = imregcorr(imsum(500:end-500,500:end-500),imsumtemplate(500:end-500,500:end-500),'rigid','Window',0);
%Rfixed=imref2d(size(imsumtemplate));
 
%align with ecc
par.transform = 'euclidean'; 
par.levels = 7;
par.iterations = 100;
ransacWarp=iat_ecc(imsum(500:end-500,500:end-500),imsumtemplate(500:end-500,500:end-500),par);
        
        [M,N]=size(imsumtemplate);
        %write aligned image to file.

%write aligned image to file.
mkdir(['../aligned',imagefolder])
cd(['../aligned',imagefolder]);
for i=1:length(files)
    %imwrite(uint16(imwarp(im(:,:,i),tform,'OutputView',Rfixed)),['aligned',files(i).name]);
    %for ecc alignment
    imwrite(uint16(iat_inverse_warping(im(:,:,i),ransacWarp,par.transform,1:N, 1:M)),['aligned',files(i).name]);
    
end


cd ..
end







