function rgbout3d(string,strength)
% write 4-channel images into 8-bit rgb images for display. This version
% accepts single-plane z-stack tif files from Volocity.

seqfolders=dir([string,'*']); %sequencing folder names

for k=1:length(seqfolders)
    cd(seqfolders);
    files=dir('*.tif');
    tcz=zeros(length(files),3);
    for i=1:length(files)
        tcz(i,:)=cell2mat(textscan(files(i).name,'T%uC%uZ%u'));
    end
    
    im=cell(1,1,length(files));
    for i=1:length(files)
        im{i}=imread(files(i).name);
    end
    im=cell2mat(im);
    
    
    
    [m,n]=size(imageC);
    imageRGB(1:m,1:n,1:3)=0;
    imageRGB(1:m,1:n,1)=(imageY+imageM+imageW)*strength;
    imageRGB(1:m,1:n,2)=(imageY+imageC+imageW)*strength;
    imageRGB(1:m,1:n,3)=(imageC+imageM+imageW)*strength;
    
    imageRGB=uint8(imageRGB);
    filename=strcat('RGB',currentfile);
    imwrite(imageRGB, filename);
    cd ..
end
end


    