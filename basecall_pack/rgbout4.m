function rgbout4(string,strength)
% write 4-channel images into 8-bit rgb images for display. The four
% channels are assigned as Cyan, Yellow, Magenta, and White.'This version
% works for non-bgn subtracted files.
fileselection=strcat(string,'*.tif');
files=dir(fullfile(fileselection));
files={files.name};
L=size(files);

for k=1:L(2)
    currentfile=files{k};
    imageC=imread(currentfile,1);
    imageY=imread(currentfile,2);
    imageM=imread(currentfile,3);
    imageW=imread(currentfile,4);
    
    imageC=imageC-imopen(imageC,strel('ball', 50, 50));
    imageY=imageY-imopen(imageY,strel('ball', 50, 50));
    imageM=imageM-imopen(imageM,strel('ball', 50, 50));
    imageW=imageW-imopen(imageW,strel('ball', 50, 50));
    
    
    [m,n]=size(imageC);
    imageRGB(1:m,1:n,1:3)=0;
    imageRGB(1:m,1:n,1)=(imageY+imageM+imageW)*strength;
    imageRGB(1:m,1:n,2)=(imageY+imageC+imageW)*strength;
    imageRGB(1:m,1:n,3)=(imageC+imageM+imageW)*strength;
    
    imageRGB=uint8(imageRGB);
    filename=strcat('RGB',currentfile);
    imwrite(imageRGB, filename);
end
end


    