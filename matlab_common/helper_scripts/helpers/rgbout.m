function rgbout(string)
% write 4-channel images into 8-bit rgb images for display. The four
% channels are assigned as Cyan, Yellow, Magenta, and White.'
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
    %imageC=max(double(imageC-mean(mean(imageC))),0)/double((max(max(imageC))-min(min(imageC))))*255;
    %imageY=max(double(imageY-mean(mean(imageY))),0)/double(max(max(imageY))-min(min(imageY)))*255;
    %imageM=max(double(imageM-mean(mean(imageM))),0)/double(max(max(imageM))-min(min(imageM)))*255;
    %imageW=max(double(imageW-mean(mean(imageW))),0)/double(max(max(imageW))-min(min(imageW)))*255;
    imageC=max(double(imageC),0)/300*255;
    imageY=max(double(imageY),0)/300*255;
    imageM=max(double(imageM),0)/300*255;
    imageW=max(double(imageW),0)/300*255;
    
    [m,n]=size(imageC);
    imageRGB(1:m,1:n,1:3)=0;
    imageRGB(1:m,1:n,1)=(imageY+imageM+imageW)/6;
    imageRGB(1:m,1:n,2)=(imageY+imageC+imageW)/6;
    imageRGB(1:m,1:n,3)=(imageC+imageM+imageW)/6;
    
    imageRGB=uint8(imageRGB);
    filename=strcat('RGB',currentfile);
    imwrite(imageRGB, filename);
end
end


    