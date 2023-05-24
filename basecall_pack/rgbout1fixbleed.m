function rgbout1fixbleed(string,strength,m,b,chprofile)
% write 4-channel images into 8-bit rgb images for display. The four
% channels are assigned as Cyan, Yellow, Magenta, and White. m indicates
% median filter size (0 = no median filtering). b indicates ball radius for
% top hat filtering (0 = no filtering). This version also fixes
% bleedthrough using chprofile
if ~exist('m','var')
    m=0;
end
if ~exist('b','var')
    b=0;
end


if m<=3&&m>0
    m=3;
end



fileselection=strcat(string,'*.tif');
files=dir(fullfile(fileselection));
files={files.name};
L=size(files);

if ~exist('chprofile','var')
    chprofile=[];
end

parfor k=1:L(2)
    currentfile=files{k};
    im=imread(currentfile,1);
    for n=2:4
        im(:,:,n)=imread(currentfile,n);
    end
    
    
    
%     imageC=imread(currentfile,1);
%     imageY=imread(currentfile,2);
%     imageM=imread(currentfile,3);
%     imageW=imread(currentfile,4);
%     %imageC=max(double(imageC-mean(mean(imageC))),0)/double((max(max(imageC))-min(min(imageC))))*255;
    %imageY=max(double(imageY-mean(mean(imageY))),0)/double(max(max(imageY))-min(min(imageY)))*255;
    %imageM=max(double(imageM-mean(mean(imageM))),0)/double(max(max(imageM))-min(min(imageM)))*255;
    %imageW=max(double(imageW-mean(mean(imageW))),0)/double(max(max(imageW))-min(min(imageW)))*255;
    %imageC=max(double(imageC),0)/300*255;
    %imageY=max(double(imageY),0)/300*255;
    %imageM=max(double(imageM),0)/300*255;
    %imageW=max(double(imageW),0)/300*255;
    if m>0

        for n=1:4
            im(:,:,n)=medfilt2(im(:,:,n),[m m]);
        end
%         
%         imageC=medfilt2(imageC,[m m]);
%         imageY=medfilt2(imageY,[m m]);
%         imageM=medfilt2(imageM,[m m]);
%         imageW=medfilt2(imageW,[m m]);
    end
    if b>0
        se=offsetstrel('ball',b,b);
        for n=1:4
            im(:,:,n)=imtophat(im(:,:,n),se);
        end
%         
%         
%         imageC=imtophat(imageC,se);
%         imageY=imtophat(imageY,se);
%         imageM=imtophat(imageM,se);
%         imageW=imtophat(imageW,se);
    end
        
        
    if ~isempty(chprofile)
        im=reshape(uint16(double(reshape(im,[],4))/chprofile),size(im,1),size(im,2),4);%subtract camera baseline and correct for bleeding
    end
        
    
    
    
    [p,n]=size(im(:,:,1));
    imageRGB(1:p,1:n,1:3)=0;
    imageRGB(1:p,1:n,1)=double(sum(im(:,:,[2 3 4]),3))*strength;
    imageRGB(1:p,1:n,2)=double(sum(im(:,:,[1 2 4]),3))*strength;
    imageRGB(1:p,1:n,3)=double(sum(im(:,:,[1 3 4]),3))*strength;
    
    
%     imageRGB(1:p,1:n,1)=double(imageY+imageM+imageW)*strength;
%     imageRGB(1:p,1:n,2)=double(imageY+imageC+imageW)*strength;
%     imageRGB(1:p,1:n,3)=double(imageC+imageM+imageW)*strength;
%     
    imageRGB=uint8(imageRGB);
    filename=strcat('RGB',currentfile);
    imwrite(imageRGB, filename);
end
end


    