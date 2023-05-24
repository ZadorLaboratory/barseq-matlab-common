function tform1=mmalignslice(files,refch,bgn,tform1)

if ~exist('tform1','var')
    
    
    %align each slice to the slice after it
    [optimizer,metric]=imregconfig('multimodal');
    optimizer.InitialRadius=optimizer.InitialRadius/5;
    optimizer.MaximumIterations=300;
    tform=cell(length(files)-1,1);
    
    %find tform of adjacent images
    parfor i=1:length(files)-1
        im=imread(files{i},refch)-bgn;
        imtemplate=imread(files{i+1},refch)-bgn;
        
        tform{i}=imregtform(im(100:end-100,100:end-100),imtemplate(100:end-100,100:end-100),'rigid',optimizer, metric);
        %register image using intensity based cross correlation allowing only
        %translation and rotation.
    end
    %calculate accumulated tform
    tform1{length(files)}=affine2d;
    for i=length(files)-1:-1:1
        tform1{i}=affine2d(tform{i}.T*tform1{i+1}.T);
    end
end


%transform images from the back
im=imread(files{end},1);
 Rfixed=imref2d(size(im));
 %im(:,:,n)=imwarp(imread(imagename,n),tform,'OutputView',Rfixed);
   
for i=1:length(files)
    info=imfinfo(files{i});
    for n=1:length(info)
        im=imwarp(imread(files{i},n),affine2d,'OutputView',Rfixed);
        
        im1=imwarp(im,tform1{i},'OutputView',Rfixed);
        
        if i==1
            imwrite(im1,['regch',num2str(n),'.tif']);
            imwrite(im,['ch',num2str(n),'.tif']);
            
        else
            imwrite(im1,['regch',num2str(n),'.tif'],'WriteMode','Append');
            imwrite(im,['ch',num2str(n),'.tif'],'WriteMode','Append');
        end
    end
end
    
    
    
    
    
        