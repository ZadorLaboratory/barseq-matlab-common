function cropconcat(positions);
files=dir(fullfile('RGB*.tif'));
files=sort({files.name});
L=size(files);
ROIcount=size(positions,1);

for n=1:ROIcount
image=imread(files{1});
concatimage=imcrop(image,positions(n,:));
for k=2:L(2)
    currentfile=files{k};
    image=imread(currentfile);
    croppedimage=imcrop(image, positions(n,:));
    concatimage=cat(2,concatimage,croppedimage);
end

outputfile=strcat('ROI',num2str(n),'.tif');
imwrite(concatimage,outputfile);
end
end

