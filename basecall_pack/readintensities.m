function sig=readintensities(r,im,convmat)
%given a list of rolony coordinates, read the intensities of four channels
%in an area given by convmat. rol is Nx2 coordinates of rolonies. im is
%QxNx4 image. convmat is a convolution matrix.
imc=zeros(size(im));
for n=1:4
    imc(:,:,n)=conv2(im(:,:,n),convmat,'same');
end

sig=zeros(size(r,1),4);
for i=1:size(r,1)
    for n=1:4
        sig(i,n)=imc(round(r(i,1)),round(r(i,2)),n);
    end
end
end
