function [h,im1rgb]=display4colorimg(im,range)
%display four color sequencing image in the usual color code. im should be
%m x n x 4. range is 2 x 1 with black/white point for display
if ~exist('range','var')
    range=[min(im(:)),max(im(:))];
end


im1=reshape(im,[],4);
im1rgb=double(im1)*[0 1 1;1 1 0;1 0 1;1 1 1];
im1rgb=reshape(im1rgb,size(im,1),size(im,2),3);
im1rgb=min(1,((im1rgb)-range(1))/(range(2)-range(1)));
h=imshow(im1rgb);