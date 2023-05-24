function mmalignBCtogenemanual(bcidx,geneidx)
%align barcode sequencing images

bcfiles=dir('*BC*.tif');
genefiles=dir('*gene*.tif');
bcfiles=sort_nat({bcfiles.name});
genefiles=sort_nat({genefiles.name});

%read first cycle images
geneim=imread(genefiles{1},geneidx);
bcim=imread(bcfiles{1},bcidx);

%register images and find transformation
%[optimizer,metric] = imregconfig('multimodal');
%tform = imregtform(bcim(100:end-100,100:end-100), geneim(100:end-100,100:end-100), 'rigid', optimizer, metric);

%%manually align pregfp to seq01
[mpoints,fpoints]=cpselect(double(bcim(100:end-100,100:end-100))/double(max(bcim(:))),double(geneim(100:end-100,100:end-100))/double(max(geneim(:))),'Wait',true);
tform=fitgeotrans(mpoints,fpoints,'nonreflectivesimilarity');


Rfixed=imref2d(size(sum(geneim,3)));
alignedbcim=uint16(imwarp(bcim,tform,'OutputView',Rfixed));
imwrite(geneim,'comp.tif');
imwrite(alignedbcim,'comp.tif','WriteMode','Append');

%transform all barcode files
for i=1:length(bcfiles)
    bcinfo=imfinfo(bcfiles{i});
    im=zeros(bcinfo(1).Height, bcinfo(1).Width,size(bcinfo,1));
    for n=1:size(bcinfo,1)
        im(:,:,n)=imread(bcfiles{i},n);
    end
    alignedim=uint16(imwarp(im,tform,'OutputView',Rfixed));
    %figure;imshowpair(min((alignedgfpim(:,:,gfpidx)/max(reshape(gfpim(:,:,gfpidx),[],1)))*1,1),min((seqim(:,:,seqidx)/max(reshape(seqim(:,:,seqidx),[],1)))*3,1));
    imwrite(alignedim(:,:,1),['reg',bcfiles{i}]);
    for n=2:size(alignedim,3)
        imwrite(alignedim(:,:,n),['reg',bcfiles{i}],'WriteMode','Append');
    end
end