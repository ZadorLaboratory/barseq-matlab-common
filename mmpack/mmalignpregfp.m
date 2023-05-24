function mmalignpregfp(slicenames,gfpidx,seqidx)
%align pre-gfp image to cycle 1 image.

for i=1:length(slicenames)
    %read pregfp
    cd pregfp
    info=imfinfo([slicenames{i},'pregfp.tif']);
    gfpim=zeros(info(1).Height, info(1).Width,size(info,1));
    for n=1:size(info,1)
        gfpim(:,:,n)=imread([slicenames{i},'pregfp.tif'],n);
    end
    % read seq01
    cd(['../',slicenames{i},'/original']);
    files=dir([slicenames{i},'*.tif']);
    files=sort_nat({files.name});
    info=imfinfo(files{1});
    seqim=zeros(info(1).Height, info(1).Width,size(info,1));
    for n=1:4
        seqim(:,:,n)=imread(files{1},n);
    end
    %register pregfp to seq01
    %[optimizer,metric] = imregconfig('multimodal');
    %tform = imregtform(gfpim(100:end-100,100:end-100,gfpidx), sum(seqim(100:end-100,100:end-100,seqidx),3), 'rigid', optimizer, metric);
    
    %manually align pregfp to seq01
    [mpoints,fpoints]=cpselect(min((gfpim(:,:,gfpidx)/max(reshape(gfpim(:,:,gfpidx),[],1)))*2,1),min((seqim(:,:,seqidx)/max(reshape(seqim(:,:,seqidx),[],1)))*5,1),'Wait',true);
    tform=fitgeotrans(mpoints,fpoints,'nonreflectivesimilarity');
    
    Rfixed=imref2d(size(sum(seqim,3)));
    alignedgfpim=uint16(imwarp(gfpim,tform,'OutputView',Rfixed));
    
    %figure;imshowpair(min((alignedgfpim(:,:,gfpidx)/max(reshape(gfpim(:,:,gfpidx),[],1)))*1,1),min((seqim(:,:,seqidx)/max(reshape(seqim(:,:,seqidx),[],1)))*3,1));
    imwrite(alignedgfpim(:,:,1),['../../pregfp/aligned',slicenames{i},'pregfp.tif']);
    if size(gfpim,3)>1
        for n=2:size(alignedgfpim,3)
            imwrite(alignedgfpim(:,:,n),['../../pregfp/aligned',slicenames{i},'pregfp.tif'],'WriteMode','Append');
        end
    end
    cd('../../');
end