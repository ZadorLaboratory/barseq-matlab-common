function [sig,logsig]=mmimgcallpoints(slicenames,ftlroi)
%given list of coordinates tflroi, read signals in pregfp images.
sig=cell(length(slicenames),1);
lroi=ftlroi;
cd pregfp
seqfiles=dir('aligned*.tif');
seqfiles=sort_nat({seqfiles.name});
radius=50;%bgn subtraction radius.
for i=1:length(slicenames)
    info=imfinfo(seqfiles{i});
    lim=zeros(info(1).Height,info(1).Width,length(info));
    for n=1:length(info)
        lim(:,:,n)=imread(seqfiles{i},n);
        %median filter and bgn subtraction
        lim(:,:,n)=medfilt2(lim(:,:,n));
        lim(:,:,n)=max(lim(:,:,n)-imopen(lim(:,:,n),strel('ball',radius,radius)),0);
    end
    %also calculate the log of images
    loglim=log2(lim+1);
    %readout signals
    currsig=ones(length(lroi{i}),length(info));
    %convolute by 7x7
    for n=1:size(lim,3)
        lim(:,:,n)=conv2(lim(:,:,n),ones(7),'same');
        loglim(:,:,n)=conv2(loglim(:,:,n),ones(7),'same');
    end
    
    for n=1:length(lroi{i}) %for each rolony, readout signal
        currsig(n,:)=reshape(lim(lroi{i}(n,2),lroi{i}(n,1),:),1,size(lim,3));
        currlogsig(n,:)=reshape(loglim(lroi{i}(n,2),lroi{i}(n,1),:),1,size(loglim,3));
    end
    
    sig{i}=currsig;
    logsig{i}=currlogsig;
end
save('pregfpsig.mat','sig','logsig');
cd ..
end