function mmfishcrop(hroi,slicenames,fishidx,score,scorethresh,rfp,rfpthresh)
%crop and concatenate fish images with 1st cycle sequencing image. This
%version produces RGB sequencing images.

cd('pregfp');

fishfiles=dir(fullfile('aligned*.tif'));
[fishfiles,~]=sort_nat({fishfiles.name});
radius=100;
ball=strel('ball', radius, radius);
cd ..

parfor i=1:length(slicenames)
    
    im=imread(['pregfp/',fishfiles{i}],fishidx(1));
    im=im-imopen(im, ball);
    if length(fishidx)>1
        for n=2:length(fishidx)
            im(:,:,n)=imread(['pregfp/',fishfiles{i}],fishidx(n));
            im(:,:,n)=im(:,:,n)-imopen(im(:,:,n), ball);
        end
    end
    
    %im=im*2;
    
    cd([slicenames{i},'/aligned']);
    seqfiles=dir('*.tif');
    seqfiles=sort_nat({seqfiles.name});
    im(:,:,length(fishidx)+1)=imread(seqfiles{1},1);
    for n=1:4
        im(:,:,length(fishidx)+n)=imread(seqfiles{1},n)*2;
        im(:,:,length(fishidx)+n)=im(:,:,length(fishidx)+n)-imopen(im(:,:,length(fishidx)+n),ball);
    end
    im=padarray(im,[40,40],0,'both');
    
    %make three colors
    imR=cat(3,im(:,:,1:length(fishidx)),im(:,:,length(fishidx)+2)+im(:,:,length(fishidx)+3)+im(:,:,length(fishidx)+4));
    imG=cat(3,im(:,:,1:length(fishidx)),im(:,:,length(fishidx)+2)+im(:,:,length(fishidx)+1)+im(:,:,length(fishidx)+4));
    imB=cat(3,im(:,:,1:length(fishidx)),im(:,:,length(fishidx)+3)+im(:,:,length(fishidx)+1)+im(:,:,length(fishidx)+4));
    
    cd('../../pregfp');
    mkdir(slicenames{i});
    minscore=min(score{i},[],2);
    for n=1:length(hroi{i})
        if minscore(n)>scorethresh&&rfp{i}(n)>rfpthresh %filter out cells that are unlikely to be successfully called.
            for m=1:size(imR,3)
                if m==1
                    cimR=imcrop(imR(:,:,1),[hroi{i}(n,1),hroi{i}(n,2),80,80]);
                    cimG=imcrop(imG(:,:,1),[hroi{i}(n,1),hroi{i}(n,2),80,80]);
                    cimB=imcrop(imB(:,:,1),[hroi{i}(n,1),hroi{i}(n,2),80,80]);
                else
                    cimR(:,:,m)=imcrop(imR(:,:,m),[hroi{i}(n,1),hroi{i}(n,2),80,80]);
                    cimG(:,:,m)=imcrop(imG(:,:,m),[hroi{i}(n,1),hroi{i}(n,2),80,80]);
                    cimB(:,:,m)=imcrop(imB(:,:,m),[hroi{i}(n,1),hroi{i}(n,2),80,80]);
                end
            end
            
            cimR=reshape(cimR,size(cimR,1),[],1);
            cimG=reshape(cimG,size(cimG,1),[],1);
            cimB=reshape(cimB,size(cimB,1),[],1);
            cim=cat(3,cimR,cimG,cimB);
            cim(:,81:81:81*size(imR,3),:)=500;  %add a white line on boundaries.
            imwrite(uint8(cim),[slicenames{i},'/',num2str(n),'.tif']);
        end
    end
    cd ..
end



