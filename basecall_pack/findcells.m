function [hroi,lroi,tr,tt,tforms]=findcells(slicenames,hthresh,lthresh)
%align 20x to 10x (2xbin) 1st cycle sequencing images.
%This version assumes that highres/lowres images were in the same orientation.
addpath('C:\Fiji.app\scripts');
javaaddpath 'C:\Program Files\MATLAB\R2018b\java\mij.jar'
Miji(false);

for i=1:length(slicenames)
    
    him=imread(['fish/aligned',slicenames{i},'highres.tif'],1);
    lfiles=dir(['10x/',slicenames{i},'/aligned/alignedfixed*seq01.tif']);
    lim=imread(['10x/',slicenames{i},'/aligned/',lfiles(1).name],1);
    
    for n=2:4
        him(:,:,n)=imread(['fish/aligned',slicenames{i},'highres.tif'],n);
        lim(:,:,n)=imread(['10x/',slicenames{i},'/aligned/',lfiles(1).name],n);
    end
    him=imgaussfilt(imresize(him,0.504),1);
    
    %manually rough align images
    %[prevpts,currpts]=cpselect(lim(:,:,1)*500, ...
    %him(:,:,1)*200,'Wait',true);
    %tforms{i}=fitgeotrans(currpts,prevpts,'nonreflectivesimilarity');
    
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.InitialRadius = optimizer.InitialRadius/3.5;
    himm = imregister(him(200:end-200,200:end-200,3), lim(200:end-200,200:end-200,3), 'translation', optimizer, metric);
    tforms{i} = imregcorr(him(200:end-200,200:end-200,3),himm,'translation','Window',0);
    
    %check alignment
    him1=[];
    %
        him1=imwarp(him(:,:,3),tforms{i},...
            'outputview',imref2d([size(lim,1),size(lim,2)])); 
    
    %    figure;imshowpair(him1*200,lim(:,:,3)*500);
        imwrite(uint8(cat(3,him1,lim(:,:,3)*2.5,zeros(size(lim,1),size(lim,2)))),['fish/seq2seq',slicenames{i},'ch1.tif']);
        
        
        
    % find rolonies using Miji  in matlab
        hrol=[];
        %find rolonies
        for n=1:4
            MIJ.createImage(him(:,:,n));
            MIJ.run('Find Maxima...', ['noise=',num2str(hthresh),' output=List']);
            hrol=[hrol;MIJ.getResultsTable]; %[x,y]
            %cleanup
            MIJ.run('Clear Results');
            MIJ.closeAllWindows
        end
        lrol=[];
        %find rolonies
        for n=1:4
            MIJ.createImage(lim(:,:,n));
            MIJ.run('Find Maxima...', ['noise=',num2str(lthresh),' output=List']);
            lrol=[lrol;MIJ.getResultsTable];
            %cleanup
            MIJ.run('Clear Results');
            MIJ.closeAllWindows
        end
        

        %combine rolonies within 7 pixels of each other
        hpeaks=zeros(size(him,1),size(him,2));
        hpeaks(sub2ind([size(him,1),size(him,2)],hrol(:,2)+1,hrol(:,1)+1))=1;
        hpeaks=hpeaks&~imdilate(hpeaks,triu(ones(15))-diag([ones(1,8),zeros(1,7)]));
        [hidxy,hidxx]=find(hpeaks);
        
        lpeaks=zeros(size(lim,1),size(lim,2));
        lpeaks(sub2ind([size(lim,1),size(lim,2)],lrol(:,2)+1,lrol(:,1)+1))=1;
        lpeaks=lpeaks&~imdilate(lpeaks,triu(ones(15))-diag([ones(1,8),zeros(1,7)]));
        [lidxy,lidxx]=find(lpeaks);
        
        %transform high-res rolonies to rough aligned image
        hrol1=[hidxx,hidxy, ones(length(hidxy),1)]*tforms{i}.T;
       
        
        %align rolonies using icp
        [tr{i},tt{i}]=icp([[lidxx,lidxy]';ones(1,length(lidxy))],hrol1');
        hroi{i}=round([hidxy,hidxx]/0.504);
        lroi{i}=(tr{i}*hrol1'+ tt{i})';
        lroi{i}(:,3)=[];
        
end
    save('rolonies.mat','hroi','lroi','tt','tr','tforms');
    
 MIJ.exit; 
end

