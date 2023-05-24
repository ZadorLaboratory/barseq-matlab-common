%fix bleedthrough in all tif images in the current folder and align them.
%updated with 20x bleedthrough profile and channel shift correction
%4/8/2018. This version uses FFT phase correlation for image alignment and
%cuts off the border for alignment.
function tforms=Ultraviewseqalignment3(filename)
%read all tif file names in dirname folder
files=dir(fullfile(['*',filename,'*.tif']));
files={files.name};
L=size(files);
fixedfilelist=cell(L(2));
mkdir original;
% Fixes bleedthrough for all files, and store output filename list
parfor k=1:L(2)
    currentfile=files{k};
    fixedfilelist{k}=fixbleed(currentfile);
    
    movefile(files{k},['original/',files{k}]);
end

copyfile(fixedfilelist{1},['aligned',fixedfilelist{1}]);
for k=2:L(2)
    currentfile=fixedfilelist{k};
    tforms{k}=alignseq(currentfile,fixedfilelist{1});
end
%rgbout('bgnsub');
%movefile('RGB*.tif','rgb/');
%movefile('p*.tif','original/');

movefile('aligned*.tif','aligned/');
delete('fixed*tif');

%rgbout('aligned');
%movefile('RGB*.tif','rgb/');
%movefile('rough*.tif','rough/');
%movefile('aligned*.tif','aligned/');
end

% align all images in fixedfilelist to fixedfilelist(1).



    
    
% This subfunction fixes bleed through between channels in illumina seq tiff file.
    function fixedimage=fixbleed(filename)
        imageC=double(imread(filename,1));
        imageY=double(imread(filename,2));
        imageM=double(imread(filename,3));
        imageW=double(imread(filename,4));
        
        %correct for channel alignment
        imageC=imtranslate(imageC,[1.9 -0.4]);
        
        %median filter, updated to go before bleedthrough correction
        imageY=medfilt2(imageY);
        imageW=medfilt2(imageW);
        imageC=medfilt2(imageC);
        imageM=medfilt2(imageM);
        
        % Assign bleed through coefficients and bgn values
        %updated with new coeff from 20x SNAP25 sequencing images
        CtoY=0.05;
        YtoC=0.35;
        MtoY=0.02;
        MtoW=0.84;
        WtoM=0.05;
        %Cbgn=600;
        %Mbgn=900;
        %Ybgn=1300;
        %Wbgn=2900;
        
        YtoCcorrection=imageY*YtoC;
        CtoYcorrection=imageC*CtoY;
        MtoYcorrection=imageM*MtoY;
        MtoWcorrection=imageM*MtoW;
        WtoMcorrection=imageW*WtoM;
        imageC=max(imageC-YtoCcorrection,0);
        imageY=max(imageY-CtoYcorrection-MtoYcorrection,0);
        imageW=max(imageW-MtoWcorrection,0);
        imageM=max(imageM-WtoMcorrection,0);
        %imageY=uint16(imageY);
        %imageW=uint16(imageW);
        
        %bgn subtraction
        radius=20;
        ball=strel('ball', radius, radius);
        backgroundC=imopen(imageC, ball);
        backgroundY=imopen(imageY, ball);
        backgroundM=imopen(imageM, ball);
        backgroundW=imopen(imageW, ball);
        imageC=max(imageC-backgroundC,0);
        imageY=max(imageY-backgroundY,0);
        imageM=max(imageM-backgroundM,0);
        imageW=max(imageW-backgroundW,0);
        
        fixedimage=strcat('fixed',filename);
        outputfile=fixedimage;
        imwrite(uint16(imageC),outputfile);
        imwrite(uint16(imageY), outputfile, 'WriteMode','append');
        imwrite(uint16(imageM), outputfile, 'WriteMode','append');
        imwrite(uint16(imageW), outputfile, 'WriteMode','append');
        
    end

    function tform=alignseq(imagename,templatename)
        image1=imread(imagename,1);
        image2=imread(imagename,2);
        image3=imread(imagename,3);
        image4=imread(imagename,4);
        %correctedimage4=max(0,image4-1000)*1.2;
        imagesum=uint16(65000*double(image1+image2+image3+image4)./double(max(max(image1+image2+image3+image4))));
        template1=imread(templatename,1);
        template2=imread(templatename,2);
        template3=imread(templatename,3);
        template4=imread(templatename,4);
        templatesum=uint16(65000*double(template1+template2+template3+template4)./double(max(max(template1+template2+template3+template4))));

        tform = imregcorr(imagesum(1000:end-1000,1000:end-1000),templatesum(1000:end-1000,1000:end-1000),'translation','Window',0);
        Rfixed=imref2d(size(templatesum));
        alignedimage1=imwarp(image1,tform,'OutputView',Rfixed);
        alignedimage2=imwarp(image2,tform,'OutputView',Rfixed);
        alignedimage3=imwarp(image3,tform,'OutputView',Rfixed);
        alignedimage4=imwarp(image4,tform,'OutputView',Rfixed);                
        
        %write aligned image to file.
        alignedfile=strcat('aligned',imagename);
        imwrite(uint16(alignedimage1), alignedfile);
        imwrite(uint16(alignedimage2),alignedfile, 'WriteMode','append');
        imwrite(uint16(alignedimage3),alignedfile, 'WriteMode','append');
        imwrite(uint16(alignedimage4),alignedfile, 'WriteMode','append');
    end
    
    

        
        
        
        
