%fix bleedthrough in all tif images in the current folder and align them.
%updated with 20x bleedthrough profile and channel shift correction
%4/8/2018
function Ultraviewfiximage(filename)
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
%find rolonies and align
mkdir RGB;
rgbout1('fixed',4);
movefile('RGB*.tif','RGB/');



%mkdir fixed;
%movefile('fixed*.tif','fixed/');
end
    
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

    

        
        
        
        
