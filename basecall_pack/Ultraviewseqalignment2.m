%fix bleedthrough in all tif images in the current folder and align them.
%updated with 20x bleedthrough profile and channel shift correction
%4/8/2018
function Ultraviewseqalignment2(filename)
%read all tif file names in dirname folder
files=dir(fullfile(['*',filename,'*.tif']));
files={files.name};
L=size(files);
fixedfilelist=cell(L(2));
mkdir original;
% Fixes bleedthrough for all files, and store output filename list
parpool(5)
parfor k=1:L(2)
    currentfile=files{k};
    fixedfilelist{k}=fixbleed(currentfile);
    
    movefile(files{k},['original/',files{k}]);
end

copyfile(fixedfilelist{1},['aligned',fixedfilelist{1}]);
parfor k=2:L(2)
    currentfile=fixedfilelist{k};
    alignseq(currentfile,fixedfilelist{1});
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

    function alignseq(imagename,templatename)
        image1=imread(imagename,1);
        image2=imread(imagename,2);
        image3=imread(imagename,3);
        image4=imread(imagename,4);
        %correctedimage4=max(0,image4-1000)*1.2;
        imagesum=uint16(65000*double(image1*5+image2+image3+image4)./double(max(max(image1*5+image2+image3+image4))));
        template1=imread(templatename,1);
        template2=imread(templatename,2);
        template3=imread(templatename,3);
        template4=imread(templatename,4);
        templatesum=uint16(65000*double(template1*5+template2+template3+template4)./double(max(max(template1*5+template2+template3+template4))));

        % set alignment conditions and align with iat, generate
        % alignedimage
        %[d1, l1]=iat_surf(imagesum);
        %[d2, l2]=iat_surf(templatesum);
        %d1=detectSURFFeatures(imagesum,'MetricThreshold',30);
        %d2=detectSURFFeatures(templatesum,'MetricThreshold',30);
        %[~, ~, imgInd, tmpInd]=iat_match_features_mex(d1,d2,ratio);
        %X1 = l1(imgInd,1:2);
        %X2 = l2(tmpInd,1:2);
        %iat_plot_correspondences(imagesum, templatesum, X1', X2');
        %X1h = iat_homogeneous_coords (X1');
        %X2h = iat_homogeneous_coords (X2');
        %[~, ransacWarp]=iat_ransac( X2h, X1h,'translation');
        %iat_plot_correspondences(imagesum,templatesum,X1(inliers,:)',X2(inliers,:)');
        par.transform = 'euclidean'; 
        par.levels = 7;
        par.iterations = 100;
        ransacWarp=iat_ecc(imagesum(:,:),templatesum(:,:),par);
        
        [M,N]=size(template4);
        [alignedimage1,~]=iat_inverse_warping(image1,ransacWarp,par.transform,1:N, 1:M);
        [alignedimage2,~]=iat_inverse_warping(image2,ransacWarp,par.transform,1:N, 1:M);
        [alignedimage3,~]=iat_inverse_warping(image3,ransacWarp,par.transform,1:N, 1:M);
        [alignedimage4,~]=iat_inverse_warping(image4,ransacWarp,par.transform,1:N, 1:M);
        %write aligned image to file.
        alignedfile=strcat('aligned',imagename);
        imwrite(uint16(alignedimage1), alignedfile);
        imwrite(uint16(alignedimage2),alignedfile, 'WriteMode','append');
        imwrite(uint16(alignedimage3),alignedfile, 'WriteMode','append');
        imwrite(uint16(alignedimage4),alignedfile, 'WriteMode','append');
    end
    
    

        
        
        
        
