function recountfish(hybthresh1)
%reset hybthresh and count again with prevously drawn cells.

%%
folders1=dir('processed*');
folders1=sort_nat({folders1.name});
%%
for o=1:numel(folders1)
    %%
    cd(folders1{o});
    hybthresh=hybthresh1;
    load temp1
    hybthresh1=hybthresh;
    %% find dots and count.
    origfiles=dir('origcrop*.tif');
    origfiles=sort_nat({origfiles.name});
%     if ~exist('hybthresh1','var')
%         hybthresh1=[30 30 30 30];
%     end
    
    rols={};
    filtrols={};
    recounts=[];
    for m=1:numel(files)
        if ~isempty(cellshapeconvex{m})
            r={};rsub=r;rsubfilt=rsub;
            for n=1:4
                a=imread(origfiles{m},n);
                CC = bwconncomp(imregionalmax(imreconstruct(max(a-hybthresh1(n),0),a)));
                r{n}=zeros(length(CC.PixelIdxList),1);
                for i=1:length(CC.PixelIdxList)
                    [~,I]=max(a(CC.PixelIdxList{i}));
                    r{n}(i)=CC.PixelIdxList{i}(I); %linear indexed peak positions
                end
                [y,x]=ind2sub(size(a),r{n});
                rsub{n}=[x,y];%rsub is in x, y
                %check if dots are within the convex hull
                rsubfilt{n}=rsub{n}(inpolygon(rsub{n}(:,1),rsub{n}(:,2),cellshape{m}(cellshapeconvex{m},1),cellshape{m}(cellshapeconvex{m},2)),:);
                recounts(m,n)=size(rsubfilt{n},1);
            end
            rols{m}=rsub;
            filtrols{m}=rsubfilt;
        else
            rols{m}={};
            filtrols{m}={};
            recounts(m,1:4)=0;
        end
    end
    
    save('results-recount.mat','rols','filtrols','recounts','cellshape','cellshapeconvex','cellpos','contours','depths','angles','ctb');
    %%
    cd ..
end