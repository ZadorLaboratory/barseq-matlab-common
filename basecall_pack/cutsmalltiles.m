function cutsmalltiles(filename)
%cut tif files into smaller tiles for Starfish. filename is the prefix of
%folders to be processed.

files=dir([filename,'*']);
mkdir('smalltiles');
for i=1:length(files)
    imfiles=dir([files(i).name,'/*.tif']);
    for m=1:length(imfiles)
        im=imread([files(i).name,'/',imfiles(m).name]);
        if i==1 && m==1 
            tilesz=round(size(im)./5);
            tileedges=(0:5)'*tilesz;
            tileedges(end,:)=size(im);
        end
        
        
        
        %cut into small pieces
        for p=1:5
            for q=1:5
                imwrite(im((tileedges(p,1)+1):tileedges(p+1,1),(tileedges(q,2)+1):tileedges(q+1,2)), ...
                    ['smalltiles/seq',num2str(i),imfiles(m).name(1:end-4),'X',num2str(p),'Y',num2str(q),'.tif']);
            end
        end
    end
end

        
        
    
    



