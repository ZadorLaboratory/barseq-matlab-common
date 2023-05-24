function convert_RGB_jpg()
%compress tif RGB to jpeg, then remove the original tif

folders=get_folders();
parfor i=1:numel(folders)
    imfiles=dir(fullfile(folders{i},'RGB','*.tif'));
    if numel(imfiles)>0
        imfiles={imfiles.name};
        for n=1:numel(imfiles)
            im=imread(fullfile(folders{i},'RGB',imfiles{n}));
            imwrite(im,fullfile(folders{i},'RGB',[imfiles{n}(1:end-3),'jpg']));
        end
        delete(fullfile(folders{i},'RGB','*.tif'));
    end
end