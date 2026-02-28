% flip the z direction of a 3-D image (tif) using imageJ, expects a full 
% file path, will overwrite the original file with the same name
function filter_masks_ccf(masks_im,ccf_im,ccf_structs)
    
    ccf_vol = read_tiff(ccf_im);
    [masks_Vol, info] = read_tiff(masks_im);
    ccf_structs = split(ccf_structs,',');

    parcID = tax2parcID(ccf_structs);

    stats = regionprops3(masks_Vol,"Centroid");
    stats = round(stats);
    
    centroid_in_tax = zeros(height(stats),1,"logical");
    for i=1:height(stats)
        ind = table2array(stats(i,:));
        if ismember(ccf_vol(ind(2),ind(1),ind(3)),parcID)
            centroid_in_tax(i) = true;
        end
    end
    centroid_in_tax = [true; centroid_in_tax];
    
    cell_ids = unique(masks_Vol);
    new_cell_ids = cell_ids(centroid_in_tax);
    
    Livol = ismember(masks_Vol,new_cell_ids);
    newMasksVol = masks_Vol;
    newMasksVol(~Livol) = 0;

    [filepath,name,~] = fileparts(masks_im);
    uniq_id = split(name,'_');
    uniq_id = uniq_id{1};
    path = fullfile(filepath,[uniq_id '_masks_ccf.tif']);
    t = Tiff(path, 'w8');

    sz = size(newMasksVol);

    tagstruct.ImageLength = sz(1);
    tagstruct.ImageWidth = sz(2);
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt; % uint
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.ImageDescription = info.ImageDescription;
    tagstruct.ExtraSamples = Tiff.ExtraSamples.Unspecified;

    for ii=1:sz(3)
        plane = newMasksVol(:,:,ii);
        setTag(t,tagstruct);
        write(t,uint16(plane));
        writeDirectory(t);
    end
    close(t)

end