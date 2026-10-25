
function integrate_masks(manual_masks_f, cp_masks_f)
    
    manual_masks = read_tiff(manual_masks_f);
    cp_masks = read_tiff(cp_masks_f);
    
    max_manual_mask = max(manual_masks,[],"all");
    
    manual_masks_nonzero = (manual_masks > 0);
    cp_masks_nonzero = (cp_masks > 0);
    
    cp_masks(cp_masks_nonzero) = cp_masks(cp_masks_nonzero) + max_manual_mask;
    
    hit_cp_masks = unique(cp_masks(manual_masks_nonzero));
    hit_cp_masks(hit_cp_masks == 0) = [];
    
    cp_masks(ismember(cp_masks, hit_cp_masks)) = 0;
    cp_masks = cp_masks + manual_masks;
    
    [masksDir,~,~] = fileparts(cp_masks_f);
    t = Tiff(fullfile(masksDir,'integrated_masks.tiff'), 'w8');
    
    sz = size(cp_masks);
    
    tagstruct.ImageLength = sz(1);
    tagstruct.ImageWidth = sz(2);
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt; % uint
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.ExtraSamples = Tiff.ExtraSamples.Unspecified;
    
    for ii=1:sz(3)
        plane = cp_masks(:,:,ii);
        setTag(t,tagstruct);
        write(t,uint16(plane));
        writeDirectory(t);
    end
    close(t)

end