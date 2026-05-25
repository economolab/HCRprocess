function processLabelSegmentation(f_masks_im, f_label_ims)
    
    % threshold: what fraction of voxels in each mask need to be present in
    % the masks of at least one of the labeled images for the mask to be
    % kept
    thresh = 0.05;

    [masks_im, info] = read_tiff(f_masks_im);

    num_label_ims = length(f_label_ims);
    label_ims = cell(num_label_ims,1);
    for i=1:num_label_ims
        label_ims{i} = read_tiff(f_label_ims{i});
        label_ims{i} = (label_ims{i} > 0);
    end

    stats = regionprops3(masks_im,"Volume");
    num_masks = height(stats);
    for i=1:num_label_ims
        stats.(strcat('LabelIm',num2str(i))) = NaN(num_masks,1);
        stats.(strcat('LabelIm',num2str(i),'_pass')) = NaN(num_masks,1);
    end

    masks_lin = masks_im(:);
    masks_lin_bin = (masks_lin > 0);
    masks_lin = masks_lin(masks_lin_bin);

    for i=1:num_label_ims
        label_im_lin = label_ims{i}(:);
        label_im_lin = label_im_lin(masks_lin_bin);
        stats.(strcat('LabelIm',num2str(i))) = accumarray(masks_lin,label_im_lin);
    end
    
    for i=1:num_label_ims
        fracs = stats.(strcat('LabelIm',num2str(i))) ./ stats.('Volume');
        stats.(strcat('LabelIm',num2str(i),'_pass')) = (fracs >= thresh);
    end

    pass = zeros(num_masks,1,'logical');
    for i=1:num_label_ims
        pass = pass | stats.(strcat('LabelIm',num2str(i),'_pass'));
    end

    masks = unique(masks_im);
    masks(1) = [];
    keep_masks = masks(pass);
    keep_masks_mask = ismember(masks_im,keep_masks);
    masks_im(~keep_masks_mask) = 0;

    [filepath,name,~] = fileparts(f_masks_im);
    uniq_id = split(name,'_');
    uniq_id = uniq_id{1};
    path = fullfile(filepath,[uniq_id '_cp_masks_label_filt.tif']);
    t = Tiff(path, 'w8');

    sz = size(masks_im);

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
        plane = masks_im(:,:,ii);
        setTag(t,tagstruct);
        write(t,uint16(plane));
        writeDirectory(t);
    end
    close(t)

end

