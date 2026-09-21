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

if isempty(app.seg_sel)
    message = ['No file selected for segmentation. Click on the ' ...
        'table to select.'];
    title = 'No segmentation file selection';
    uialert(app.UIFigure,message,title)
    return
end

% fetch the tifs from the core output folder for the
% segmentation selection
im_dir = fullfile(app.exp_dir,'post','core_output',app.seg_sel);
listing = dir(im_dir);
f_names = {listing.name};
TF_im = contains(f_names,{'.tif','.tiff'},'IgnoreCase',true);
im_files = f_names(TF_im);

% figure out which are registration images
TF_reg_channels = contains(im_files,app.reg_channel_token,'IgnoreCase',true);
TF_bin = ~TF_reg_channels;
f_label_channels_ims = im_files(TF_bin);

[indx,~] = listdlg('ListString',f_label_channels_ims, ...
    'SelectionMode','multiple', ...
    'PromptString',{'Select label channels to segment','(Ctrl+Click to select multiple)'});
f_label_channels_ims = f_label_channels_ims(indx);

% segment the registration image (same as normal segmentation)
input_dir = fullfile(app.exp_dir,'post','core_output',app.seg_sel);
listing = fetch_dir(input_dir,'file');
split_str = split(app.reg_round,'_');
roundstr = split_str{2};
strs_filt = filt_kw(listing,roundstr,'keep');
strs_filt = filt_kw(strs_filt,app.reg_channel_token,'keep');
input_f = fullfile(input_dir,strs_filt{1});
app.seg_image(input_f)

f_label_ims_masks = cell(length(f_label_channels_ims),1);
% segment the label images
for i=1:length(f_label_channels_ims)
    app.seg_image(fullfile(im_dir,f_label_channels_ims{i}),strcat('_',f_label_channels_ims{i},'_masks.tif'))
    f_label_ims_masks{i} = strcat('_',f_label_channels_ims{i},'_masks.tif');
end

% get the masks directory for this section
[filepath,~,~] = fileparts(input_f);
[filepath,uniq_id,~] = fileparts(filepath);
[filepath,~,~] = fileparts(filepath);
masks_dir = fullfile(filepath,'masks',uniq_id);

f_masks_im = fullfile(masks_dir,strcat(uniq_id,'_cp_masks.tif'));
for i=1:length(f_label_ims_masks)
    f_label_ims_masks{i} = fullfile(masks_dir,strcat(uniq_id,f_label_ims_masks{i}));
end

% to do: figure out all the file naming bs here
processLabelSegmentation(f_masks_im,f_label_ims_masks)

