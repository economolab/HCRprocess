masks_im_f = 'D:\2026-09-12_FN_SNr_2\post\masks\s01L\s01L_masks_qc.tif';
bin_genes_table_f = 'D:\2026-09-12_FN_SNr_2\post\core_output\s01L\s01L_bin_genes.txt';
[masksDir,~,~] = fileparts(masks_im_f);
[core_dir,~,~] = fileparts(bin_genes_table_f);

bin_genes_table = readtable(bin_genes_table_f);
masks_table = readtable('D:\2026-09-12_FN_SNr_2\post\core_output\s01L\s01L_final_qc_masks.txt');
masks_im = read_tiff(masks_im_f);

label_genes = {'egfp', 'dtom'};
good_masks = strcmp(masks_table.('Result'),'good');

labeled_masks_bool = zeros(length(good_masks),1,'logical');
for i=1:length(label_genes)
    labeled_masks_bool = (labeled_masks_bool | bin_genes_table.(label_genes{i}));
end

labeled_masks_bool = (labeled_masks_bool & good_masks);

sub_bin_genes_table = bin_genes_table(labeled_masks_bool,:);
sub_bin_genes_table.Properties.VariableNames{1} = 'Mask ID';
sub_bin_genes_table.('Labels Assessed') = repmat("null", height(sub_bin_genes_table), 1);
sub_bin_genes_table.('Genes Assessed') = repmat("null", height(sub_bin_genes_table), 1);
[~,table_f,~] = fileparts(bin_genes_table_f);
table_f = strcat(table_f, '_labeled.xls');
writetable(sub_bin_genes_table, fullfile(core_dir, table_f))

masks = masks_table.('MaskID');
labeled_masks = masks(labeled_masks_bool);
masks_im(~ismember(masks_im,labeled_masks)) = 0;

t = Tiff(fullfile(masksDir,'check_labeled_masks.tiff'), 'w8');
sz = size(masks_im);

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
    plane = masks_im(:,:,ii);
    setTag(t,tagstruct);
    write(t,uint16(plane));
    writeDirectory(t);
end
close(t)





