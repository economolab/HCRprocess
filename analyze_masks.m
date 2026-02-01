vol = read_tiff('D:\2025-10-18_contra_1\post\core_output\s04\s04annotation.tif');
[imVol, info] = read_tiff('C:\Users\jpv88\Documents\GitHub\HCRprocess\s04_with_overlay_cp_masks.tif');
% [imVol, info] = read_tiff('D:\2025-10-18_contra_1\post\core_output\s04\s04_nobg_smooth_cp_masks.tif');

%%

tax = {'IRN', 'PARN'};
parcID = tax2parcID(tax);

stats = regionprops3(imVol,"Centroid");
stats = round(stats);

centroid_in_tax = zeros(height(stats),1,"logical");
for i=1:height(stats)
    ind = table2array(stats(i,:));
    if ismember(vol(ind(2),ind(1),ind(3)),parcID)
        centroid_in_tax(i) = true;
    end
end
centroid_in_tax = [true; centroid_in_tax];

cell_ids = unique(imVol);
new_cell_ids = cell_ids(centroid_in_tax);

Livol = ismember(imVol,new_cell_ids);
newImVol = imVol;
newImVol(~Livol) = 0;

%%
writeTiffOverlay('s04_IRNPARN_masks.tif',newImVol,info)

%%

volMask = zeros(size(vol),'logical');
for i=1:size(vol,3)
    BW = edge(vol(:,:,i));
    BW = bwmorph(BW,"thicken",2);
    BW = bwmorph(BW,"bridge");
    BW = bwmorph(BW,"majority",2);
    volMask(:,:,i) = BW;
end

volExport = zeros(size(newImVol), "uint16");
for i=1:size(newImVol,3)
    overlayIntensity = max(newImVol(:,:,i),[],"all");
    overlayPlane = zeros(size(newImVol,1),size(imVol,2),"uint16");
    overlayPlane(volMask(:,:,i)) = overlayIntensity;
    volExport(:,:,i) = imblend(overlayPlane,newImVol(:,:,i),ForegroundOpacity=0.1);
end


%%
writeTiffOverlay('s04_filt_mask_CPSAM_curated.tif',volExport,info)

%%

tax = {'IRN', 'PARN'};
parcID = tax2parcID(tax);

%%

Livol = ismember(vol,parcID);
imVol(~Livol) = 0;

tax_cell_ids = unique(imVol);
new_cell_ids = transpose(0:length(tax_cell_ids)-1);
new_cell_ids = uint16(new_cell_ids);
d = dictionary(tax_cell_ids,new_cell_ids);
newImVol = zeros(size(imVol),"uint16");

for i=1:length(tax_cell_ids)
    newImVol(imVol == tax_cell_ids(i)) = d(tax_cell_ids(i));
    disp(i)
end

%%

stats = regionprops3(newImVol,"Centroid");
stats = round(stats);

centroid_in_tax = zeros(height(stats),1,"logical");
for i=1:height(stats)
    ind = table2array(stats(i,:));
    if ismember(vol(ind(1),ind(2),ind(3)),parcID)
        centroid_in_tax(i) = true;
    end
end

centroid_in_tax = [true; centroid_in_tax];
cell_id_to_keep = new_cell_ids(centroid_in_tax);

Livol = ismember(newImVol,cell_id_to_keep);
finalImVol = newImVol;
finalImVol(~Livol) = 0;
%%
writeTiffOverlay('s04_filt_mask_raw.tif',finalImVol,info)

%% 
function parcID = tax2parcID(tax,tax_level)

    arguments
        tax
        tax_level = 'structure'
    end

    abcDir = params('abcDir');
    colorDir = fullfile(abcDir,'metadata','Allen-CCF-2020','20230630','views');
    ccfTerms = fullfile(colorDir,'parcellation_to_parcellation_term_membership_acronym.csv');
    Tterms = readtable(ccfTerms);
    
    parcID = zeros(length(tax),1);
    for i=1:length(tax)
        mask = strcmp(Tterms.(tax_level),tax{i});
        parcID(i) = Tterms.parcellation_index(mask); 
    end

end

function writeTiffOverlay(path, im, info)

    t = Tiff(path, 'w8');

    sz = size(im);

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
        plane = im(:,:,ii);
        setTag(t,tagstruct);
        write(t,uint16(plane));
        writeDirectory(t);
    end
    close(t)

end