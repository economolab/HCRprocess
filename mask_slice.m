vol = read_tiff('D:\2025-10-18_contra_1\post\core_output\s04\s04annotation.tif');
[imVol, info] = read_tiff('D:\2025-10-18_contra_1\post\core_output\s04\r1_neur445.tif');

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

function grid = makeGrid(im,rects,grid_line_thickness)

    arguments
        im
        rects
        grid_line_thickness = 3
    end
    
    half_thickness = round(grid_line_thickness / 2);
    im_height = size(im,1);
    im_width = size(im,2);
    grid = zeros(im_height,im_width,"logical");
    rect_width = im_width / rects;
    rect_height = im_height / rects;

    for i = 1:rects-1
        x1 = uint16(i*rect_height - half_thickness);
        x2 = uint16(i*rect_height + half_thickness);
        y1 = uint16(i*rect_width - half_thickness);
        y2 = uint16(i*rect_width + half_thickness);
        grid(x1:x2,:) = 1;
        grid(:,y1:y2) = 1;
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

%% overlay just IRN/PARN and grid

tax = {'IRN', 'PARN'};
parcID = tax2parcID(tax);

grid = makeGrid(imVol,10,5);
Livol = ismember(vol,parcID);
volMask = zeros(size(Livol),'logical');
for i=1:size(Livol,3)
    BW = bwperim(Livol(:,:,i));
    BW = bwmorph(BW,"thicken",2);
    BW = bwmorph(BW,"bridge");
    BW = bwmorph(BW,"majority",2);
    BW = BW | grid;
    volMask(:,:,i) = BW;
end

volExport = zeros(size(imVol), "uint16");
for i=1:size(imVol,3)
    overlayIntensity = max(imVol(:,:,i),[],"all");
    overlayPlane = zeros(size(imVol,1),size(imVol,2),"uint16");
    overlayPlane(volMask(:,:,i)) = overlayIntensity;
    volExport(:,:,i) = imblend(overlayPlane,imVol(:,:,i),ForegroundOpacity=0.1);
end

%% overlay just anno boundaries

volMask = zeros(size(vol),'logical');
for i=1:size(vol,3)
    BW = edge(vol(:,:,i));
    BW = bwmorph(BW,"thicken",2);
    BW = bwmorph(BW,"bridge");
    BW = bwmorph(BW,"majority",2);
    volMask(:,:,i) = BW;
end

volExport = zeros(size(imVol), "uint16");
for i=1:size(imVol,3)
    overlayIntensity = max(imVol(:,:,i),[],"all");
    overlayPlane = zeros(size(imVol,1),size(imVol,2),"uint16");
    overlayPlane(volMask(:,:,i)) = overlayIntensity;
    volExport(:,:,i) = imblend(overlayPlane,imVol(:,:,i),ForegroundOpacity=0.1);
end

%%

writeTiffOverlay('s04_with_full_anno_overlay.tif',volExport,info)
