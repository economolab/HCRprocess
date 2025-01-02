% parse number of channels and number of planes in a tif

function write_tiff(path, im, info)
    
    t = Tiff(path, 'w8');

    sz = size(im);

    tagstruct.ImageLength = sz(1);
    tagstruct.ImageWidth = sz(2);
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt; % uint
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = sz(3);
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.ImageDescription = info.ImageDescription;
    tagstruct.ExtraSamples = Tiff.ExtraSamples.Unspecified;

    for ii=1:sz(4)
        plane = im(:,:,:,ii);
        setTag(t,tagstruct);
        write(t,uint16(plane));
        writeDirectory(t);
    end
    close(t)
    
end