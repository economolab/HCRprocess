% parse number of channels and number of planes in a tif

function [im, info] = read_tiff(path)
    
    [chans, planes, info] = parse_im_info(path);
    im = tiffreadVolume(path);
    sz = size(im);
    im = reshape(im, sz(1), sz(2), chans, planes);
    
end