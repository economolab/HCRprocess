% read multi-z tiff
% If target channel is specified, returns image of dimensions: 
% (height, width, planes)
% Otherwise, returns image of dimensions:
% (height, width, channels, planes)

function [im, info] = read_tiff(path,targetChan)

    arguments
        path
        targetChan = []
    end

    [chans, planes, info] = parse_im_info(path);
    
    % single IFD for the whole image, saved by ImageJ
    if length(info) == 1

        im = tiffreadVolume(path);
        sz = size(im);
        im = reshape(im, sz(1), sz(2), chans, planes);

    % IFD for every plane, saved by MATLAB
    else

        im = zeros(info(1).Height,info(1).Width,chans,planes);
        warning('off','imageio:tiffutils:libtiffWarning')
        t = Tiff(path,"r");
        
        % read 1 plane at a time
        for i=1:planes
            
            im(:,:,:,i) = read(t);

            % increment directory, unless you've just loaded the last plane
            if i < planes
                nextDirectory(t)
            end

        end
        
    end
    
    % restrict output image to target channel, if it was provided
    if ~isempty(targetChan)
        im = im(:,:,targetChan,:);
        im = squeeze(im);
    end
 
end