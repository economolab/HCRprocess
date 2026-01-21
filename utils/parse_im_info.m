% parse number of channels and number of planes in a tif

function [chans, planes, info, data_found_dict] = parse_im_info(path)
    
    info = imfinfo(path);
    data_found_dict = dictionary;
    
    if isfield(info,'ImageDescription')

        id = info.ImageDescription;

        tok = regexp(id, 'channels=(\d*)', 'tokens');
        if isempty(tok)
            disp('Could not read number of channels in tif.  Assuming 1 channel');
            chans = 1;
            data_found_dict('chans_found') = false;
        else
            chans = str2double(tok{1}{1});
            data_found_dict('chans_found') = true;
        end
        
        tok = regexp(id, 'slices=(\d*)', 'tokens');
        if isempty(tok)
            tok = regexp(id, 'frames=(\d*)', 'tokens');
        end
    
        if isempty(tok)
            disp('Could not read number of z planes in tif.  Assuming planes = number of IFDs');
            planes = length(info);
            data_found_dict('planes_found') = false;
        else
            planes = str2double(tok{1}{1});
            data_found_dict('planes_found') = true;
        end

    else

        disp(['Could not find Image Description field in metadata, assuming 1 channel and' ...
            ' planes = number of IFDs'])
        
        chans = 1;
        data_found_dict('chans_found') = false;
        planes = length(info);
        data_found_dict('planes_found') = false;

    end
    
end