% fetch a cell array of strings containing all the files or subfolders
% in a directory
% type = 'file' for files, type = 'dir' for subfolders

function listing = fetch_dir(path, type)

    listing = dir(path);
    listing = listing(~ismember({listing.name}, {'.','..'}));

    switch type
        case 'file'
            listing = listing(~[listing.isdir]);
        case 'dir'
            listing = listing([listing.isdir]);
    end

    listing = {listing.name};

end