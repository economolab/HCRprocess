% find the full filename of an image given its uniq_id and target directory

function fullfn = find_fn_uniq_id(uniq_id,targetDir)

    listing = fetch_dir(targetDir,'file');
    
    % make sure there is a file in the target directory that contains the
    % unique id
    if sum(contains(listing,uniq_id)) == 1
        fn = listing{contains(listing,uniq_id)};
        fullfn = fullfile(targetDir,fn);
    else
        fullfn = [];
    end
  
end