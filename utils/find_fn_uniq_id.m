% find the full filename of an image given its uniq_id and target directory

function fullfn = find_fn_uniq_id(uniq_id,targetDir)

    listing = fetch_dir(targetDir,'file');
    fn = listing{contains(listing,uniq_id)};
    fullfn = fullfile(targetDir,fn);
  
end