% check if a directory exists, and if it doesn't, make it

function check_build_dir(dir_path)

    list_dir = dir(dir_path);
    if isempty(list_dir)
        mkdir(dir_path)
    end
    
end