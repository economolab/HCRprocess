function writeTileConfig(path, names, series, coordStr, dim)

    fid = fopen(path, 'w');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '# Define the number of dimensions we are working on\n');
    fprintf(fid, 'dim = %d\n\n', dim);
    fprintf(fid, '# Define the image coordinates\n');

    for i = 1:numel(names)
        fprintf(fid, "%s; %s; (%s)\n", names(i), series(i), coordStr(i));
    end
    
end