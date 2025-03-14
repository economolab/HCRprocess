function stitch(file, d, file_frac)
   
    java.lang.Runtime.getRuntime.gc;
    
    ImageJ

    ij.IJ.run("Grid/Collection stitching", "type=[Positions from file] " + ...
        "order=[Defined by image metadata] " + ...
        "browse=" + ...
        file + ...
        " multi_series_file=" + ...
        file + ...
        " fusion_method=[Max. Intensity] " + ...
        "regression_threshold=0.30 " + ...
        "max/avg_displacement_threshold=2.50 " + ...
        "absolute_displacement_threshold=3.50 " + ...
        "compute_overlap increase_overlap=0 " + ...
        "invert_x " + ...
        "computation_parameters=[Save memory (but be slower)] " + ...
        "image_output=[Fuse and display]");

    [filepath,name,~] = fileparts(file);
    savef = filepath(1:end-4);
    savef = fullfile(savef,'stitch',name);
    savef = strcat(savef,'.tif');

    ij.IJ.saveAs("Tiff", savef);
    ij.IJ.run("Close All");
    ij.IJ.run("Quit","");
    
    java.lang.Runtime.getRuntime.gc;

    junk = fullfile(filepath,'TileConfiguration.registered.txt');
    delete(junk)

    d.Message = strcat(file_frac, ' Deleting empty slices...');

    [V, info] = read_tiff(savef);

    sz = size(V);
    xy_size = sz(1) * sz(2) * sz(3);
    maskZ = ones(1, sz(4),'logical');

    for i=1:length(maskZ)
        if (sum(V(:,:,:,i) == 0, "all") / xy_size) > 0.5
            maskZ(i) = 0;
        end
    end
    
    V = V(:,:,:,maskZ);
    sz = size(V);

    id = info.ImageDescription;

    replace = 'slices=' + string(sum(maskZ));
    id = regexprep(id, 'slices=(\d*)', replace);
    replace = 'images=' + string(sz(3)*sz(4));
    id = regexprep(id, 'images=(\d*)', replace);

    [info.ImageDescription] = deal(id);

    d.Message = strcat(file_frac, ' Saving...');

    write_tiff(savef, V, info)
    % write_tiff_ij(savef)

end