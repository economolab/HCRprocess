function stitch(file, shadingCorrection, fusionMethod, stripeCorrection, d, file_frac)
    
    arguments
        file
        shadingCorrection = false
        fusionMethod = 'Max. Intensity'
        stripeCorrection = []
        d = []
        file_frac = []
    end
    
    if shadingCorrection
        pyenv('ExecutionMode', 'OutOfProcess');
        escapedFile = strrep(file, '\', '\\');
        inputStr = strcat(escapedFile, ',', stripeCorrection);
        cmd = sprintf("shadingCorrection.py '%s'", inputStr);
        pyrunfile(cmd)
        [filepath, name, ~] = fileparts(file);
        file = fullfile(filepath, [name, '.tiff']);
    end

    java.lang.Runtime.getRuntime.gc;

    switch fusionMethod
        case 'Max. Intensity'
            fusion_str = " fusion_method=[Max. Intensity] ";
        case 'Linear Blending'
            fusion_str = " fusion_method=[Linear Blending] ";
    end

    ImageJ

    ij.IJ.run("Grid/Collection stitching", "type=[Positions from file] " + ...
        "order=[Defined by image metadata] " + ...
        "browse=" + ...
        file + ...
        " multi_series_file=" + ...
        file + ...
        fusion_str + ...
        "regression_threshold=0.30 " + ...
        "max/avg_displacement_threshold=2.50 " + ...
        "absolute_displacement_threshold=3.50 " + ...
        "compute_overlap subpixel_accuracy increase_overlap=0 " + ...
        "invert_x " + ...
        "computation_parameters=[Save memory (but be slower)] " + ...
        "image_output=[Fuse and display]");

    % if shadingCorrection
    % 
    %     [filepath,name,ext] = fileparts(file);
    %     file_reg_ch = fullfile(filepath,strcat(name,'_reg_ch',ext));
    % 
    %     ImageJ
    %     ij.IJ.run("Grid/Collection stitching", "type=[Positions from file] " + ...
    %         "order=[Defined by image metadata] " + ...
    %         "browse=" + ...
    %         file_reg_ch + ...
    %         " multi_series_file=" + ...
    %         file_reg_ch + ...
    %         fusion_str + ...
    %         "regression_threshold=0.30 " + ...
    %         "max/avg_displacement_threshold=2.50 " + ...
    %         "absolute_displacement_threshold=3.50 " + ...
    %         "compute_overlap subpixel_accuracy increase_overlap=0 " + ...
    %         "invert_x " + ...
    %         "computation_parameters=[Save memory (but be slower)] " + ...
    %         "image_output=[Fuse and display]");
    % 
    %     ij.IJ.run("Close All");
    %     ij.IJ.run("Quit","");
    % 
    %     file_in  = fullfile(filepath, 'TileConfiguration.registered.txt');
    %     file_out = fullfile(filepath, 'TileConfiguration.txt');   % don't overwrite the registered one
    % 
    %     [~, ~, coordStr, dim] = readTileConfig(file_in);
    %     n      = numel(coordStr);
    %     names  = compose("tile_%02d.tif", (0:n-1).');
    %     series = strings(n,1);   % empty series field
    % 
    %     writeTileConfig(file_out, names, series, coordStr, dim);
    % 
    %     for i = 0:n-1   % 0-based series index, matches your tile index column
    %         opts = loci.plugins.in.ImporterOptions();
    %         opts.setId(file);
    %         opts.setSeriesOn(i, true);
    %         if i ~= 0, opts.setSeriesOn(0, false); end   % series 0 is on by default
    % 
    %         imps = loci.plugins.BF.openImagePlus(opts);
    %         imp  = imps(1);
    %         ij.IJ.saveAs(imp, 'Tiff', fullfile(filepath, sprintf('tile_%02d.tif', i)));
    %         imp.close();
    %     end
    % 
    %     ImageJ
    %     ij.IJ.run("Grid/Collection stitching", "type=[Positions from file] " + ...
    %         "order=[Defined by TileConfiguration] " + ...
    %         "directory=[" + filepath + "] " + ...
    %         "layout_file=TileConfiguration.txt " + ...
    %         fusion_str + ...
    %         "regression_threshold=0.30 " + ...
    %         "max/avg_displacement_threshold=2.50 " + ...
    %         "absolute_displacement_threshold=3.50 " + ...
    %         "computation_parameters=[Save memory (but be slower)] " + ...
    %         "image_output=[Fuse and display]");
    % 
    % else
    % 
    %     ImageJ
    % 
    %     ij.IJ.run("Grid/Collection stitching", "type=[Positions from file] " + ...
    %         "order=[Defined by image metadata] " + ...
    %         "browse=" + ...
    %         file + ...
    %         " multi_series_file=" + ...
    %         file + ...
    %         fusion_str + ...
    %         "regression_threshold=0.30 " + ...
    %         "max/avg_displacement_threshold=2.50 " + ...
    %         "absolute_displacement_threshold=3.50 " + ...
    %         "compute_overlap increase_overlap=0 " + ...
    %         "invert_x " + ...
    %         "computation_parameters=[Save memory (but be slower)] " + ...
    %         "image_output=[Fuse and display]");
    % 
    % end
   
    [filepath,name,~] = fileparts(file);
    savef = filepath(1:end-4);
    savef = fullfile(savef,'stitch',name);
    if strcmp(fusionMethod,'Linear Blending')
        savef = strcat(savef,'_lin');
    end
    savef = strcat(savef,'.tif');

    ij.IJ.saveAs("Tiff", savef);
    ij.IJ.run("Close All");
    ij.IJ.run("Quit","");
    
    java.lang.Runtime.getRuntime.gc;

    junk = fullfile(filepath,'TileConfiguration.registered.txt');
    delete(junk)

    if shadingCorrection
        delete(file)
    end

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

end