function process_histo(file)

    [~,~,ext] = fileparts(file);

    switch ext

        case '.nd2'
   
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
        
            ij.IJ.run("Z Project...", "projection=[Max Intensity]");
        
            [filepath,name,~] = fileparts(file);
            savef = filepath(1:end-4);
            savef = fullfile(savef,'processed',name);
            savef = strcat(savef,'.tif');
        
            ij.IJ.saveAs("Tiff", savef);
            ij.IJ.run("Close All");
            ij.IJ.run("Quit","");
            
            java.lang.Runtime.getRuntime.gc;
        
            junk = fullfile(filepath,'TileConfiguration.registered.txt');
            delete(junk)

        case '.vsi'
            
            java.lang.Runtime.getRuntime.gc;
            
            ImageJ

            ij.IJ.run("Bio-Formats", strcat("open=[",file) + ...
                "] color_mode=Default open_all_series " + ...
                "rois_import=[ROI manager] " + ...
                "view=Hyperstack " + ...
                "stack_order=XYCZT")

            ij.IJ.run("Close");
            ij.IJ.run("Close");
            ij.IJ.run("Close");
            ij.IJ.run("Close");
            ij.IJ.run("Close");
            ij.IJ.run("Close");

            ij.IJ.run("Make Composite");

            [filepath,name,~] = fileparts(file);
            savef = filepath(1:end-4);
            savef = fullfile(savef,'processed',name);
            savef = strcat(savef,'.tif');
        
            ij.IJ.saveAs("Tiff", savef);
            ij.IJ.run("Close All");
            ij.IJ.run("Quit","");
            
            java.lang.Runtime.getRuntime.gc;
           
    end

end