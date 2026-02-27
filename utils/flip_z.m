% flip the z direction of a 3-D image (tif) using imageJ, expects a full 
% file path, will overwrite the original file with the same name
function flip_z(file)

    java.lang.Runtime.getRuntime.gc;
    
    ImageJ

    ij.IJ.run("Bio-Formats", strcat("open=",file) + ...
        " color_mode=Default open_all_series " + ...
        "rois_import=[ROI manager] " + ...
        "view=Hyperstack " + ...
        "stack_order=XYCZT")

    ij.IJ.run("Flip Z");

    ij.IJ.saveAs("Tiff", file);
    ij.IJ.run("Close All");
    ij.IJ.run("Quit","");
    
    java.lang.Runtime.getRuntime.gc;

end