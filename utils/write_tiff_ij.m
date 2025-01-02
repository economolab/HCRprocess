% parse number of channels and number of planes in a tif

function write_tiff_ij(fullfn)
    
    java.lang.Runtime.getRuntime.gc;
    
    ImageJ

    ij.IJ.run("Bio-Formats", strcat("open=",fullfn) + ...
        " color_mode=Default open_all_series " + ...
        "rois_import=[ROI manager] " + ...
        "view=Hyperstack " + ...
        "stack_order=XYCZT")

    ij.IJ.saveAs("Tiff", fullfn);
    ij.IJ.run("Close All");
    ij.IJ.run("Quit","");

    

    java.lang.Runtime.getRuntime.gc;
    
end