function process_final(file, radius)
   
    [filepath,name,ext] = fileparts(file);
    orig_file = strcat(name,ext);

    newStr = split(name,"__");
    slice = newStr{1};
    channels = newStr{2};
    channels = split(channels,"_");

    [filepath,~,~] = fileparts(filepath);
    [~,name,~] = fileparts(filepath);
    newStr = split(name,"_");
    
    for i=1:length(newStr)
        if contains(newStr{i},'r')
            round = newStr{i};
        end
    end

    for i=1:length(channels)
        channels{i} = strcat(round,'_',channels{i},'.tif');
    end
   
    [filepath,~,~] = fileparts(filepath);
    savepath = fullfile(filepath,'post','core_output',slice);

    java.lang.Runtime.getRuntime.gc;

    ImageJ

    ij.IJ.run("Bio-Formats", strcat("open=",file) + ...
        " color_mode=Default open_all_series " + ...
        "rois_import=[ROI manager] " + ...
        "view=Hyperstack " + ...
        "stack_order=XYCZT")
    
    spec = strcat("rolling=",radius);
    ij.IJ.run("Subtract Background...", spec);
    ij.IJ.run("Split Channels");
    
    for i=length(channels):-1:1
        ij.IJ.resetMinAndMax()
        ij.IJ.saveAs("Tiff",fullfile(savepath,channels{i}))
        ij.IJ.run("Close");
    end

    ij.IJ.run("Quit","");

    java.lang.Runtime.getRuntime.gc;

end