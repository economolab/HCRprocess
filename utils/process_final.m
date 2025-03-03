function process_final(file, subBG, radius)
   
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

    if subBG == 1
        ij.IJ.run("Subtract Background...", spec);
    end
    
    ij.IJ.run("Split Channels");
    pause(1)
    
    for i=length(channels):-1:1
        ij.IJ.resetMinAndMax()
        pause(1)
        ij.IJ.saveAs("Tiff",fullfile(savepath,channels{i}))
        pause(1)
        ij.IJ.run("Close");
        pause(1)
    end

    ij.IJ.run("Quit","");
    pause(10)

    java.lang.Runtime.getRuntime.gc;

end