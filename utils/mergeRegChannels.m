% parse number of channels and number of planes in a tif

function outf = mergeRegChannels(fname, outFolder)

    listing = dir(outFolder);
    listing = {listing.name};
    out = regexp(listing,'channel_(\d*)','forceCellOutput');
    mask = zeros(1, length(out),'logical');

    for i=1:length(out)
        if out{i} == 1
            mask(i) = 1;
        end
    end

    chan_files = listing(mask);
    
    ImageJ

    for i=1:length(chan_files)
        ij.IJ.run("Bio-Formats", strcat("open=",fullfile(outFolder,chan_files{i})) + ...
        " color_mode=Default open_all_series " + ...
        "rois_import=[ROI manager] " + ...
        "view=Hyperstack " + ...
        "stack_order=XYCZT")
    end
    
    input_token = "";
    for i=1:length(chan_files)
        chan = sscanf(chan_files{i},'channel_%u.tif');
        formatSpec = 'c%u=channel_%u.tif ';
        str = sprintf(formatSpec,chan,chan);
        input_token = append(input_token,str);
    end

    input_token = append(input_token,"create");
    ij.IJ.run("Merge Channels...",input_token);
    
    [~,name,ext] = fileparts(fname);
    savef = fullfile(fileparts(outFolder),strcat(name,ext));
    ij.IJ.saveAs("Tiff", savef);
    ij.IJ.run("Close All");
    ij.IJ.run("Quit","");

    outf = savef;

end