function numChannels = getNumChannelsTifDirectory(directory)
files = dir(directory);
for i = 1:length(files)
    [~, ~, fExt] = fileparts(files(i).name); 
    if strcmp(fExt,'.tif')
        numChannels = getNumChannelsTif(fullfile(files(i).folder,files(i).name));
        return
    end
end
numChannels=1;