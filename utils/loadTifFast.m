function im = loadTifFast(FileTif,targetChannel)
 if ~exist('targetChannel','var')
     % third parameter does not exist, so default it to something
      targetChannel = 1;
 end
disp(['Loading ' FileTif '...']);
warning off;

numChannels = getNumChannelsTif(FileTif);

if (targetChannel>numChannels)
    disp(['only detected ' num2str(numChannels) ' channels, but target channel was ' num2str(targetChannel) '. Setting target to 1.']);
    targetChannel=1;
end

im = tiffreadVolume(FileTif,"PixelRegion",{[1 1 inf],[1 1 inf],[targetChannel numChannels inf]});
% im = tiffreadVolume(FileTif);

warning on;
disp('Done loading Tif');