function numChannels = getNumChannelsTif(fileTif)
file_info=imfinfo(fileTif);
if (file_info(1).SamplesPerPixel>1)
    numChannels = file_info(1).SamplesPerPixel; return
end

if isfield(file_info(1),"ImageDescription")
    a = split(file_info(1).ImageDescription,newline);
else
    numChannels=1;
    return
end
a = a(contains(a,'='));
a=split(a,'=');
[row,col] = find(strcmp(a,'channels'));
if (~isempty(row) && ~isempty(col))
    numChannels= str2double(a(row,col+1)); return
end

numChannels=1;