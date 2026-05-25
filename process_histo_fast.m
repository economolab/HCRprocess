[filepath,~,~] = fileparts(mfilename('fullpath'));
cd(filepath)
addpath(['./' 'utils'])
folderName = params('fijiScriptsDir');
addpath(folderName)

title = 'Select directory containing histology images';
selpath = uigetdir(pwd,title);

if ~exist(fullfile(selpath,'raw'),'dir')
    
    mkdir(fullfile(selpath,'raw'))
    
    listing = fetch_dir(selpath,'file');
    
    if ~isempty(listing)
        for i=1:length(listing)
            movefile(fullfile(selpath,listing{i}), fullfile(selpath,'raw'))
        end
    end
    
    listing = fetch_dir(fullfile(selpath),'dir');
    listing(ismember(listing,'raw')) = [];
    
    if ~isempty(listing)
        for i=1:length(listing)
            movefile(fullfile(selpath,listing{i}), fullfile(selpath,'raw'))
        end
    end

end

check_build_dir(fullfile(selpath,'processed'))

listing = fetch_dir(fullfile(selpath,'raw'), 'file');

for i=1:length(listing)
    process_histo(fullfile(selpath,'raw',listing{i}))
end


