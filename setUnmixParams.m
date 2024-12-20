function setUnmixParams(inputf, parentApp, mode)

parent_directory = pwd;
addpath(genpath(fullfile(parent_directory, 'utils')));

initFig(inputf, parentApp,mode);




function initFig(inputf, parentApp,mode)

h.mode = mode;
h.parentApp = parentApp;
[filepath,name,ext] = fileparts(inputf);
h.parent = filepath;
h.fn = strcat(name,ext);
autoload = 1;

bcol = [1 1 1];
pcol = [0.9 0.9 0.9];

h.fig = figure(427);
h.ax = axes();

set(h.fig, 'Units', 'pixels', 'Position', [100 100 1800 850], 'Color', bcol);



guidata(h.fig, h);
loadTifStack(h.fig)
cropStack(h.fig, autoload);
h = guidata(h.fig);

set(h.ax, 'Units', 'normalized', 'Position', [0.0272 0.0576 0.403 0.853]);
title(h.ax, h.fn);

h.vis = imagesc(h.ax, zeros(h.Nrow, h.Ncol, 3)); 
axis(h.ax, 'off');
axis(h.ax, 'image');

h.Nchan = size(h.stack, 3);

h.lipo.modelchans =  [1 4 5];
h.lipo.targetchans = [1 4 5];

h.lipo.bg = zeros(h.Nchan, 1);
h.lipo.amp = zeros(h.Nchan, 1);

possibleclrs = [...
    1 0 0;...
    0 1 0;...
    0 0 1;...
    1 0 1;...
    0 1 1;...
    1 1 0;...
    1 1 1];
defclrs = possibleclrs(1:h.Nchan, :);

h.str = {};

fn_split = split(h.fn,'__');
chans = fn_split{2};
h.str = split(chans,'_');

% in case additional image info was not included after a second double
% underscore
last_chan = split(h.str{end},'.');
h.str{end} = last_chan{1};

%Main panels
h.clim_pan = uipanel('Title','Channels', 'FontSize',12, 'Units', 'Normalized', ...
    'BackgroundColor',pcol, 'Position',[0.45 0.65 0.53 0.26]);

for i = 1:h.Nchans
    cmax = 255;
    yoff = 10;
    xoff = 130*(i-1);
    
    h.chanLabel_edit(i) = uicontrol('parent', h.clim_pan, 'Style', 'edit', 'String', h.str{i},'Units','pixels' ...
        ,'BackgroundColor', [1 1 1],'Position',[xoff+15 yoff+150 110 20],'HorizontalAlignment','Center', 'Fontsize', 7, 'Fontweight', 'Bold');
    
    h.color_button(i) = uicontrol('parent', h.clim_pan,'Style', 'pushbutton',  'Units', 'Pixels', 'Position', ...
    [xoff+15, yoff+115, 110, 30], 'String', '', 'Callback', {@colorButtonFcn,gcf, i}, 'BackgroundColor', defclrs(i, :));
    
    uicontrol('parent', h.clim_pan, 'Style', 'text', 'String', 'Min: ','Units','pixels' ...
        ,'BackgroundColor',pcol,'Position',[xoff+14 yoff+66 40 24],'HorizontalAlignment','Right');
    uicontrol('parent', h.clim_pan, 'Style', 'text', 'String', 'Max: ','Units','pixels' ...
        ,'BackgroundColor',pcol,'Position',[xoff+14 yoff+37 40 24],'HorizontalAlignment','Right');
    
    h.cmin_edit(i) = uicontrol('parent', h.clim_pan, 'Style', 'edit', 'String', '0','Units','pixels' ...
        ,'BackgroundColor',[1 1 1],'Position',[xoff+61 yoff+72 56 24],'HorizontalAlignment','Right', ...
        'Callback', {@updateImage,gcf}, 'Enable', 'Off');
    
    h.cmax_edit(i) = uicontrol('parent', h.clim_pan, 'Style', 'edit', 'String', num2str(cmax),'Units','pixels' ...
        ,'BackgroundColor',[1 1 1],'Position',[xoff+61 yoff+42 56 24],'HorizontalAlignment','Right', ...
        'Callback', {@updateImage,gcf}, 'Enable', 'Off');
    
    h.autolim(i) = uicontrol('parent', h.clim_pan, 'Style', 'radiobutton',  'Units', 'pixels', 'Position', ...
        [xoff+78 yoff+7 70 24], 'String', 'Auto','FontWeight','Bold', 'Min', 0, 'Max', 1, ...
        'BackgroundColor', pcol, 'Value', 1, 'Callback', {@enableLims,gcf,i});
    
    h.enable_checkbox(i) = uicontrol('parent', h.clim_pan, 'Style', 'checkbox', 'Value', i<4,'Units','pixels' ...
    ,'BackgroundColor',pcol,'Position',[xoff+28 yoff+7 46 24],'HorizontalAlignment','Right', ...
    'Callback', {@updateImage,gcf}, 'BackgroundColor', pcol, 'String', 'On', 'FontSize', 10);

    uicontrol('parent', h.clim_pan, 'Style', 'text', 'String', ' ','Units','pixels' ...
        ,'BackgroundColor',pcol-0.1,'Position',[xoff+132 yoff+10 5 150],'HorizontalAlignment','Right');
end

h.hilo_checkbox = uicontrol('parent', h.clim_pan, 'Style', 'checkbox', 'Value', 0,'Units','pixels' ...
    ,'BackgroundColor',pcol,'Position',[800 150 200 24],'HorizontalAlignment','Right', ...
    'Callback', {@updateImage,gcf}, 'BackgroundColor', pcol, 'String', 'HiLo LUT', 'FontSize', 10);

h.save_button = uicontrol('Style', 'pushbutton',  'Units', 'Normalized', 'Position', ...
    [0.82 0.14 0.075 0.05], 'String', 'Save settings','FontWeight','Bold', 'Callback', ...
    {@saveButtonFcn,gcf});


uicontrol('Style', 'text', 'Units', 'Normalized', 'String', 'Spectral unmixing', ...
    'BackgroundColor',bcol-0.1,'Position',[0.45 0.61 0.225, 0.025],'HorizontalAlignment','Center', ...
    'FontSize', 12, 'FontWeight', 'Bold');



data = cell(10, 5);
data(:,:) = {0};
data(1,:) = {1, 2, 0, 1, false};
data(:, 5) = {false};

h.unmixTable = uitable(h.fig, 'Units', 'Normalized', 'Position', [0.45 0.475 0.225 0.135], ...
    'ColumnWidth', {80, 80, 70, 70, 70}, 'ColumnName', {'Unmix from', 'Unmix to', ...
    'Offset', 'Scale', 'Enable'}, 'Data', data, 'ColumnEditable', [true, true, true, true, true]);



uicontrol('Style', 'text', 'Units', 'Normalized', 'String', 'Image shifts', ...
    'BackgroundColor',bcol-0.1,'Position',[0.68 0.61 0.135, 0.025],'HorizontalAlignment','Center', ...
    'FontSize', 12, 'FontWeight', 'Bold');

data = cell(5, 3);
data(1,:) = {1, 0, 0};

h.shiftTable = uitable(h.fig, 'Units', 'Normalized', 'Position', [0.68 0.475 0.135 0.135], ...
    'ColumnWidth', {70, 70, 70}, 'ColumnName', {'Channel', 'X shift', ...
    'Yshift'}, 'Data', data, 'ColumnEditable', [true, true, true]);


h.opts_pan = uipanel('Title','Options', 'FontSize',12, 'Units', 'Normalized', ...
    'BackgroundColor',pcol, 'Position',[0.82 0.26 0.16 0.375]);

xx = -790;
yy = -50;
h.lipo_checkbox = uicontrol('parent', h.opts_pan,'Style', 'checkbox', 'Value', 0,'Units','pixels' ...
    ,'BackgroundColor',[1 1 1],'Position',[800+xx 120+yy 150 20],'HorizontalAlignment','Right', ...
    'FontWeight', 'Bold', 'BackgroundColor', pcol, 'String', 'Remove lipofuscin');

h.crosshairs_checkbox = uicontrol('parent', h.opts_pan,'Style', 'checkbox', 'Value', 0,'Units','pixels' ...
    ,'BackgroundColor',[1 1 1],'Position',[950+xx 120+yy 150 20],'HorizontalAlignment','Right', ...
    'FontWeight', 'Bold', 'BackgroundColor', pcol, 'String', 'Show ROIs', ...
    'Callback', {@lipoCheckboxCallback,gcf});

h.lipoAx = axes('parent', h.opts_pan);
set(h.lipoAx, 'Units', 'pixels', 'Position', [820+xx 160+yy 220 85]);
plot(h.lipoAx, h.lipo.bg); hold on; plot(h.lipoAx, h.lipo.amp);

h.lipomodelchans_edit = uicontrol('parent', h.opts_pan,'Style', 'edit', 'String', num2str(h.lipo.modelchans),'Units','pixels' ...
    ,'BackgroundColor',[1 1 1],'Position',[800+xx 315+yy 70 20],'HorizontalAlignment','Right');
uicontrol('parent', h.opts_pan,'Style', 'text', 'String', 'Channels to build model from','Units','pixels' ...
    ,'BackgroundColor',pcol,'Position',[880+xx 315+yy 150 17],'HorizontalAlignment','Left');

h.liporemovechans_edit = uicontrol('parent', h.opts_pan,'Style', 'edit', 'String', num2str(h.lipo.targetchans),'Units','pixels' ...
    ,'BackgroundColor',[1 1 1],'Position',[800+xx 290+yy 70 20],'HorizontalAlignment','Right');
uicontrol('parent', h.opts_pan,'Style', 'text', 'String', 'Channels to remove lipo from','Units','pixels' ...
    ,'BackgroundColor',pcol,'Position',[880+xx 290+yy 150 17],'HorizontalAlignment','Left');

h.updateLipoModel_button = uicontrol('parent', h.opts_pan,'Style', 'pushbutton',  'Units', 'pixels', 'Position', ...
    [800+xx 255+yy 260 25], 'String', 'Update lipo model','FontWeight','Bold', 'Callback', ...
    {@updateLipoModel,gcf});



xx = -790;
yy = -60;
h.LPF_checkbox = uicontrol('parent', h.opts_pan,'Style', 'checkbox', 'Value', 0,'Units','pixels' ...
    ,'BackgroundColor',[1 1 1],'Position',[800+xx 100+yy 120 20],'HorizontalAlignment','Right', ...
    'FontWeight', 'Bold', 'BackgroundColor', pcol, 'String', 'Spatial Filter');

h.radius_edit = uicontrol('parent', h.opts_pan,'Style', 'edit', 'String', '1.25','Units','pixels' ...
    ,'BackgroundColor',[1 1 1],'Position',[800+xx 75+yy 50 20],'HorizontalAlignment','Right');

uicontrol('parent', h.opts_pan,'Style', 'text', 'String', 'width (pixels)','Units','pixels' ...
    ,'BackgroundColor',pcol,'Position',[860+xx 75+yy 100 17],'HorizontalAlignment','Left');


h.updateImage_button = uicontrol('Style', 'pushbutton',  'Units', 'Normalized', 'Position', ...
    [0.82 0.2 0.16 0.05], 'String', 'Update Image','FontWeight','Bold', 'Callback', ...
    {@getImageStack,gcf});

h.batchProcess_button = uicontrol('Style', 'pushbutton',  'Units', 'Normalized', 'Position', ...
    [0.82 0.02 0.16 0.05], 'String', 'Batch process','FontWeight','Bold', 'Callback', ...
    {@batchProcess,gcf});

h.bgROI = images.roi.Rectangle(h.ax,'Position',[100 100 25 25], 'LineWidth', 3, 'Visible', 'Off');
h.lipoCrosshair = drawcrosshair('Parent', h.ax, 'Position', [50, 50], 'Visible', 'Off');


imFns = findImageFiles(h.parent, '.tif');

uicontrol('Style', 'text', 'Units', 'Normalized', 'String', 'Image files', ...
    'BackgroundColor',bcol-0.1,'Position',[0.45 0.42 0.365 0.025],'HorizontalAlignment','Center', ...
    'FontSize', 12, 'FontWeight', 'Bold');
h.file_listbox = uicontrol('Style', 'listbox', 'Units', 'normalized', 'Position', ...
    [0.45 0.14 0.365 0.28], 'String', imFns, 'Value', 1, 'BackgroundColor', [1 1 1], ...
    'Max', 3, 'Min', 1);

h.switchFile_button = uicontrol('Style', 'pushbutton',  'Units', 'Normalized', 'Position', ...
    [0.45 0.08 0.16 0.05], 'String', 'Switch file','FontWeight','Bold', 'Callback', ...
    {@switchFile,gcf});

h.makeComposite_button = uicontrol('Style', 'pushbutton',  'Units', 'Normalized', 'Position', ...
    [0.63 0.08 0.16 0.05], 'String', 'Composite figure','FontWeight','Bold', 'Callback', ...
    {@makeComposite,gcf});

guidata(h.fig, h);
updateImage([], [], h.fig);

if strcmp(h.mode,'write')
    loadButtonFcn([],[],gcf)
    getImageStack([],[],gcf)
    exportImageStack([],[],gcf)
    h.parentApp.finish_unmix()
    close(h.fig)
end





function makeComposite(~, ~, fig)

h = guidata(fig);

Ncol = 3;
Nrow = ceil((h.Nchans+1)./Ncol);

dr = 1./Nrow;
dc = 1./Ncol;

figure(643);
ax(1) = axes;
set(ax(1), 'Units', 'Normalized', 'Position', [0 (1-dr) dc dr])
data = get(h.vis, 'CData');
imagesc(data);
axis square off;




ax = zeros(h.Nchans, 1);
for i = 1:h.Nchans
    loc = i+1;
    row = ceil(loc/Ncol);
    col = mod(loc-1, Ncol)+1;

    ax(i) = axes;
    set(ax(i), 'Units', 'Normalized', 'Position', [dc*(col-1) (1-dr)-(row-1)*dr dc dr]);

    data = zeros(h.Nrow, h.Ncol, 3);
    clr = get(h.color_button(i), 'BackgroundColor');
    scim = getScaledImage(h, i);

    for j = 1:3
        data(:,:,j) = data(:,:,j)+scim.*clr(j);
    end

    imagesc(data);
    axis square off;

end
    


function switchFile(~, ~, fig)
h = guidata(fig);
filelist = get(h.file_listbox, 'String');
num = get(h.file_listbox, 'Value');
if numel(num)>1
    num = num(1);
end

h.fn = filelist{num};
guidata(fig, h);

loadTifStack(fig);
cropStack(fig, 0);
getImageStack([], [], fig)






function imFns = findImageFiles(parent, token)

contents = dir(parent);
good = false(numel(contents), 1);
for i = 1:numel(good)
   if ~isempty(strfind(contents(i).name, token))
       good(i) = true;
   end
end

contents = contents(good);
imFns = cell(numel(contents), 1);
for i = 1:numel(contents)
    imFns{i} = contents(i).name;
end






function loadTifStack(fig)

h = guidata(fig);
fullfn = fullfile(h.parent, h.fn);
h.info = imfinfo(fullfn);   

h.stack = imread(fullfn);
h.Nrow = h.info(1).Height;
h.Ncol = h.info(1).Width;
h.ymax = h.Nrow;
h.xmax = h.Ncol;

[h.Nchans, h.Nplanes] = parseInfo(h.info(1));
h.stack = tiffreadVolume(fullfn);
% h.stack = reshape(h.stack, h.Nrow, h.Ncol, h.Nchans, h.Nplanes);
h.stack = permute(h.stack, [1 2 4 3]);

% cnt = 0;
% for i = 1:h.Nplanes
%     for j = 1:h.Nchans 
%         cnt = cnt+1;
%         h.stack(:,:,j,i) = imread(fullfn, 'tif', 'Info', h.info, 'Index', cnt);
%     end
% end

h.tiffstack = h.stack;
guidata(fig, h);




function [chans, planes] = parseInfo(info)
id = info.ImageDescription;

tok = regexp(id, 'channels=(\d*)', 'tokens');
if isempty(tok)
    disp('Could not read number of channels in tif.  Assuming 1 channel');
    chans = 1;
else
    chans = str2double(tok{1}{1});
end

tok = regexp(id, 'slices=(\d*)', 'tokens');

if isempty(tok)
    tok = regexp(id, 'frames=(\d*)', 'tokens');
end
if isempty(tok)
    disp('Could not read number of z planes in tif.  Assuming 1 plane');
    planes = 1;
else
    planes = str2double(tok{1}{1});
end




function cropStack(fig, autoload)

h = guidata(fig);

if autoload
    maxpix = inf;
    h.plane = ceil(h.Nplanes/2);
else
    defplane = ceil(h.Nplanes/2);
    prompt = {['Enter plane to load (of ' num2str(h.Nplanes) '):'],'Image size (''inf'' for maximum)'};
    dlgtitle = 'Input';
    dims = [1 35];

    definput = {num2str(defplane), 'inf'};

    answer = inputdlg(prompt,dlgtitle,dims,definput);
    h.plane = str2double(answer{1});

    if strcmpi(answer{2}, 'inf')
        maxpix = inf;
    else
        maxpix = str2double(answer{2});
    end
end

h.xmax = min(maxpix, h.Ncol);
h.ymax = min(maxpix, h.Nrow);
h.stack = h.stack(1:h.ymax, 1:h.xmax, :, h.plane);


guidata(fig, h);



function loadButtonFcn(~, ~, fig)

h = guidata(fig);

[filepath,~,~] = fileparts(h.parent);
filepath = fullfile(filepath,'internal','uparams');
[~,name,~] = fileparts(h.fn);
loadf = fullfile(filepath,strcat(name,'.mat'));
load(loadf);

h.str = params.str;
for i = 1:h.Nchans
    set(h.chanLabel_edit(i), 'String', h.str{i});
end

for i = 1:h.Nchan
    set(h.color_button(i), 'BackgroundColor', params.clr{i});
    set(h.cmin_edit(i), 'String', params.cmin{i});
    set(h.cmax_edit(i), 'String', params.cmax{i});
    set(h.autolim(i), 'Value', params.autolim(i));
    set(h.enable_checkbox(i), 'Value', params.enable(i));
    enableLims([], [], h.fig, i);
end

h.shiftTable.Data = params.shiftDat;
h.unmixTable.Data = params.unmixDat;

h.lipo = params.lipo;
set(h.lipo_checkbox, 'Value', params.lipocheck);

set(h.lipomodelchans_edit, 'String', num2str(params.lipo.modelchans));
set(h.liporemovechans_edit, 'String', num2str(params.lipo.targetchans));

hold (h.lipoAx, 'Off');
plot(h.lipoAx, h.lipo.bg);
hold (h.lipoAx, 'On');
plot(h.lipoAx, h.lipo.amp);


set(h.LPF_checkbox, 'Value', params.lpf);
set(h.radius_edit, 'String', params.radius);


guidata(fig, h);
getImageStack([], [], h.fig);
updateImage([], [], h.fig);



function saveButtonFcn(~, ~, fig)

h = guidata(fig);

[filepath,~,~] = fileparts(h.parent);
filepath = fullfile(filepath,'internal','uparams');
[~,name,~] = fileparts(h.fn);
savef = fullfile(filepath,strcat(name,'.mat'));

params.str = cell(h.Nchans, 1);
for i = 1:h.Nchans
    params.str{i} = get(h.chanLabel_edit(i), 'String');
end

for i = 1:h.Nchans
    params.clr{i}  = get(h.color_button(i), 'BackgroundColor');
    params.cmin{i} = get(h.cmin_edit(i), 'String');
    params.cmax{i} = get(h.cmax_edit(i), 'String');
    params.autolim(i) = get(h.autolim(i), 'Value');
    params.enable(i) = get(h.enable_checkbox(i), 'Value');
end

params.unmixDat = h.unmixTable.Data;
params.shiftDat = h.shiftTable.Data;

params.lpf = get(h.LPF_checkbox, 'Value');
params.radius = get(h.radius_edit, 'String');

params.lipo           = h.lipo;
params.lipocheck = get(h.lipo_checkbox, 'Value');

save(savef, 'params');
h.parentApp.finish_unmix_params()
close(fig)




function batchProcess(~, ~, fig)
h = guidata(fig);

% defpath = fileparts(h.info(1).Filename);
% [fns, h.parent] = uigetfile('*.tif', 'Select files to process', defpath, 'MultiSelect', 'On');
% if ~iscell(h.fn)
%     if numel(h.fn)==1
%         return;
%     else
%         h.fn = {h.fn};
%     end
% end

fns = get(h.file_listbox, 'String');
nums = get(h.file_listbox, 'Value');


for i = 1:numel(nums)
    h.fn = fns{nums(i)};
    guidata(fig, h);
    loadTifStack(fig);
    exportImageStack([], [], fig);
end
disp('Done Batch Processing!!!!')




function exportImageStack(~, ~, fig)
h = guidata(fig);

[filepath,~,~] = fileparts(h.parent);
filepath = fullfile(filepath,'unmix');
[~,name,~] = fileparts(h.fn);
savef = fullfile(filepath,strcat(name,'.tif'));

t = Tiff(savef, 'w8');
tagstruct.ImageLength = h.Nrow;
tagstruct.ImageWidth = h.Ncol;
tagstruct.SampleFormat = Tiff.SampleFormat.UInt; % uint
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SamplesPerPixel = h.Nchans;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
for ii=1:h.Nplanes
    h.plane = ii;
    guidata(fig, h);
    getImageStack([], [], fig)
    h = guidata(fig);
    im = h.stack;
    setTag(t,tagstruct);
    warning('off','MATLAB:imagesci:Tiff:missingExtraSamples')
    write(t,uint16(im));
    writeDirectory(t);
end
close(t)

% tmp.plane = h.plane;
% tmp.ymax = h.ymax;
% tmp.xmax = h.xmax;
% 
% h.ymax = h.Nrow;
% h.xmax = h.Ncol;
% 
% cnt = 0;
% for i = 1:h.Nplanes
%     h.plane = i;
% 
%     guidata(fig, h);
%     getImageStack([], [], fig)
% 
%     h = guidata(fig);
%     im = h.stack;
% 
%     for j = 1:h.Nchans
%         cnt = cnt+1;
% 
%         if cnt == 1
%             imwrite(uint16(squeeze(im(:,:,j))), fn, 'tif', 'writemode', 'overwrite', ...
%                 'Compression', 'None', 'Description', h.info(1).ImageDescription);
%         else
%             imwrite(uint16(squeeze(im(:,:,j))), fn, 'tif','writemode', 'append', ...
%                 'Compression', 'None', 'Description', h.info(1).ImageDescription);
%         end
%     end
%     updateImage([], [], fig);
% end



% h.plane = tmp.plane;
% h.ymax = tmp.ymax;
% h.xmax = tmp.xmax;

guidata(fig, h);
getImageStack([],[], fig)




function getImageStack(~, ~, fig)

refreshStack(fig);
doShifts(fig);
lowPassFilterStack(fig)
doUnmixing(fig);
removeLipo(fig);

updateImage([], [], fig);





function refreshStack(fig)

h = guidata(fig);
h.stack = h.tiffstack(1:h.ymax, 1:h.xmax, :, h.plane);
guidata(fig, h);




function doShifts(fig)
h = guidata(fig);

newstack = zeros(h.ymax, h.xmax, h.Nchans);

xshift = zeros(h.Nchans, 1);
yshift = zeros(h.Nchans, 1);

shiftDat = h.shiftTable.Data;
for i = 1:size(shiftDat, 1)
    if ~isempty(shiftDat{i, 1})
        xshift(shiftDat{i, 1}) = shiftDat{i, 2};
        yshift(shiftDat{i, 1}) = shiftDat{i, 3};
    end
end

for i = 1:h.Nchans 
    
    y1 = 1:h.ymax;
    y2 = 1:h.ymax;
    if yshift(i)>0
        y1 = 1+yshift(i):h.ymax;
        y2 = 1:h.ymax-yshift(i);
%         disp('Shifting up');
    elseif yshift(i)<0
        y1 = 1:h.ymax+yshift(i);
        y2 = 1-yshift(i):h.ymax;
%         disp('Shifting down');
    end
    
    x1 = 1:h.xmax;
    x2 = 1:h.xmax;
    if xshift(i)>0
        x1 = 1+xshift(i):h.xmax;
        x2 = 1:h.xmax-xshift(i);
%         disp('Shifting right');
    elseif xshift(i)<0
        x1 = 1:h.xmax+xshift(i);
        x2 = 1-xshift(i):h.xmax;
%         disp('Shifting left');
    end
    
    newstack(y1,x1,i) = h.stack(y2, x2, i);
end
h.stack = newstack;
guidata(fig, h);



function lowPassFilterStack(fig)

h = guidata(fig);
if get(h.LPF_checkbox, 'Value')
    pix = str2double(get(h.radius_edit, 'String'));
    h.stack = gaussFilt(h.stack, pix);
end
guidata(fig, h);



function doUnmixing(fig)
h = guidata(fig);

% 5x5 cell array representing values in spectral unmixing panel
tabDat = h.unmixTable.Data;

% Y by X by channels
newstack = h.stack;

% iterate through rows
for i = 1:size(tabDat, 1)

    % if spectral unmixing row is enabled
    if ~isempty(tabDat{i, 5}) && tabDat{i, 5}
        
        from = tabDat{i, 1};
        to   = tabDat{i, 2};
        
        % sourceIm = image to be subtracted from image where unmixing is
        % desired
        sourceIm = double(newstack(:,:,from));
        offset = tabDat{i, 3};
        scale = tabDat{i, 4};

        sourceIm = (sourceIm - offset).*scale;
        sourceIm(sourceIm<0) = 0;

        unmixedIm = newstack(:,:,to) - sourceIm;
        unmixedIm(unmixedIm < 0) = 0;

        newstack(:,:,to) = unmixedIm;
    end
end
h.stack = newstack;
guidata(fig, h);





function removeLipo(fig)
h = guidata(fig);
newstack = h.stack;

if get(h.lipo_checkbox, 'Value')


    h.lipo.modelchans = str2num(get(h.lipomodelchans_edit, 'String'));
    h.lipo.targetchans = str2num(get(h.liporemovechans_edit, 'String'));

    lim = zeros(size(h.stack, 1), size(h.stack, 2), numel(h.lipo.modelchans));
    for i = 1:numel(h.lipo.modelchans)
        amp = h.lipo.amp(h.lipo.modelchans(i));
        bg = h.lipo.bg(h.lipo.modelchans(i));

        lim(:,:,i) = (h.stack(:,:,h.lipo.modelchans(i)) - bg)./amp;  %lipo has amplitude of 1
    end
    bgim = min(lim, [], 3);
    bgim(bgim<0) = 0;
    
    for i = 1:numel(h.lipo.targetchans)
        chan = h.lipo.targetchans(i);
        
        bg = bgim.*h.lipo.amp(chan);
        newstack(:,:,chan) = h.stack(:,:,chan) - bg;
    end
end
h.stack = newstack;
guidata(fig, h)























































