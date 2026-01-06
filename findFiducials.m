function findFiducials(fixed_fn, moving_fn, fixed_reg_channel, moving_reg_channel, MainApp, load_fid)

parent_directory = pwd;
 
fn = {};
im = loadData(fn, parent_directory);
addpath(genpath(fullfile(parent_directory, 'reg_utils')));
im.fn = fn;

% fixed_fn = 'F:\proto_dir\2024-12-15_r0\stitch\s02__NT445_GFP488_dTom547.tif';
% moving_fn = 'F:\proto_dir\2024-12-15_r0\stitch\s02__NT445_GFP488_dTom547.tif';
% fixed_reg_channel = 1;
% moving_reg_channel = 1;

im.fixed_fn = fixed_fn;
im.moving_fn = moving_fn;
im.fixed_reg_channel = fixed_reg_channel;
im.moving_reg_channel = moving_reg_channel;
im.MainApp = MainApp;

initGUI(im, load_fid);




function initGUI(h, load_fid)



h.fig = figure(124);
h.fig.WindowState = 'maximized';
BackCol = [1 1 1];
set(h.fig, 'CloseRequestFcn', @close_req);
set(h.fig, 'Units', 'Normalized', 'Position', [0.01 0.1 0.87 0.7]);
set(h.fig, 'DefaultAxesLineWidth', 2, 'DefaultAxesFontSize', 12, 'Color', BackCol);
set(h.fig, 'WindowScrollWheelFcn', {@figScroll, h.fig});
set(h.fig, 'WindowButtonDownFcn', {@figClick, h.fig});
set(h.fig, 'Renderer', 'openGL');

fig = uifigure;
d = uiprogressdlg(fig,'Title','Loading...', ...
    'Message',strcat("Loading registration images and fiducials for: ",h.moving_fn), ...
    'Indeterminate','on');

h.ax(1) = axes;
h.axim(1) = imagesc(h.fixed(:,:,1)); hold on;
colormap(gray);
a = prctile(h.fixed(1:(size(h.fixed,1)*size(h.fixed,2))), [0.01, 99.9]);

cmin(1) = a(1);
cmax(1) = a(2);

if (cmax(1)>cmin(1))
    clim([cmin(1) cmax(1)]);
end
axis off;
axis image;

h.ax(2) = axes;  
h.axim(2) = imagesc(h.moving(:,:,1));  hold on;
colormap(gray);
a = prctile(h.moving(1:(size(h.moving,1)*size(h.moving,2))), [0.01, 99.9]);
cmin(2) = a(1);
cmax(2) = a(2);

clim([cmin(2) cmax(2)]);
axis off;
axis image;

h.ax(3) = axes;
h.axim(3) = imagesc(squeeze(h.fixed(:,1,:))); hold on;
colormap(gray);
if (cmax(1)>cmin(1))
    clim([cmin(1) cmax(1)]);
end
axis off;
axis image;

h.ax(4) = axes;
h.axim(4) = imagesc(squeeze(h.moving(:,1,:))); hold on;
colormap(gray);
clim([cmin(2) cmax(2)]);
axis off;
axis image;

num = 1;
h.point(num) = drawpoint(h.ax(num), 'Color', 'g', 'Position', [round(size(h.fixed, 2)./2), round(size(h.fixed, 1)./2)]);
% setColor(h.point(num), 'g');
num = 2;
h.point(num) = drawpoint(h.ax(num), 'Color', 'm', 'Position', [round(size(h.moving, 2)./2), round(size(h.moving, 1)./2)]);
% setColor(h.point(num), 'm');
h.point(3) = drawpoint(h.ax(num),'Position', [round(size(h.moving, 2)./2), round(size(h.moving, 1)./2)]);
delete(h.point(3));

addlistener(h.point(1),"MovingROI",@(src,pos) movebothpoints(gcf, 1));
addlistener(h.point(2),"MovingROI",@(src,pos) movebothpoints(gcf, 2));

h.plus(1) = plot(h.ax(3), round(size(h.fixed, 3)./2), round(size(h.fixed, 1)./2), 'g+');
h.plus(2) = plot(h.ax(4), round(size(h.moving, 3)./2), round(size(h.moving, 1)./2), 'm+');

set(h.ax(1), 'Position', [0.005 0.375 0.44 0.625]); 
set(h.ax(2), 'Position', [0.45 0.375 0.44 0.625]);
set(h.ax(3), 'Position', [0.07 0.05 0.3 0.25]);
set(h.ax(4), 'Position', [0.52 0.05 0.3 0.25]);

h.temppts{1} = [];
h.temppts{2} = [];

h.tree.mov = [];
h.tree.cur = [];
h.treepts = [];

h.max_depth = 5.51;
% xx = -10*h.max_depth:1:10*h.max_depth;
% h.treecmap = hsv2rgb([0.5.*ones(numel(xx), 1), ones(numel(xx), 1), 1 - abs(xx')./(15.*max_depth)]);


h.pts.fix = [];
h.pts.mov = [];
h.pts.cur = [];

h.haveaf = 0;
h.havebs = 0;
h.haveont = isfield(h, 'ont')&isfield(h, 'an');
h.doAfWarpfast = @exposedDoAfWarpfast;


pcol = [0.9 0.9 0.9];
h.panfix = uipanel('Title','Fixed','FontSize',12, 'Units', 'Normalized', ...
    'BackgroundColor',pcol, 'Position',[0.005 0.3 0.44 0.07]);
h.panmov = uipanel('Title','Moving','FontSize',12, 'Units', 'Normalized', ...
    'BackgroundColor',pcol, 'Position',[0.45 0.3 0.44 0.07]);
h.panaxes = uipanel('Title','Axes control','FontSize',12, 'Units', 'Normalized', ...
    'BackgroundColor',[0.9 0.9 0.9], 'Position',[0.82 0.13 0.07 0.165]);
h.panauto = uipanel('Title','AutoReg','FontSize',12, 'Units', 'Normalized', ...
             'BackgroundColor',pcol, 'Position',[0.82 0.01 0.07 0.11]);
h.clim1_pan = uipanel('Title','B/C','FontSize',12, 'Units', 'Normalized', ...
    'BackgroundColor',pcol, 'Position',[0.005 0.16 0.08 0.125]);
h.clim2_pan = uipanel('Title','B/C','FontSize',12, 'Units', 'Normalized', ...
    'BackgroundColor',pcol, 'Position',[0.45 0.16 0.08 0.125]);

frames = size(h.fixed,3);
h.frame_slider(1) = uicontrol('Parent', h.panfix, 'Style', 'Slider', 'Units', 'normalized', 'Position', ...
    [0.13 0.15 0.85 0.3], 'Min', 1, 'Max', frames, 'Value', 1, 'SliderStep', 1./(max(frames-1, 1))+[0 0], ...
    'Callback', {@updateImage,gcf}, 'BackgroundColor', [0.75 0.75 0.75]);
uicontrol('parent', h.panfix, 'Style', 'text', 'String', 'Frame: ','Units','normalized' ...
    ,'BackgroundColor',pcol,'Position',[0.01 0.15 0.05 0.3],'HorizontalAlignment','Right', ...
    'FontSize', 10);
h.framenumtext(1) = uicontrol('parent', h.panfix, 'Style', 'text', 'String', '1/1 ','Units','normalized' ...
    ,'BackgroundColor',pcol,'Position',[0.065 0.15 0.06 0.3],'HorizontalAlignment','Left', ...
    'FontSize', 10);

frames = size(h.moving,3);
h.frame_slider(2) = uicontrol('Parent', h.panmov, 'Style', 'Slider', 'Units', 'normalized', 'Position', ...
    [0.13 0.15 0.85 0.3], 'Min', 1, 'Max', frames, 'Value', 1, 'SliderStep', 1./(max(frames-1, 1))+[0 0], ...
    'Callback', {@updateImage,gcf}, 'BackgroundColor', [0.75 0.75 0.75]);

uicontrol('parent', h.panmov, 'Style', 'text', 'String', 'Frame: ','Units','normalized' ...
    ,'BackgroundColor',pcol,'Position',[0.01 0.15 0.05 0.3],'HorizontalAlignment','Right', ...
    'FontSize', 10);
h.framenumtext(2) = uicontrol('parent', h.panmov, 'Style', 'text', 'String', '1/1 ','Units','normalized' ...
    ,'BackgroundColor',pcol,'Position',[0.065 0.15 0.06 0.3],'HorizontalAlignment','Left', ...
    'FontSize', 10);

h.link_checkbox = uicontrol('parent', h.panaxes, 'Style', 'checkbox', 'Value', 0,'Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.1 0.825 0.8 0.1],'HorizontalAlignment','Right', ...
    'Callback', {@link_checkbox_fcn,gcf}, 'BackgroundColor', pcol, 'String', 'Scroll both axes', 'FontSize', 10);        

h.overlay_checkbox = uicontrol('parent', h.panaxes, 'Style', 'checkbox', 'Value', 0,'Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.1 0.625 0.8 0.1],'HorizontalAlignment','Right', ...
    'Callback', {@overlay_button_fcn,gcf}, 'BackgroundColor', pcol, 'String', 'Overlay', 'FontSize', 10, ...
    'Enable', getEn(all(size(h.moving)==size(h.fixed))));

h.ShowOnt_checkbox = uicontrol('parent', h.panaxes, 'Style', 'checkbox', 'Value', 0,'Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.1 0.425 0.8 0.1],'HorizontalAlignment','Right', ...
    'Callback', {@updateImage,gcf}, 'BackgroundColor', pcol, 'String', 'Show ontology', 'FontSize', 10, ...
    'Enable', getEn(h.haveont));  

h.HideSWC_checkbox = uicontrol('parent', h.panaxes, 'Style', 'checkbox', 'Value', 0,'Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.1 0.225 0.8 0.1],'HorizontalAlignment','Right', ...
    'Callback', {@updateImage,gcf}, 'BackgroundColor', pcol, 'String', 'Hide SWCs', 'FontSize', 10, ...
    'Enable', 'On');  
         
uicontrol('parent', h.clim1_pan, 'Style', 'text', 'String', 'Min Value: ','Units','normalized' ...
    ,'BackgroundColor',pcol,'Position',[0.1 0.715 0.35 0.2],'HorizontalAlignment','Right');

uicontrol('parent', h.clim1_pan, 'Style', 'text', 'String', 'Max value: ','Units','normalized' ...
    ,'BackgroundColor',pcol,'Position',[0.1 0.465 0.35 0.2],'HorizontalAlignment','Right');


h.cmin_edit(1) = uicontrol('parent', h.clim1_pan, 'Style', 'edit', 'String', '0','Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.5 0.75 0.4 0.2],'HorizontalAlignment','Right', ...
    'Callback', {@updateImage,gcf}, 'Enable', 'Off');

h.cmax_edit(1) = uicontrol('parent', h.clim1_pan, 'Style', 'edit', 'String', num2str(cmax(1)),'Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.5 0.5 0.4 0.2],'HorizontalAlignment','Right', ...
    'Callback', {@updateImage,gcf}, 'Enable', 'Off');

h.autolim(1) = uicontrol('parent', h.clim1_pan, 'Style', 'radiobutton',  'Units', 'normalized', 'Position', ...
    [0.25 0.2 0.5 0.2], 'String', 'Auto scale','FontWeight','Bold', 'Min', 0, 'Max', 1, ...
    'BackgroundColor', pcol, 'Value', 1, 'Callback', {@enableLims,gcf, 1});

         
uicontrol('parent', h.clim2_pan, 'Style', 'text', 'String', 'Min Value: ','Units','normalized' ...
    ,'BackgroundColor',pcol,'Position',[0.1 0.715 0.35 0.2],'HorizontalAlignment','Right');

uicontrol('parent', h.clim2_pan, 'Style', 'text', 'String', 'Max value: ','Units','normalized' ...
    ,'BackgroundColor',pcol,'Position',[0.1 0.465 0.35 0.2],'HorizontalAlignment','Right');

h.cmin_edit(2) = uicontrol('parent', h.clim2_pan, 'Style', 'edit', 'String', '0','Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.5 0.75 0.4 0.2],'HorizontalAlignment','Right', ...
    'Callback', {@updateImage,gcf}, 'Enable', 'Off');

h.cmax_edit(2) = uicontrol('parent', h.clim2_pan, 'Style', 'edit', 'String', num2str(cmax(2)),'Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.5 0.5 0.4 0.2],'HorizontalAlignment','Right', ...
    'Callback', {@updateImage,gcf}, 'Enable', 'Off');


h.autolim(2) = uicontrol('parent', h.clim2_pan, 'Style', 'radiobutton',  'Units', 'normalized', 'Position', ...
    [0.25 0.2 0.5 0.2], 'String', 'Auto scale','FontWeight','Bold', 'Min', 0, 'Max', 1, ...
    'BackgroundColor', pcol, 'Value', 1, 'Callback', {@enableLims,gcf, 2});

h.mark = uicontrol('Style', 'pushbutton',  'Units', 'normalized', 'Position', ...
    [.895 0.22 0.1 .04], 'String', 'Add point','FontWeight','Bold', 'Callback', ...
    {@mark_button_fcn,gcf});


ty = 0.00;
h.afWarp = uicontrol('Style', 'pushbutton',  'Units', 'normalized', 'Position', ...
    [.895 0.165-ty 0.1 .03], 'String', 'Affine Warp','FontWeight','Bold', 'Callback', ...
    {@afwarp_button_fcn,gcf});

h.bsWarp = uicontrol('Style', 'pushbutton',  'Units', 'normalized', 'Position', ...
    [.895 0.1325-ty 0.1 .03], 'String', 'b-spline Warp','FontWeight','Bold', 'Callback', ...
    {@bswarp_button_fcn,gcf});

h.SavePts = uicontrol('Style', 'pushbutton',  'Units', 'normalized', 'Position', ...
    [.895 0.04 0.1 .04], 'String', 'Save points','FontWeight','Bold', 'Callback', ...
    {@save_button_fcn,gcf});

h.DeletePts = uicontrol('Style', 'pushbutton',  'Units', 'normalized', 'Position', ...
    [.895 0.095-ty 0.1 .03], 'String', 'Delete point','FontWeight','Bold', 'Callback', ...
    {@delete_button_fcn,gcf});

h.AutoReg_button = uicontrol('Parent', h.panauto, 'Style', 'pushbutton',  'Units', 'normalized', 'Position', ...
    [0.1 0.05 0.8 0.4], 'String', 'Auto register (slow)','FontWeight','Bold', 'Callback', ...
    {@AutoReg_button_fcn,gcf});
h.usePts_checkbox = uicontrol('Parent', h.panauto, 'Style', 'checkbox',  'Units', 'normalized', 'Position', ...
    [0.1 0.75 0.8 0.2], 'String', 'Use points', 'Backgroundcolor', pcol, 'FontSize', 10);

h.msgText = uicontrol('Style', 'text', 'String', 'Ready!','Units','normalized' ...
    ,'BackgroundColor',[1 1 1],'Position',[0.005 0.01 0.44 0.02],'HorizontalAlignment','Center', ...
    'FontSize', 14, 'FontWeight', 'Bold');


cname = {'x fixed'; 'y fixed'; 'z fixed';};
cw = 55;
h.pts_table = uitable(gcf, 'Data', [], 'ColumnName', cname, 'Units', 'Normalized', ...
    'Position', [0.895 0.32 0.1 0.675], 'ColumnEditable', [true true true], 'ColumnFormat', {'Numeric', 'Numeric', 'Numeric'});
set(h.pts_table, 'ColumnWidth', {cw; cw; cw;}', 'CellSelectionCallback', {@CellSelect,gcf});

max_channel=9;
guidata(h.fig, h);

load_fixed_button_fcn(gcf,h.fixed_fn,h.fixed_reg_channel)
load_moving_button_fcn(gcf,h.moving_fn,h.moving_reg_channel)

if load_fid == 1
    load_button_fcn(gcf)
end

close(d)
pause(0.01)
close(fig)


function en = getEn(val)
if val
    en = 'On';
else
    en = 'Off';
end


function ont = importOntology(fn)

fid = fopen(fn);
cnt = 0;
while ~feof(fid)
    [~] = fgetl(fid);
    cnt = cnt + 1;
end
fclose(fid);

id = zeros(cnt, 1);
info = cell(cnt, 1);

fid = fopen(fn);
cnt = 0;
while ~feof(fid)
    str = fgetl(fid);
    cnt = cnt + 1;
    
    inds = strfind(str, ',');
    
    if ~isempty(inds)
    
    id(cnt) = str2double(str(1:inds(1)-1));
    temp = str(inds(1)+1:end);
    
    temp = temp(temp~='"');
    info{cnt} = temp;
    end

end
fclose(fid);

ont.name = info;
ont.id = id;


function enableLims(~, ~, fig, num)

h = guidata(fig);

if (get(h.autolim(num), 'Value')==1)
    set(h.cmin_edit(num), 'Enable', 'Off');
    set(h.cmax_edit(num), 'Enable', 'Off');
else
    set(h.cmin_edit(num), 'Enable', 'On');
    set(h.cmax_edit(num), 'Enable', 'On');
end
updateImage([], [], fig);
guidata(fig, h);


function setFigLims(~, ~, fig, num, im)

h = guidata(fig);

if (get(h.autolim(num), 'Value')==0)
    cmin = str2double(get(h.cmin_edit(num), 'String'));
    cmax = str2double(get(h.cmax_edit(num), 'String'));
else
    a = prctile(double(im(:)), [0.1 99.9]);

    cmin = a(1);
    cmax = a(2);
    
    if cmax<=cmin
        cmax = cmin+1;
    end
    
    set(h.cmin_edit(num), 'String', num2str(cmin));
    set(h.cmax_edit(num), 'String', num2str(cmax));
    
end
set(h.ax(num), 'Clim', [cmin cmax]);
set(h.ax(num+2), 'Clim', [cmin cmax]);
guidata(fig, h);


function delete_button_fcn(~, ~, fig)

h = guidata(fig);

d = get(h.pts_table, 'Data');
inds = setxor(1:size(d,1), h.SelectedRows);
set(h.pts_table, 'Data', d(inds, :));

h.pts.fix = h.pts.fix(inds,:);
h.pts.mov = h.pts.mov(inds,:);
h.pts.cur = h.pts.cur(inds,:);

guidata(fig, h);


function CellSelect(~, eventdata, fig)

h = guidata(fig);
h.SelectedRows = unique(eventdata.Indices(:,1));

if numel(h.SelectedRows)==1
    frame1 = round(h.pts.fix(h.SelectedRows, 3));
    frame2 = round(h.pts.cur(h.SelectedRows, 3));
    
    set(h.frame_slider(1), 'value', frame1);
    set(h.frame_slider(2), 'value', frame2);
    set(h.point(1), "Position", [h.pts.fix(h.SelectedRows, 1),h.pts.fix(h.SelectedRows, 2)]);
    set(h.point(2), "Position", [h.pts.cur(h.SelectedRows, 1),h.pts.cur(h.SelectedRows, 2)]);
end


guidata(fig, h);
updateImage([], [], fig);


function save_button_fcn(~, ~, fig,outfile)
h = guidata(fig);

if numel(h.pts.fix)<3
    disp('No points to save');
    return;
end

[~,name,~] = fileparts(h.moving_fn);
name = split(name,'__');
name = name{1};

def_fn = [name '_Pts' '.mat'];

if ~exist(fullfile(pwd, 'PTS'), 'dir')
    mkdir(fullfile(pwd, 'PTS'));
end

[filepath, ~, ~] = fileparts(h.moving_fn);
[filepath, ~, ~] = fileparts(filepath);
outfile = fullfile(filepath,"internal","fiducials",def_fn);

d = h.pts;
image_files = h.fn;
fiximsz = size(h.fixed);

disp(['Saving ' num2str(size(d.fix, 1)) ' points']);
save(outfile, 'd', 'fiximsz', 'image_files');

h.MainApp.finish_fiducials()

guidata(fig, h);
close(fig)

function save_button_fcn_close_req(~, ~, fig,outfile)
h = guidata(fig);

if numel(h.pts.fix)<3
    disp('No points to save');
    return;
end

[~,name,~] = fileparts(h.moving_fn);
name = split(name,'__');
name = name{1};

def_fn = [name '_Pts' '.mat'];

if ~exist(fullfile(pwd, 'PTS'), 'dir')
    mkdir(fullfile(pwd, 'PTS'));
end

[filepath, ~, ~] = fileparts(h.moving_fn);
[filepath, ~, ~] = fileparts(filepath);
outfile = fullfile(filepath,"internal","fiducials",def_fn);

d = h.pts;
image_files = h.fn;
fiximsz = size(h.fixed);

disp(['Saving ' num2str(size(d.fix, 1)) ' points']);
save(outfile, 'd', 'fiximsz', 'image_files');

h.MainApp.finish_fiducials()

guidata(fig, h);

function load_fixed_button_fcn(fig, fixed_fn, fixed_reg_channel)
load_button_helper(fig, fixed_fn,0,fixed_reg_channel)

function load_moving_button_fcn(fig, moving_fn, moving_reg_channel)
load_button_helper(fig, moving_fn,1, moving_reg_channel)

function load_button_helper(fig, fileName, isMoving, regChannel)
h = guidata(fig);
% loader.fig = figure(125);
% loader.panel = uipanel('Title','Load options','FontSize',12, 'Units', 'Normalized', ...
%     'BackgroundColor',[0.9 0.9 0.9], 'Position',[0.005 0 1 1]);
% 
% loader.HideSWC_checkbox = uicontrol('parent', h.panaxes, 'Style', 'checkbox', 'Value', 0,'Units','normalized' ...
%     ,'BackgroundColor',[1 1 1],'Position',[0.1 0.225 0.8 0.1],'HorizontalAlignment','Right', ...
%     'Callback', {@updateImage,gcf}, 'BackgroundColor', pcol, 'String', 'Hide SWCs', 'FontSize', 10, ...
%     'Enable', 'On');  

index = isMoving+1; %index in all arrays, 1 is fixed, 2 is moving

if length(fileName)<=0
    return
end
if ~exist("regChannel","var")
    regChannel = round(get(h.numChannelsSlider, 'value'));
end
image = read_tiff(fileName,regChannel);
numChannels = getNumChannelsTif(fileName);
% if size(image,3)==1
%     image = cat(3,zeros(size(image)),image,zeros(size(image)));
% end

if isMoving
    h.moving = image;
    h.current = image;
    h.fn.moving = fileName;
    h.fn.movingChan = regChannel;
    h.movingChannels = numChannels;
    h.movingRegChannel = regChannel;
else
    h.fixed = image;
    h.haveont = false;
    h.fn.fixed = fileName;
    h.fn.fixedChan = regChannel;
    h.fixedChannels = numChannels;
    h.fixedRegChannel = regChannel;
    if isfield(h, 'ont')
        h = rmfield(h, "ont");
    end
    if isfield(h, 'an')
        h = rmfield(h, "an");
    end
end

frames = size(image,3);
h.frame_slider(index).Max =frames;
h.frame_slider(index).Value = 1;
h.frame_slider(index).SliderStep = 1./(max(frames-1, 1))+[0 0];

delete(h.axim(index));
h.axim(index) = imagesc(h.ax(index),image(:,:,1)); hold on;

colormap(gray);

a = prctile(image(find(image>0,size(image,1)*size(image,2))), [0.01, 99.9]);
cmin(index) = a(1);
cmax(index) = a(2);

set(h.ax(index),'clim',[cmin(index) cmax(index)]);

if (~isempty(h.pts.fix))
    anyPointOutOfRange=0;
    for i = 1:3
       anyPointOutOfRange =  any(h.pts.fix(:,i)>size(h.fixed,i)) || any(h.pts.mov(:,i)>size(h.moving,i)) || anyPointOutOfRange;
    end
    if (anyPointOutOfRange)
        
        h.tree.mov = [];
        h.treepts = [];
        
        h.pts.fix = [];
        h.pts.mov = [];
        h.pts.cur = [];
        
        set(h.pts_table, 'Data', round(h.pts.fix));
    end
end

delete(h.temppts{1})
delete(h.temppts{2})

h.tree.cur = [];

h.haveaf = 0;
h.havebs = 0;
h.haveont = isfield(h, 'ont')&isfield(h, 'an');
guidata(h.fig, h);
h.point(index).Position = [round(size(image, 2)./2), round(size(image, 1)./2)];
uistack(h.point(index),'top')


if ndims(h.moving) == ndims(h.fixed)
    h.is3D = ndims(h.moving)==3;
else
    h.is3D = 0;
end
guidata(fig, h);
updateImage([],[],fig);




function load_button_fcn(fig)

h = guidata(fig);

% find the directory where the fiducials are stored
[filepath,name,~] = fileparts(h.moving_fn);
[filepath,~,~] = fileparts(filepath);
targetDir = fullfile(filepath,'internal','fiducials');

% obtain uniq id of file
uniq_id = split(name,'__');
uniq_id = uniq_id{1};

% obtain full file name of fiducials
fullfn = find_fn_uniq_id(uniq_id,targetDir);

S = load(fullfn);
if isfield(S, 'd');     d = S.d;        end


if isfield(S,'image_files')
    image_files = S.image_files;
    if isfield (image_files, 'moving') && isfield (image_files, 'movingChan')
        try
            load_button_helper([],[],fig,image_files.moving,1,image_files.movingChan);
        catch
        end
    end
    if isfield (image_files, 'fixed') && isfield (image_files, 'fixedChan')
        try
            load_button_helper([],[],fig,image_files.fixed,0,image_files.fixedChan);
        catch
        end
    end
end
% load_button_helper also edits h and those edits will get overwritten if h
% is retrieved before changing images 

h = guidata(fig);
h.pts = d;
h.pts.cur = h.pts.mov;

set(h.pts_table, 'Data', round(h.pts.fix));

guidata(fig, h);


function tform = loadTransform(fn) %#ok<DEFNU> 
tform.fn = fn;
disp(['Loading ' tform.fn])
fid = fopen(tform.fn);

str =  fscanf(fid, '%s');
oritok = regexp(str, 'ox:(\d*)oy:(\d*)oz:(\d*)', 'tokens');
sztok  = regexp(str, 'sx:(.*)sy:(.*)sz:(.*)', 'tokens');

tform.ox = str2double(oritok{1}{1});
tform.oy = str2double(oritok{1}{2});
tform.oz = str2double(oritok{1}{3});

tform.sx = str2double(sztok{1}{1});
tform.sy = str2double(sztok{1}{2});
tform.sz = str2double(sztok{1}{3});


function bswarp_button_fcn(~, ~, fig)
h = guidata(fig);
lastsz = size(h.current);
tic;

if isempty(h.pts.fix)
    set(h.msgText, 'String', 'No points for warp!!!'); pause(0.1);
    return;
end

% if ~get(h.SkipAf_checkbox, 'Value')
skipAf = true;
if not(skipAf)
    set(h.msgText, 'String', 'Doing affine warp first', 'ForegroundColor', [1 1 1], ...
        'Backgroundcolor', [1 0 0]); pause(0.1);

    [h.pts, h.tree, h.current, h.af] = doAfWarpfast(h.pts, h.tree, h.moving, size(h.fixed));
    
    h.haveaf = 1;
else
    h.tree.cur = h.tree.mov;
    h.pts.cur = h.pts.mov;
    h.current = h.moving;
end

curpts = h.pts.cur(:, [2 1 3]);
fixpts = h.pts.fix(:, [2 1 3]);

set(h.msgText, 'String', 'Solving for b-spline warp', 'ForegroundColor', [1 1 1], ...
        'Backgroundcolor', [1 0 0]); pause(0.1);
    
% options.Verbose=false;
% options.MaxRef = 5;
options.Verbose=true;
options.MaxRef = 3;

[h.O_trans,h.Spacing]         = point_registration(size(h.fixed),fixpts, curpts, options);
[h.O_trans_inv,h.Spacing_inv] = point_registration(size(h.current),curpts, fixpts, options);

h.havebs = 1;

set(h.msgText, 'String', 'Applying b-spline warp', 'ForegroundColor', [1 1 1], ...
        'Backgroundcolor', [1 0 0]); pause(0.1);

h.current = bspline_transform(h.O_trans,h.current,h.Spacing,3);
curpts    = bspline_trans_points_double(h.O_trans_inv, h.Spacing_inv, curpts);
h.pts.cur = curpts(:, [2 1 3]);

if ~isempty(h.tree.cur)
    curtree   = bspline_trans_points_double(h.O_trans_inv, h.Spacing_inv, h.tree.cur(:,[2 1 3]));
    h.tree.cur = curtree(:, [2 1 3]);
end
guidata(fig, h);


if any(lastsz~=size(h.current))
    h.point(2).Position=[size(h.current, 2), size(h.current, 1)]./2;
end

frames = size(h.current, 3);
set(h.frame_slider(2), 'SliderStep', [1./(frames-1) 1./(frames-1)]);
set(h.frame_slider(2), 'Max', frames);

if get(h.frame_slider(2), 'Value')>frames
    set(h.frame_slider(2), 'Value', frames);
end



updateImage([], [], fig);
link_checkbox_fcn([], [], fig);

t = toc;
set(h.msgText, 'String', ['Done Warping (' num2str(t) ' sec)'], 'ForegroundColor', [0 0 0], ...
    'Backgroundcolor', [1 1 1]);

function reloadPy()
warning('off','MATLAB:ClassInstanceExists')
clear classes %#ok<CLCLS> 
mod = py.importlib.import_module('elastixBspline');
py.importlib.reload(mod);

function plot_z_correlations(~,~,main_fig,itx_fig)
itx_h = guidata(itx_fig);
main_h = guidata(main_fig);

% Iterate through every tif file in the directory, check if it's a 3D tif
% and if it is, load it and plot the correlation of each z-plane with the
% provided 2D image.
% Inputs:
%   correlation_image_path: path to the 2D image for correlation
%   tif_folder_path: path to the directory containing the tif files
% Outputs:
%   none

% Load the correlation image
corr_img = main_h.fixed;
tif_folder_path = get(itx_h.outputFolder, 'String');
% Get a list of all the tif files in the directory
tif_files = dir(fullfile(tif_folder_path, '*.tif'));

% Iterate through each tif file
for i = 1:numel(tif_files)
    % Load the tif file
    tif_path = fullfile(tif_files(i).folder, tif_files(i).name);
    tif_info = imfinfo(tif_path);
    
    % Check if the tif file has 3 dimensions
    if numel(tif_info) == 1
        % This is a 2D tif file, skip it
        continue;
    end
    
    % Load the tif file as a 3D array
    tif_data = loadTifFast(tif_path);
    
    % Plot the correlation of each z-plane with the correlation image
    corr_data = zeros(numel(tif_info),1);
    for z = 1:numel(tif_info)
        corr_data(z) = corr2(tif_data(:,:,z), corr_img);
%         imshowpair(tif_data(:,:,z), corr_img, 'montage');
    end
    figure; plot(1:numel(tif_info),corr_data);
    xlabel("plane #");
    ylabel("correlation");
end

guidata(main_fig,main_h)
guidata(itx_fig,itx_h)



function updateChannelSliderLoad(~,~,fig)
h = guidata(fig);

numChannel = round(get(h.numChannelsSlider, 'Value'));
set(h.numChannelText, 'String', num2str(numChannel));


function updateItxSliders(~,~,fig)
itx_h = guidata(fig);

numRes = round(get(itx_h.resSlider, 'Value'));
numIter = round(get(itx_h.iterSlider, 'Value'));
numFrames = round(get(itx_h.framesSlider, 'Value'));
set(itx_h.numResText, 'String', [num2str(numRes) '/' num2str(round(get(itx_h.resSlider, 'Max')))]);
set(itx_h.numIterText, 'String', [num2str(numIter) '/' num2str(round(get(itx_h.iterSlider, 'Max')))]);
set(itx_h.numFramesText, 'String', [num2str(numFrames) '/' num2str(round(get(itx_h.framesSlider, 'Max')))]);

function updateOutputFolder(~,~,fig)
itx_h = guidata(fig);
selected_dir = uigetdir();
if (selected_dir~=0)
    set(itx_h.outputFolder, 'String', selected_dir)
end

function updateApplyFolder(~,~,fig)
itx_h = guidata(fig);
selected_dir = uigetdir();
if (selected_dir~=0)
    set(itx_h.applyFolder, 'String', selected_dir)
end


function outString = addQuotes(inString)
outString = sprintf('"%s"', inString);


function intensity_bswarp_fcn(~,~,main_fig,itx_fig,matlab_python_working)
h = guidata(main_fig);
itx_h = guidata(itx_fig);
tic;



if ~get(h.SkipAf_checkbox, 'Value')
    set(h.msgText, 'String', 'Doing affine warp first', 'ForegroundColor', [1 1 1], ...
        'Backgroundcolor', [1 0 0]); pause(0.1);

    [h.pts, h.tree, h.current, h.af] = doAfWarpfast(h.pts, h.tree, h.moving, size(h.fixed));
    
    h.haveaf = 1;
else
    h.tree.cur = h.tree.mov;
    h.pts.cur = h.pts.mov;
    h.current = h.moving;
end


% try
initialTransformFileName = "initial_transform_vals.txt";
fixedImFilename = h.fn.fixed;
movingImFilename = h.fn.moving;
numRes = round(get(itx_h.resSlider, 'Value'));
numIter = round(get(itx_h.iterSlider, 'Value'));
outFolder = get(itx_h.outputFolder, 'String');
inputFolder = get(itx_h.applyFolder, 'String');
outFile = get(itx_h.outputFileName, 'String');
numChannelsInToTransform = getNumChannelsTifDirectory(inputFolder);

numChannelsInFixed = h.fixedChannels;
numChannelsInMoving = h.movingChannels;
regChannelFixed = h.fixedRegChannel;
regChannelMoving = h.movingRegChannel;

save_button_fcn([],[],main_fig,fullfile(outFolder,['Pts_' date '.mat']));

if any(outFile(end-3:end)~='.tif')
    outFile = [outFile '.tif'];
end

affine = inv(h.af);
if h.is3D
    affine = affine([1:3 5:7 9:11 4 8 12]);
else
    affine = affine([1:2 4:5 3 6]);
end

mkdir(inputFolder,'points');
pointsDir = fullfile(inputFolder,'points');
fixedFilepts = fullfile(pointsDir,"ptsFixed.txt");
writePointsFile(fixedFilepts,h.pts.fix)
movFilepts = fullfile(pointsDir,"ptsMov.txt");
writePointsFile(movFilepts,h.pts.mov)

% List .mat files in the input directory
files = dir(fullfile(inputFolder, '*.mat'));

% Check if there's a single .mat file in the directory
if length(files) == 1
    % Full path to the .mat file
    matFile = fullfile(files(1).folder, files(1).name);
    
    % Load the .mat file into 'rois' variable
    rois = load(matFile).roi;
    for z = 1:size(rois, 3)
        % Label each separate ROI in the current z-plane
        labeledImage = logical(rois(:, :, z));

        % Get properties of each labeled region
        props = regionprops(labeledImage, 'Centroid');

        % Store centroids in the cell array
        centroids{z} = cat(1, props.Centroid);
    end
    % figure;imagesc(rois(:,:,1));hold on; plot(centroids{1}(1),centroids{1}(2),'.');
end


if matlab_python_working
    im = h.pythonModule.doIntensityBspline(initialTransformFileName,fixedImFilename,movingImFilename,inputFolder,...
        affine,fullfile(outFolder,outFile),num2str(numRes),num2str(numIter),numChannelsInFixed,...
        numChannelsInMoving,numChannelsInToTransform,regChannelFixed,regChannelMoving);
    
else
    argsIn = {addQuotes(initialTransformFileName)};
    argsIn(end+1) = {addQuotes(fixedImFilename)};
    argsIn(end+1) = {addQuotes(movingImFilename)};
    argsIn(end+1) = {addQuotes(inputFolder)};
    argsIn(end+1) = {addQuotes(regexprep(num2str(affine),' +',','))};
    argsIn(end+1) = {addQuotes(fullfile(outFolder,outFile))};
    argsIn(end+1) = {addQuotes(num2str(numRes))};
    argsIn(end+1) = {addQuotes(num2str(numIter))};
    argsIn(end+1) = {addQuotes(num2str(numChannelsInFixed))};
    argsIn(end+1) = {addQuotes(num2str(numChannelsInMoving))};
    argsIn(end+1) = {addQuotes(num2str(numChannelsInToTransform))};
    argsIn(end+1) = {addQuotes(num2str(regChannelFixed))};
    argsIn(end+1) = {addQuotes(num2str(regChannelMoving))};
    
    system("python callElastixBspline.py " + strjoin(argsIn, ' '));

    
    
%     system("python callElastixBspline.py " + join(argsIn, ' '));
end


if exist(fixedFilepts, 'file')
    delete(fixedFilepts);
end

if exist(movFilepts, 'file')
    delete(movFilepts);
end

if exist(pointsDir, 'dir')
    rmdir(pointsDir);
end

% catch
% end


set(h.msgText, 'String', 'Solving for b-spline warp', 'ForegroundColor', [1 1 1], ...
        'Backgroundcolor', [1 0 0]); pause(0.1);





link_checkbox_fcn([], [], main_fig);

t = toc;
set(h.msgText, 'String', ['Done Warping (' num2str(t) ' sec)'], 'ForegroundColor', [0 0 0], ...
    'Backgroundcolor', [1 1 1]);



function afwarp_button_fcn(~, ~, fig)
tic;
h = guidata(fig);

if isempty(h.pts.fix)
    set(h.msgText, 'String', 'No points for warp!!!'); pause(0.1);
    return;
end

set(h.msgText, 'String', 'Doing affine warp', 'ForegroundColor', [1 1 1], ...
        'Backgroundcolor', [1 0 0]); pause(0.1);
fiximsz = size(h.fixed);
[h.pts, h.tree, h.current, h.af] = doAfWarpfast(h.pts, h.tree, h.moving, fiximsz);

h.haveaf = 1;

% setPosition(h.point(2),[size(h.current, 1), size(h.current, 2)]./2);

if h.is3D
    frames = size(h.current, 3);
    set(h.frame_slider(2), 'SliderStep', [1./(frames-1) 1./(frames-1)]);
    set(h.frame_slider(2), 'Max', frames);
    
    if get(h.frame_slider(2), 'Value')>frames
        set(h.frame_slider(2), 'Value', frames);
    end
end
 
guidata(fig, h);

updateImage([], [], fig);
link_checkbox_fcn([], [], fig);

t=toc;
set(h.msgText, 'String', ['Done Warping (' num2str(t) ' sec)'], 'ForegroundColor', [0 0 0], ...
    'Backgroundcolor', [1 1 1]); pause(0.1);


function AutoReg_button_fcn(~, ~, fig)

tic;
h = guidata(fig);

set(h.msgText, 'String', 'Doing auto warp', 'ForegroundColor', [1 1 1], ...
        'Backgroundcolor', [1 0 0]); pause(0.1);

Options.Penalty = 1e-5;
Options.Verbose = 2;
Options.Similarity = 'mi';

if get(h.usePts_checkbox, 'Value')
    if isempty(h.pts.fix)
        set(h.msgText, 'String', 'No points for warp!!!'); pause(0.1);
        return;
    end
    
    Options.Points1 = h.pts.cur;
    Options.Points2 = h.pts.fix;
end
    
[h.current,h.auto.O_trans,h.auto.Spacing, h.auto.M, h.auto.B,h.auto.F] = image_registration(h.current, h.moving, Options);

guidata(fig, h);
updateImage([], [], fig);
t=toc;
set(h.msgText, 'String', ['Done Warping (' num2str(t) ' sec)'], 'ForegroundColor', [0 0 0], ...
    'Backgroundcolor', [1 1 1]); pause(0.1);   

function exposedDoAfWarpfast(fig,outFolder)
h = guidata(fig);

if ~get(h.SkipAf_checkbox, 'Value')
    set(h.msgText, 'String', 'Doing affine warp first', 'ForegroundColor', [1 1 1], ...
        'Backgroundcolor', [1 0 0]); pause(0.1);

    [h.pts, h.tree, h.current, h.af] = doAfWarpfast(h.pts, h.tree, h.moving, size(h.fixed));
    
    h.haveaf = 1;
else
    h.tree.cur = h.tree.mov;
    h.pts.cur = h.pts.mov;
    h.current = h.moving;
end
save_button_fcn([],[],fig,fullfile(outFolder,['Pts_' date '.mat']));
set(h.msgText, 'String', 'Done Warping', 'ForegroundColor', [0 0 0], ...
                'Backgroundcolor', [1 1 1]);
guidata(fig,h);

function [pts, tree, im, af] = doAfWarpfast(pts, tree, im, fiximsz)

if (size(pts.fix,1)<3)
    return
end
imsz = size(im);

imold = im;
if numel(fiximsz)<3
    [pts, tree, im, af] = doAfWarpfast2D(pts, tree, im, fiximsz);
    return;
end
im = zeros(max(fiximsz(1), imsz(1)), max(fiximsz(2), imsz(2)), max(fiximsz(3), imsz(3)));
im(1:imsz(1), 1:imsz(2), 1:imsz(3)) = imold;
clear imold;

posf = pts.fix;
posm = pts.mov;


posm = [posm(:, [2 1 3]) ones(size(posm, 1), 1)];
posf = [posf(:, [2 1 3]) ones(size(posf, 1), 1)];

mid = [size(im)./2 1];

posm(:,1) = posm(:,1) - mid(1);
posm(:,2) = posm(:,2) - mid(2);
posm(:,3) = posm(:,3) - mid(3);

af = posm\posf;

af(4,1) = af(4,1) - mid(1);
af(4,2) = af(4,2) - mid(2);
af(4,3) = af(4,3) - mid(3);

im = affine_transform(im, inv(af'),1);

padx = round(fiximsz(2) - imsz(2));
pady = round(fiximsz(1) - imsz(1));
padz = round(fiximsz(3) - imsz(3));

dxn = round(padx.*(padx<0));
dyn = round(pady.*(pady<0));
dzn = round(padz.*(padz<0));

im = im(1:end+dyn, 1:end+dxn, 1:end+dzn);

posf = pts.fix;
posm = pts.mov;

posm = [posm ones(size(posm, 1), 1)];
posf = [posf ones(size(posf, 1), 1)];

af = posm\posf;

af(abs(af)<1e-9) = 0;
af(abs(af-1)<1e-9) = 1;

a = [pts.mov, ones(size(pts.cur,1), 1)];
a = a*af;
pts.cur = a(:,1:3);

if ~isempty(tree.mov)
    a = [tree.mov, ones(size(tree.mov,1), 1)];
    a = a*af;
    tree.cur = a(:,1:3);
end






function [pts, tree, im, af] = doAfWarpfast2D(pts, tree, im, fiximsz)

if (length(pts.fix)<3)
    return
end
imsz = size(im);

% imold = im;
% im = zeros(max(fiximsz(1), imsz(1)), max(fiximsz(2), imsz(2)));
% im(1:imsz(1), 1:imsz(2)) = imold;
% clear imold;

posf = pts.fix(:,1:2);
posm = pts.mov(:,1:2);

posm = [posm(:, [1 2]) ones(size(posm, 1), 1)];
posf = [posf(:, [1 2]) ones(size(posf, 1), 1)];

% mid = [fiximsz(1)/2 fiximsz(2)/2 1];

% posm(:,1) = posm(:,1) - mid(1);
% posm(:,2) = posm(:,2) - mid(2);


af = posm\posf;

% af(3,1) = af(3,1) - mid(1);
% af(3,2) = af(3,2) - mid(2);


% im = affine_transform(im, inv(af'),1);

im = imwarp(im,affinetform2d(af'),OutputView=imref2d(fiximsz));

% padx = round(fiximsz(2) - imsz(2));
% pady = round(fiximsz(1) - imsz(1));
% % padz = round(fiximsz(3) - imsz(3));
% 
% dxn = round(padx.*(padx<0));
% dyn = round(pady.*(pady<0));
% dzn = round(padz.*(padz<0));

% im = im(1:end+dyn, 1:end+dxn);
posf = pts.fix;
posm = pts.mov;

posm = [posm(:, [1 2]) ones(size(posm, 1), 1)];
posf = [posf(:, [1 2]) ones(size(posf, 1), 1)];


af = posm\posf;

af(abs(af)<1e-9) = 0;
af(abs(af-1)<1e-9) = 1;

a = [pts.mov(:,1:2), ones(size(pts.cur,1), 1)];
a = a*af;
pts.cur = [a(:,1:2) ones(size(pts.cur,1),1)];

if ~isempty(tree.mov)
    a = [tree.mov, ones(size(tree.mov,1), 1)];
    a = a*af;
    tree.cur = a(:,1:3);
end

function link_checkbox_fcn(~, ~, fig)
h = guidata(fig);

if ndims(h.fixed)==ndims(h.current) && all(size(h.fixed)==size(h.current))
    set(h.overlay_checkbox, 'Enable', 'on');
end

if get(h.link_checkbox, 'Value')
    frame = get(h.frame_slider(1), 'Value');
    if frame>size(h.current,3)
        frame = size(h.current,3);
    end
    
    set(h.frame_slider(2), 'Value', frame);
    
end
% guidata(fig, h);
updateImage([], [], fig);


function mark_button_fcn(~,~,fig)
h = guidata(fig);

p1 = h.point(1).Position;
p2 = h.point(2).Position;

z1 = round(get(h.frame_slider(1), 'Value'));
z2 = round(get(h.frame_slider(2), 'Value'));

h.pts.fix(end+1, :) = [p1 z1];
h.pts.cur(end+1, :) = [p2 z2];

a = [p2 z2 1];
if h.havebs
    a(1:3) = bspline_trans_points_double(h.O_trans, h.Spacing, a([2 1 3]));
    a(1:3) = a([2 1 3]);
end
if h.haveaf
%     afinv = inv(h.af);
%     a = a([1:2 4])/h.af;
    a = a/h.af;
end
h.pts.mov(end+1, :) = a(1:3);

% tp1= plot(h.ax(1), p1(1), p1(2), 'g.', 'MarkerSize', 10);
% tp2= plot(h.ax(2), p2(1), p2(2), 'g.', 'MarkerSize', 10);

% if get(h.overlay_checkbox, 'Value')
%     c = get(h.ax(2), 'Children');
%     c = [c(2:3); tp2; c(4:end)];
%     set(h.ax(2), 'Children', c);
% else
%     c = get(h.ax(2), 'Children');
%     c = [c(2); tp2; c(3:end)];
%     set(h.ax(2), 'Children', c);
% end
% 
% c = get(h.ax(1), 'Children');
% c = [c(2); tp1; c(3:end)];
% set(h.ax(1), 'Children', c);

% h.temppts{1}(end+1) = tp1;
% h.temppts{2}(end+1) = tp2;

% d = get(h.pts_table, 'Data');
% newd = [d; p1(1) p2(1) p1(2) p2(2) z1 z2];
set(h.pts_table, 'Data', round(h.pts.fix));

%for ref atlas stuff
if (h.haveont)
    query = [round(p1) round(z1)];
    val = h.an(query(2), query(1), query(3));
    id = find(h.ont.id==val);
    if (numel(id)==1)
        loc = h.ont.name{id};
        
    else
        loc = 'outer space';
    end
    set(h.msgText, 'String', ['Marked point in #' num2str(size(h.pts.mov,1)) ' in ' loc]);
else
    set(h.msgText, 'String', ['Marked point #' num2str(size(h.pts.mov,1))]);
end
    
  
guidata(fig,h);
updateImage([], [], fig);


function contour_button_fcn(~,~,fig, yval)
h = guidata(fig);

hand = impoly(h.ax(2), 'Closed', 0);
pos = getPosition(hand);
delete(hand);

z2 = round(get(h.frame_slider(2), 'Value'));

movpos = [pos z2.*ones(size(pos, 1), 1)];

fixpos = zeros(size(movpos));
xpos = 0;
for i = 1:size(movpos, 1)
    fixpos(i,:) = [xpos, yval, z2];
    
    if i<size(movpos, 1)
        xpos = xpos + sqrt((movpos(i+1, 1)-movpos(i, 1)).^2 + (movpos(i+1, 2)-movpos(i, 2)).^2);
    end
end

h.pts.fix(end+1:end+size(fixpos, 1), :) = fixpos;
h.pts.cur(end+1:end+size(fixpos, 1), :) = movpos;

a = [movpos ones(size(movpos, 1), 1)];
if h.havebs
    a(:,1:3) = bspline_trans_points_double(h.O_trans, h.Spacing, a(:,[2 1 3]));
    a(:,1:3) = a(:,[2 1 3]);
end
if h.haveaf
%     afinv = inv(h.af);
    a = a/h.af;
end
h.pts.mov(end+1:end+size(fixpos, 1), :) = a(:,1:3);


set(h.pts_table, 'Data', round(h.pts.fix));

%for ref atlas stuff
if (h.haveont)
    query = [round(p1) round(z1)];
    val = h.an(query(2), query(1), query(3));
    id = find(h.ont.id==val);
    if (numel(id)==1)
        loc = h.ont.name{id};
        
    else
        loc = 'outer space';
    end
    set(h.msgText, 'String', ['Marked point in #' num2str(size(h.pts.mov,1)) ' in ' loc]);
else
    set(h.msgText, 'String', ['Marked point #' num2str(size(h.pts.mov,1))]);
end
    
  
guidata(fig,h);
updateImage([], [], fig);


function mid = findThresh(x, y)

% zeroinds = abs(y<1);
% y(zeroinds) = -10*std(y);

n10 = ceil(numel(x)./5);
amp1 = mean(y(1:n10));
amp2 = mean(y(end-n10+1:end));

% smallamp = min(amp1, amp2);
yf = y;
flip=0;
if amp1>amp2
    yf = yf(end:-1:1);
    flip=1;
end

amp1 = mean(yf(1:n10));
amp2 = mean(yf(end-n10+1:end));

amp = amp2-amp1;
offset = amp1;
mid = numel(x)./2;

p0 = [amp mid];

pf = nlinfit(x, yf, @(A, x)sigfunc(A, x, offset), p0);
mid = find(yf>offset+pf(1)./2, 1, 'first');

if flip
    mid = numel(x)-mid+1;
end

% figure; plot(x, y);
% yy = sigfunc(pf, x, offset);
% hold on; plot(x, yy, 'r', x(round(mid)), y(round(mid)), 'k.')



function y = sigfunc(A, x, offset)
y = offset + A(1)./ (1 + exp(-(x-A(2))./3));


function figScroll(~,evnt,fig)

h = guidata(fig);

[pos1(1), pos1(2)] = gpos(h.ax(1));
[pos2(1), pos2(2)] = gpos(h.ax(2));

dim1 = size(h.fixed);
dim2 = size(h.current);

xhit(1) = (pos1(1)>0)&pos1(1)<dim1(2);
yhit(1) = (pos1(2)>0)&pos1(2)<dim1(1);

xhit(2) = (pos2(1)>0)&pos2(1)<dim2(2);
yhit(2) = (pos2(2)>0)&pos2(2)<dim2(1);

hits = (xhit.*yhit)>0.5;
if (sum(hits)~=1)
    return;
end
axind = find(hits);

if axind==1
    maxframe = size(h.fixed,3);
end
if axind==2
    maxframe = size(h.current,3);
end


frame = round(get(h.frame_slider(axind), 'Value'));
frame = frame+evnt.VerticalScrollCount;
if frame<1
    frame = 1;
end
if frame>maxframe
    frame = maxframe;
end

set(h.frame_slider(axind), 'value', frame);
movebothpoints(fig, axind);

    
    
    


if get(h.link_checkbox, 'Value')
    axind = find(~hits);
    movebothpoints(fig, axind);
    if axind==1
        maxframe = size(h.fixed,3);
    end
    if axind==2
        maxframe = size(h.current,3);
    end

    if frame>maxframe
        frame = maxframe;
    end
    set(h.frame_slider(axind), 'value', frame);
    
end

updateImage([], [], fig);


function figClick(~,~,fig)


h = guidata(fig);

if ~h.haveont
    return;
end

[pos1(1), pos1(2)] = gpos(h.ax(1));
[pos2(1), pos2(2)] = gpos(h.ax(2));

newpoint=get(fig,'CurrentPoint');

dim1 = size(h.fixed);
dim2 = size(h.current);

xhit(1) = (pos1(1)>0)&pos1(1)<dim1(2);
yhit(1) = (pos1(2)>0)&pos1(2)<dim1(1);

xhit(2) = (pos2(1)>0)&pos2(1)<dim2(2);
yhit(2) = (pos2(2)>0)&pos2(2)<dim2(1);

hits = (xhit.*yhit)>0.5;
if (sum(hits)~=1)
    return;
end
axind = find(hits);

axpos = plotboxpos(h.ax(axind));

newpoint(1) = (newpoint(1)-axpos(1))./axpos(3);
newpoint(2) = (axpos(2)+axpos(4)-newpoint(2))./axpos(4);

chck = any(newpoint>1)|any(newpoint<0);
if chck
    return;
end

sz = (size(h.current,1)==size(h.fixed,1))&(size(h.fixed,2)==size(h.current,2));

if (axind==2)&&(~sz)
    return;
end


if axind==1
    %     query = pos1;
    newpoint(1) = newpoint(1).*size(h.fixed, 2);
    newpoint(2) = newpoint(2).*size(h.fixed, 1);
    
    query = newpoint;
end
if axind==2
    %     query = pos2;
    newpoint(1) = newpoint(1).*size(h.current, 2);
    newpoint(2) = newpoint(2).*size(h.current, 1);
    query = newpoint;
    
end

frame = round(get(h.frame_slider(axind), 'Value'));

query = [ceil(query(2)) ceil(query(1)) frame];
%for ref atlas stuff

sz = size(h.an);
query(query>sz) = sz(query>sz);

val = h.an(query(1), query(2), query(3));
id = find(h.ont.id==val);
if (numel(id)==1)
    loc = h.ont.name{id};
    
else
    loc = 'outer space';
end
set(h.msgText, 'String', ['Mouse click in ' loc ' (# ' num2str(val) ')' ]);


function updateImage(~, ~, fig)
h = guidata(fig);

if ndims(h.fixed)==3
    frame1 = round(get(h.frame_slider(1), 'Value'));
else
    frame1=1;
end
pos1 = round(h.point(1).Position);

if get(h.ShowOnt_checkbox, 'Value')
    tempfix = h.an(:,:,frame1);
    tempfixsag = squeeze(h.an(:, pos1(1), :));
else
    tempfix = h.fixed(:,:,frame1);
    tempfixsag = squeeze(h.fixed(:, pos1(1), :));
end
set(h.axim(1), 'Cdata', tempfix);
set(h.axim(3), 'Cdata', tempfixsag);
set(h.framenumtext(1), 'String', [num2str(frame1) '/' num2str(size(h.fixed,3))]);

setFigLims([], [], fig, 1, tempfix);

if ndims(h.current)==3
    frame2 = round(get(h.frame_slider(2), 'Value'));
else
    frame2=1;
end
pos2 = round(h.point(2).Position);
set(h.framenumtext(2), 'String', [num2str(frame2) '/' num2str(size(h.current,3))]);

tempmov = double(h.current(:,:,frame2));
if (h.is3D)
    tempmovsag = double(squeeze(h.current(:,pos2(1),:)));
end

over = get(h.overlay_checkbox, 'Value')&&all(size(h.current)==size(h.fixed));
setFigLims([], [], fig, 2, tempmov);


if (over)
    
    tempfix = double(tempfix);
    tempfixsag = double(tempfixsag);
    
    clim(1) = str2double(get(h.cmin_edit(1), 'String'));
    clim(2) = str2double(get(h.cmax_edit(1), 'String'));
    
    tempfix = tempfix-clim(1);
    tempfix = tempfix./(clim(2)-clim(1));
    tempfix(tempfix>1) = 1;
    tempfix(tempfix<0) = 0;
    
    tempfixsag = tempfixsag-clim(1);
    tempfixsag = tempfixsag./(clim(2)-clim(1));
    tempfixsag(tempfixsag>1) = 1;
    tempfixsag(tempfixsag<0) = 0;
    
    clim(1) = str2double(get(h.cmin_edit(2), 'String'));
    clim(2) = str2double(get(h.cmax_edit(2), 'String'));
    
    tempmov = tempmov-clim(1);
    tempmov = tempmov./(clim(2)-clim(1));
    tempmov(tempmov>1) = 1;
    tempmov(tempmov<0) = 0;
    
    if (h.is3D)
        tempmovsag = tempmovsag-clim(1);
        tempmovsag = tempmovsag./(clim(2)-clim(1));
        tempmovsag(tempmovsag>1) = 1;
        tempmovsag(tempmovsag<0) = 0;
        set(h.axim(4), 'Cdata', cat(3, tempmovsag, tempfixsag, tempmovsag));
    end
    set(h.axim(2), 'Cdata', cat(3, tempmov, tempfix, tempmov));
else
    set(h.axim(2), 'Cdata', tempmov);
    if (h.is3D)
        set(h.axim(4), 'Cdata', tempmovsag);
    end
end


% axis(h.ax, 'image'); I don't think this call needs to be here, maybe when
% loading things for the first time, but it resets zoom annoyingly

if ~isempty(h.pts.fix)
    inds{1} = find(abs(h.pts.fix(:,3)-frame1)<5.51);
    inds{2} = find(abs(h.pts.cur(:,3)-frame2)<5.51);

    ov = get(h.overlay_checkbox, 'Value');
    
    h.temppts{1} = plotPointsOnAxis(h.ax(1), h.temppts{1}, inds{1}, h.pts.fix, 1, 0, 30, '.', zeros(numel(inds{1}), 1));  
    h.temppts{2} = plotPointsOnAxis(h.ax(2), h.temppts{2}, inds{2}, h.pts.cur, 1+ov, 0, 30, '.', zeros(numel(inds{2}), 1)); 
end

if ~isempty(h.tree.cur)
    ov = get(h.overlay_checkbox, 'Value');
    
    if get(h.HideSWC_checkbox, 'Value')
        inds = [];
    else
        inds = find(abs(h.tree.cur(:,3)-frame2)<h.max_depth);
    end
    
    dist = abs(h.tree.cur(inds,3) -  frame2);
    
    %         if ~isempty(inds);  h.treepts = plotPointsOnAxis(h.ax(2), h.treepts, inds, h.tree, 1+ov, [0 1 1], 10, '+'); end
    h.treepts = plotPointsOnAxis(h.ax(2), h.treepts, inds, h.tree.cur, 1+ov, 0.5, 20, '.', dist);    
end

guidata(fig, h);


function tp = plotPointsOnAxis(ax, tp, inds, pts, numdots, hue, sz, mkr, dist) %#ok<INUSD,INUSL> 
% tp = scatter(ax, pts(inds, 1), pts(inds, 2), 5.*ones(size(inds)), repmat(clr, numel(inds), 1));

if ishandle(tp)
%     numchildbefore = numel(get(ax, 'Children'));
%     set(tp, 'Xdata', pts(inds, 1), 'Ydata', pts(inds,2), 'MarkerSize', sz, 'MarkerFaceColor', clr, 'MarkerEdgeColor', clr);
%     a = hsv2rgb([hue*ones(numel(inds), 1), ones(numel(inds), 1), 1-dist./8]);
    a = [1 0 0];
    set(tp, 'Xdata', pts(inds, 1), 'Ydata', pts(inds,2), 'CData', a, 'SizeData', sz);
%     numchildafter = numel(get(ax, 'Children'))
%     numel(get(tp, 'Children'))
%     numel(inds)
%     
%     
%     c = cell(numel(inds), 1);
%     for i = 1:numel(c)
%         c{i} = 1 - dist(i)./8;
%     end
%     get(tp)
%     set(get(tp, 'Children'), {'FaceAlpha'}, c);
else
%     tp = plot(ax, pts(inds, 1), pts(inds, 2), mkr, 'MarkerSize', sz, 'MarkerFaceColor', clr, 'MarkerEdgeColor', clr);
%     a = hsv2rgb([hue*ones(numel(inds), 1), ones(numel(inds), 1), 1-dist./8]);
    a = [ones(numel(inds, 1)), ones(numel(inds, 1)), ones(numel(inds, 1))];
    tp = scatter(ax, pts(inds, 1), pts(inds, 2), sz.*ones(size(inds)), a, 'filled');
    numchild = 1;
    c = get(ax, 'Children');
    c = [c(numchild+1:numchild+1+numdots-1); c(1:numchild); c(numchild+numdots+1:end)];
    set(ax, 'Children', c);

end


function overlay_button_fcn(~, ~, fig)
h = guidata(fig);
if get(h.overlay_checkbox, 'Value')
    h.point(3) = drawpoint(h.ax(2), 'Color', 'g', 'Position', h.point(1).Position);
    addlistener(h.point(3),"MovingROI",@(src,pos) movebothpoints(fig, 3));
else
    delete(h.point(3));
end
h.point(2).Position = h.point(1).Position;
bringToFront(h.point(2));
guidata(fig, h);
updateImage([], [], fig);


function movebothpoints(fig, num)
h = guidata(fig);
pos = h.point(num).Position;
if num==1 || num==2
    frame = round(get(h.frame_slider(num), 'Value'));
    set(h.plus(num), 'XData', frame, 'Ydata', pos(2));
end

if num==1 || num==3
    other_point = 4-num; %corresponding mirror of point
    if isvalid(h.point(other_point))
        h.point(other_point).Position = pos;
    end
end
% updateImage([], [], fig);


function writeTiff16(im, fn)
imwrite(uint16(squeeze(im(:,:,1, :))), fn, 'tif', 'writemode', 'overwrite', 'Compression', 'None');
for k = 2:size(im,3)
    imwrite(uint16(squeeze(im(:,:,k, :))), fn, 'tif','writemode', 'append', 'Compression', 'None');
end


function writeTiff32(im, fn)
imwrite(squeeze(im(:,:,1, :)), fn, 'tif', 'writemode', 'overwrite', 'Compression', 'None');
for k = 2:size(im,3)
    imwrite(squeeze(im(:,:,k, :)), fn, 'tif','writemode', 'append', 'Compression', 'None');
end


function writeTiff8(im, fn)
imwrite(uint8(squeeze(im(:,:,1, :))), fn, 'tif', 'writemode', 'overwrite', 'Compression', 'None');
for k = 2:size(im,3)
    imwrite(uint8(squeeze(im(:,:,k, :))), fn, 'tif','writemode', 'append', 'Compression', 'None');
end


function im = loadData(fn, par)

% [~ , im.an, im.fib, ~] = getAllenData(0, 1, 1, 1);  im.an = permute(im.an, [2 3 1]);  im.fib = permute(im.fib, [2 3 1]); 

if isfield(fn, 'ont')
    im.ont = importOntology(fullfile(par, fn.ont));
end
if isfield(fn, 'an')
    im.an = loadTifFast(fullfile(par, fn.an));
end


if isfield(fn, 'moving')
    im.moving = loadTifFast(fullfile(par, fn.moving));
else
    im.moving = rand(100,100,25);
end

if isfield(fn, 'fixed')
    im.fixed = loadTifFast(fullfile(par, fn.fixed));
else
    im.fixed = zeros(size(im.moving));
end

im.current = im.moving;

function [x,y]=gpos(h_axes)
%GPOS Get current position of cusor and return its coordinates in axes with handle h_axes
% h_axes - handle of specified axes
% [x,y]  - cursor coordinates in axes h_aexs
%
% -------------------------------------------------------------------------
% Note:
%  1. This function should be called in the figure callback WindowButtonMotionFcn.
%  2. It works like GINPUT provided by Matlab,but it traces the position
%       of cursor without click and is designed for 2-D axes.
%  3. It can also work even the units of figure and axes are inconsistent,
%       or the direction of axes is reversed.
% -------------------------------------------------------------------------

% Written by Kang Zhao,DLUT,Dalian,CHINA. 2003-11-19
% E-mail:kangzhao@student.dlut.edu.cn

h_figure=gcf;

units_figure = get(h_figure,'units');
units_axes   = get(h_axes,'units');

if_units_consistent = 1;

if ~strcmp(units_figure,units_axes)
    if_units_consistent=0;
    set(h_axes,'units',units_figure); % To be sure that units of figure and axes are consistent
end

% Position of origin in figure [left bottom]
pos_axes_unitfig    = get(h_axes,'position');
width_axes_unitfig  = pos_axes_unitfig(3);
height_axes_unitfig = pos_axes_unitfig(4);

xDir_axes=get(h_axes,'XDir');
yDir_axes=get(h_axes,'YDir');

% Cursor position in figure
pos_cursor_unitfig = get( h_figure, 'currentpoint'); % [left bottom]

if strcmp(xDir_axes,'normal')
    left_origin_unitfig = pos_axes_unitfig(1);
    x_cursor2origin_unitfig = pos_cursor_unitfig(1) - left_origin_unitfig;
else
    left_origin_unitfig = pos_axes_unitfig(1) + width_axes_unitfig;
    x_cursor2origin_unitfig = -( pos_cursor_unitfig(1) - left_origin_unitfig );
end

if strcmp(yDir_axes,'normal')
    bottom_origin_unitfig     = pos_axes_unitfig(2);
    y_cursor2origin_unitfig = pos_cursor_unitfig(2) - bottom_origin_unitfig;
else
    bottom_origin_unitfig = pos_axes_unitfig(2) + height_axes_unitfig;
    y_cursor2origin_unitfig = -( pos_cursor_unitfig(2) - bottom_origin_unitfig );
end

xlim_axes=get(h_axes,'XLim');
width_axes_unitaxes=xlim_axes(2)-xlim_axes(1);

ylim_axes=get(h_axes,'YLim');
height_axes_unitaxes=ylim_axes(2)-ylim_axes(1);

x = xlim_axes(1) + x_cursor2origin_unitfig / width_axes_unitfig * width_axes_unitaxes;
y = ylim_axes(1) + y_cursor2origin_unitfig / height_axes_unitfig * height_axes_unitaxes;

% Recover units of axes,if original units of figure and axes are not consistent.
if ~if_units_consistent
    set(h_axes,'units',units_axes); 
end


function varargout = csvimport( fileName, varargin )
% CSVIMPORT reads the specified CSV file and stores the contents in a cell array or matrix
%
% The file can contain any combination of text & numeric values. Output data format will vary
% depending on the exact composition of the file data.
%
% CSVIMPORT( fileName ):         fileName     -  String specifying the CSV file to be read. Set to
%                                                [] to interactively select the file.
%
% CSVIMPORT( fileName, ... ) : Specify a list of options to be applied when importing the CSV file.
%                              The possible options are:
%                                delimiter     - String to be used as column delimiter. Default
%                                                value is , (comma)
%                                columns       - String or cell array of strings listing the columns
%                                                from which data is to be extracted. If omitted data
%                                                from all columns in the file is imported. If file
%                                                does not contain a header row, the columns
%                                                parameter can be a numeric array listing column
%                                                indices from which data is to be extracted.
%                                outputAsChar  - true / false value indicating whether the data
%                                                should be output as characters. If set to false the
%                                                function attempts to convert each column into a
%                                                numeric array, it outputs the column as characters
%                                                if conversion of any data element in the column
%                                                fails. Default value is false.
%                                uniformOutput - true / false value indicating whether output can be
%                                                returned without encapsulation in a cell array.
%                                                This parameter is ignored if the columns / table
%                                                cannot be converted into a matrix.
%                                noHeader      - true / false value indicating whether the CSV
%                                                file's first line contains column headings. Default
%                                                value is false.
%                                ignoreWSpace  - true / false value indicating whether to ignore
%                                                leading and trailing whitespace in the column
%                                                headers; ignored if noHeader is set to true.
%                                                Default value is false.
%
% The parameters must be specified in the form of param-value pairs, parameter names are not
% case-sensitive and partial matching is supported.
%
% [C1 C2 C3] = CSVIMPORT( fileName, 'columns', {'C1', 'C2', C3'}, ... )
%   This form returns the data from columns in output variables C1, C2 and C3 respectively, the
%   column names are case-sensitive and must match a column name in the file exactly. When fetching
%   data in column mode the number of output columns must match the number of columns to read or it
%   must be one. In the latter case the data from the columns is returned as a single cell matrix.
%
% [C1 C2 C3] = CSVIMPORT( fileName, 'columns', [1, 3, 4], ,'noHeader', true, ... )
%   This form returns the data from columns in output variables C1, C2 and C3 respectively, the
%   columns parameter must contain the column indices when the 'noHeader' option is set to true.

%
% Notes:  1. Function has not been tested on badly formatted CSV files.
%         2. Created using R2007b but has been tested on R2006b.
%
% Revisions:
%   04/28/2009: Corrected typo in an error message
%               Added igonoreWSpace option
%   08/16/2010: Replaced calls to str2num with str2double, the former uses eval leading to unwanted
%               side effects if cells contain text with function names
%

if ( nargin == 0 ) || isempty( fileName )
  [fileName, filePath] = uigetfile( '*.csv', 'Select CSV file' );
  if isequal( fileName, 0 )
    return;
  end
  fileName = fullfile( filePath, fileName );
else
  if ~ischar( fileName )
    error( 'csvimport:FileNameError', 'The first argument to %s must be a valid .csv file', ...
      mfilename );
  end
end

%Setup default values
p.delimiter       = ',';
p.columns         = [];
p.outputAsChar    = false;
p.uniformOutput   = true;
p.noHeader        = false;
p.ignoreWSpace    = false;

validParams     = {     ...
  'delimiter',          ...
  'columns',            ...
  'outputAsChar',       ...
  'uniformOutput',      ...
  'noHeader',           ...
  'ignoreWSpace'        ...
  };

%Parse input arguments
if nargin > 1
  if mod( numel( varargin ), 2 ) ~= 0
    error( 'csvimport:InvalidInput', ['All input parameters after the fileName must be in the ' ...
      'form of param-value pairs'] );
  end
  params  = lower( varargin(1:2:end) );
  values  = varargin(2:2:end);

  if ~all( cellfun( @ischar, params ) )
    error( 'csvimport:InvalidInput', ['All input parameters after the fileName must be in the ' ...
      'form of param-value pairs'] );
  end

  lcValidParams   = lower( validParams );
  for ii =  1 : numel( params )
    result        = strmatch( params{ii}, lcValidParams );
    %If unknown param is entered ignore it
    if isempty( result )
      continue
    end
    %If we have multiple matches make sure we don't have a single unambiguous match before throwing
    %an error
    if numel( result ) > 1
      exresult    = strmatch( params{ii}, validParams, 'exact' );
      if ~isempty( exresult )
        result    = exresult;
      else
        %We have multiple possible matches, prompt user to provide an unambiguous match
        error( 'csvimport:InvalidInput', 'Cannot find unambiguous match for parameter ''%s''', ...
          varargin{ii*2-1} );
      end
    end
    result      = validParams{result};
    p.(result)  = values{ii};
  end
end

%Check value attributes
if isempty( p.delimiter ) || ~ischar( p.delimiter )
  error( 'csvimport:InvalidParamType', ['The ''delimiter'' parameter must be a non-empty ' ...
    'character array'] );
end
if isempty( p.noHeader ) || ~islogical( p.noHeader ) || ~isscalar( p.noHeader )
  error( 'csvimport:InvalidParamType', ['The ''noHeader'' parameter must be a non-empty ' ...
    'logical scalar'] );
end
if ~p.noHeader
  if ~isempty( p.columns )
    if ~ischar( p.columns ) && ~iscellstr( p.columns )
      error( 'csvimport:InvalidParamType', ['The ''columns'' parameter must be a character array ' ...
        'or a cell array of strings for CSV files containing column headers on the first line'] );
    end
    if p.ignoreWSpace
      p.columns = strtrim( p.columns );
    end
  end
else
  if ~isempty( p.columns ) && ~isnumeric( p.columns )
    error( 'csvimport:InvalidParamType', ['The ''columns'' parameter must be a numeric array ' ...
      'for CSV files containing column headers on the first line'] );
  end
end
if isempty( p.outputAsChar ) || ~islogical( p.outputAsChar ) || ~isscalar( p.outputAsChar )
  error( 'csvimport:InvalidParamType', ['The ''outputAsChar'' parameter must be a non-empty ' ...
    'logical scalar'] );
end
if isempty( p.uniformOutput ) || ~islogical( p.uniformOutput ) || ~isscalar( p.uniformOutput )
  error( 'csvimport:InvalidParamType', ['The ''uniformOutput'' parameter must be a non-empty ' ...
    'logical scalar'] );
end

%Open file
[fid, msg] = fopen( fileName, 'rt' );
if fid == -1
  error( 'csvimport:FileReadError', 'Failed to open ''%s'' for reading.\nError Message: %s', ...
    fileName, msg );
end

colMode         = ~isempty( p.columns );
if ischar( p.columns )
  p.columns     = cellstr( p.columns );
end
nHeaders        = numel( p.columns );

if colMode
  if ( nargout > 1 ) && ( nargout ~= nHeaders )
    error( 'csvimport:NumOutputs', ['The number of output arguments must be 1 or equal to the ' ...
      'number of column names when fetching data for specific columns'] );
  end
end

%Read first line and determine number of columns in data
rowData         = fgetl( fid );
rowData         = regexp( rowData, p.delimiter, 'split' );
nCols           = numel( rowData );

%Check whether all specified columns are present if used in column mode and store their indices
if colMode
  if ~p.noHeader
    if p.ignoreWSpace
      rowData     = strtrim( rowData );
    end
    colIdx        = zeros( 1, nHeaders );
    for ii = 1 : nHeaders
      result      = strmatch( p.columns{ii}, rowData );
      if isempty( result )
        fclose( fid );
        error( 'csvimport:UnknownHeader', ['Cannot locate column header ''%s'' in the file ' ...
          '''%s''. Column header names are case sensitive.'], p.columns{ii}, fileName );
      elseif numel( result ) > 1
        exresult  = strmatch( p.columns{ii}, rowData, 'exact' );
        if numel( exresult ) == 1
          result  = exresult;
        else
          warning( 'csvimport:MultipleHeaderMatches', ['Column header name ''%s'' matched ' ...
            'multiple columns in the file, only the first match (C:%d) will be used.'], ...
            p.columns{ii}, result(1) );
        end
      end
      colIdx(ii)  = result(1);
    end
  else
    colIdx        = p.columns(:);
    if max( colIdx ) > nCols
      fclose( fid );
      error( 'csvimport:BadIndex', ['The specified column index ''%d'' exceeds the number of ' ...
        'columns (%d) in the file'], max( colIdx ), nCols );
    end
  end
end

%Calculate number of lines
pos             = ftell( fid );
if pos == -1
  msg = ferror( fid );
  fclose( fid );
  error( 'csvimport:FileQueryError', 'FTELL on file ''%s'' failed.\nError Message: %s', ...
    fileName, msg );
end
data            = fread( fid );
nLines          = numel( find( data == sprintf( '\n' ) ) ) + 1;
%Reposition file position indicator to beginning of second line
if fseek( fid, pos, 'bof' ) ~= 0
  msg = ferror( fid );
  fclose( fid );
  error( 'csvimport:FileSeekError', 'FSEEK on file ''%s'' failed.\nError Message: %s', ...
    fileName, msg );
end

data            = cell( nLines, nCols );
data(1,:)       = rowData;
emptyRowsIdx    = [];
%Get data for remaining rows
for ii = 2 : nLines
  rowData       = fgetl( fid );
  if isempty( rowData )
    emptyRowsIdx = [emptyRowsIdx(:); ii];
    continue
  end
  rowData       = regexp( rowData, p.delimiter, 'split' );
  nDataElems    = numel( rowData );
  if nDataElems < nCols
    warning( 'csvimport:UnevenColumns', ['Number of data elements on line %d (%d) differs from ' ...
      'that on the first line (%d). Data in this line will be padded.'], ii, nDataElems, nCols );
    rowData(nDataElems+1:nCols) = {''};
  elseif nDataElems > nCols
    warning( 'csvimport:UnevenColumns', ['Number of data elements on line %d (%d) differs from ' ...
      'that one the first line (%d). Data in this line will be truncated.'], ii, nDataElems, nCols );
    rowData     = rowData(1:nCols);
  end
  data(ii,:)    = rowData;
end
%Close file handle
fclose( fid );
data(emptyRowsIdx,:)   = [];

%Process data for final output
uniformOutputPossible  = ~p.outputAsChar;
if p.noHeader
  startRowIdx          = 1;
else
  startRowIdx          = 2;
end
if ~colMode
  if ~p.outputAsChar
    %If we're not outputting the data as characters then try to convert each column to a number
    for ii = 1 : nCols
      colData     = cellfun( @str2double, data(startRowIdx:end,ii), 'UniformOutput', false );
      %If any row contains an entry that cannot be converted to a number then return the whole
      %column as a char array
      if ~any( cellfun( @isnan, colData ) )
        if ~p.noHeader
          data(:,ii)= cat( 1, data(1,ii), colData{:} );
        else
          data(:,ii)= colData;
        end
      end
    end
  end
  varargout{1}    = data;
else
  %In column mode get rid of the headers (if present)
  data            = data(startRowIdx:end,colIdx);
  if ~p.outputAsChar
    %If we're not outputting the data as characters then try to convert each column to a number
    for ii = 1 : nHeaders
      colData     = cellfun( @str2double, data(:,ii), 'UniformOutput', false );
      %If any row contains an entry that cannot be converted to a number then return the whole
      %column as a char array
      if ~any( cellfun( @isnan, colData ) )
        data(:,ii)= colData;
      else
        %If any column cannot be converted to a number then we cannot convert the output to an array
        %or matrix i.e. uniform output is not possible
        uniformOutputPossible = false;
      end
    end
  end
  if nargout == nHeaders
    %Loop through each column and convert to matrix if possible
    varargout = cell(nHeaders,1);
    for ii = 1 : nHeaders
      if p.uniformOutput && ~any( cellfun( @ischar, data(:,ii) ) )
        varargout{ii} = cell2mat( data(:,ii) );
      else
        varargout{ii} = data(:,ii);
      end
    end
  else
    %Convert entire table to matrix if possible
    if p.uniformOutput && uniformOutputPossible
      data        =  cell2mat( data );
    end
    varargout{1}  = data;
  end
end


function pos = plotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h

% Copyright 2010 Kelly Kearney

% Check input

if nargin < 1
    h = gca;
end

if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end

% Get position of axis in pixels

currunit = get(h, 'units');
set(h, 'units', 'pixels');
axisPos = get(h, 'Position');
set(h, 'Units', currunit);

% Calculate box position based axis limits and aspect ratios

darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    
    pos = axisPos;
    
else

    dx = diff(get(h, 'XLim'));
    dy = diff(get(h, 'YLim'));
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

% Convert plot box position to the units used by the axis

temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', get(h, 'parent'));
set(temp, 'Units', currunit);
pos = get(temp, 'position');
delete(temp);

function close_req(fig,~)

selection = questdlg('Save fiducials before closing?', ...
    'Close Request Function', ...
    'Yes','No','Yes'); 
switch selection 
    case 'Yes'
        save_button_fcn_close_req([],[],fig,[]);
        delete(gcf)
    case 'No'
        delete(gcf) 
end


