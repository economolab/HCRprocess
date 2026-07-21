function updatePixelIntensity(fig)
%get intensity for each pixel on mouseover

    h = guidata(fig);
    if ~isfield(h, 'pixelInfo') || ~isvalid(h.pixelInfo)
        return;
    end
  
    warning off MATLAB:modes:mode:InvalidPropertySet

    ax = get(fig, 'CurrentAxes');
    if isempty(ax), return; end

    %get cursor position
    cp = get(ax, 'CurrentPoint');
    x  = round(cp(1,1));
    y  = round(cp(1,2));

    [rows, cols, nChan] = size(h.stack);

    %do nothing if cursor not in image bounds
    if x < 1 || x > cols || y < 1 || y > rows
        return;
    end

    %get channel wavelengths
    sp = split(h.fn, '__');
    chansstr = sp(2);
    chans = split(chansstr, '_');
    h.lambdas = strings(1,h.Nchans);
    for i=1:nChan
        chan = chans(i);
        cstr = chan{1};
        lambda = cstr(end-2:end);
        h.lambdas(1,i) = lambda;
    end

    str = strings(1,6);
    str(1) = sprintf('(%d, %d)', x, y);

    try
        for i=1:nChan
            val = round(h.stack(y+1, x+1, i));
            str(i+1) = sprintf('  %s intensity: %d', h.lambdas(i), val);
        end
    catch
    end

    set(h.pixelInfo, 'String', str);

    fig.WindowButtonDownFcn = @(f,~)saveCoords(f);
    fig.WindowButtonUpFcn = @(f,~)getAvgIntensity(f);
    guidata(fig, h)
end