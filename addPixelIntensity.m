function addPixelIntensity(fig)
%add pixel info bar showing 12bit pixel intensity for each channel on mouse
%click and drag... currently interferes a bit with the built-in zoom
%function but ill deal with that later
    h = guidata(fig);
    
    title = 'Click scroll wheel and drag for average intensity';
    %create info bar
    h.pixelInfo = uicontrol(h.fig, ...
        'Style',               'text', ...
        'String',   title, ...
        'Units',               'normalized', ...
        'Position',            [0.53 0.65 0.05 0.26], ...      
        'BackgroundColor',     [0.9 0.9 0.9], ...
        'ForegroundColor',     [0 0 0], ...
        'HorizontalAlignment', 'left', ...
        'FontSize',            9);
    
    %get channel wavelengths
    sp = split(h.fn, '__');
    chansstr = sp(2);
    chans = split(chansstr, '_');
    h.lambdas = strings(1,h.Nchans);
    for i=1:h.Nchans
        chan = chans(i);
        cstr = chan{1};
        lambda = cstr(end-2:end);
        h.lambdas(1,i) = lambda;
    end
    
    guidata(fig, h)
    
    

    % fig.WindowButtonDownFcn = @(f,~)saveCoords(f);
    % fig.WindowButtonUpFcn = @(f,~)getAvgIntensity(f);
    fig.WindowButtonMotionFcn = @(f,~)updatePixelIntensity(f);

end