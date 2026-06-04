function getAvgIntensity(fig)
    
    h = guidata(fig);

    ax = h.ax;
    lastx = h.lastClick(1);
    lasty = h.lastClick(2);
    
    cp = get(ax, 'CurrentPoint');
    x  = round(cp(1,1));
    y  = round(cp(1,2));
    
    xvals = [lastx, x];
    yvals = [lasty, y];

    [rows, cols, ~] = size(h.stack);

    if x < 1 || x > cols || y < 1 || y > rows
        return;
    end

    visChans = zeros(1, h.Nchans);
    
    for i = 1:h.Nchans
        uibox = h.enable_checkbox(i);
        if uibox.Value == 1
            visChans(i) = 1;
        end
    end

    sumvis = sum(visChans);
    str = strings(1, sumvis+1);
    str(1) = sprintf('(%d, %d)', x, y);
    
    cnt = 2;
    for i = 1:h.Nchans
        if visChans(i) == 1
            vals = h.stack(min(yvals):max(yvals), min(xvals):max(xvals), i);
            meanIntensity = round(mean(vals, 'all'));
            str(cnt) = sprintf('  %s average: %d', h.lambdas(i), meanIntensity);
            cnt = cnt + 1;
        end
    end

    set(h.pixelInfo, 'String', str);

    %set motion callback back to updatePixelIntensity, delete box on button release
    fig.WindowButtonMotionFcn = @(f,~)outsideBounds(f);
    if isfield(h, 'selectionBox') && isvalid(h.selectionBox)
        delete(h.selectionBox);
    end
    
    guidata(fig, h)

    
end