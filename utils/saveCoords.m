function saveCoords(fig)
    if ~strcmp(get(fig, 'SelectionType'), 'extend')
        return;
    end

    h = guidata(fig);
    
    ax = h.ax;

    cp = get(ax, 'CurrentPoint');
    x  = round(cp(1,1));
    y  = round(cp(1,2));

    [rows, cols, ~] = size(h.stack);

    if x < 1 || x > cols || y < 1 || y > rows
        return;
    end

    h.lastClick = [x, y];

    % delete any previous selection box
    if isfield(h, 'selectionBox') && isvalid(h.selectionBox)
        delete(h.selectionBox);
    end

    %clear current button press functions
    % fig.WindowButtonDownFcn = [];
    % fig.WindowButtonUpFcn = [];
    % fig.WindowButtonMotionFcn = [];

    % h.selectionBox = drawpolygon(ax);
    % 
    % wait(h.selectionBox)

    % create new rectangle at click point with size 1x1 to start
    h.selectionBox = rectangle(h.ax, ...
        'Position',  [x, y, 1, 1], ...
        'EdgeColor', [0.6 0.6 0.6], ...
        'LineWidth', 1.5, ...
        'LineStyle', '--');


    guidata(fig, h)

    fig.WindowButtonMotionFcn = @(f,~) updateSelectionBox(f);
end