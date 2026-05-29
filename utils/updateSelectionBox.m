function updateSelectionBox(fig)

    h = guidata(fig);

    if ~isfield(h, 'lastClick') || ~isfield(h, 'selectionBox') || ~isvalid(h.selectionBox)
        return;
    end

    [rows, cols, ~] = size(h.stack);

    cp = get(h.ax, 'CurrentPoint');
    x  = round(cp(1,1));
    y  = round(cp(1,2));

    % clamp to image bounds
    x = max(1, min(cols, x));
    y = max(1, min(rows, y));

    lastx = h.lastClick(1);
    lasty = h.lastClick(2);

    xmin  = min(lastx, x);
    ymin  = min(lasty, y);
    rectw = max(1, abs(x - lastx));   % rectangle needs positive width/height
    recth = max(1, abs(y - lasty));

    set(h.selectionBox, 'Position', [xmin, ymin, rectw, recth]);
end