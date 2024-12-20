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