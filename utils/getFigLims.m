function lims = getFigLims(fig, num, im)

h = guidata(fig);

if (get(h.autolim(num), 'Value')==0)
    cmin = str2double(get(h.cmin_edit(num), 'String'));
    cmax = str2double(get(h.cmax_edit(num), 'String'));
else
    a = qprctile(double(im(:)), [0.1 99.995]);

    cmin = a(1);
    cmax = a(2);
    
    if cmax<=cmin
        cmax = cmin+1;
    end
    
    set(h.cmin_edit(num), 'String', num2str(round(cmin)));
    set(h.cmax_edit(num), 'String', num2str(round(cmax)));
end

guidata(fig, h);
lims = [cmin cmax];