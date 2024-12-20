function updateImage(~, ~, fig)
h = guidata(fig);

im = zeros(h.ymax, h.xmax, 3);

for i = 1:h.Nchans
    if ~get(h.enable_checkbox(i), 'Value')
        continue;
    end
    
    clr = get(h.color_button(i), 'BackgroundColor');
    scim = getScaledImage(h, i);
    
    for j = 1:3
        im(:,:,j) = im(:,:,j)+scim.*clr(j);
    end

        
end



if get(h.hilo_checkbox, 'Value')

    mask = ones(size(im,1), size(im,2));
    for j = 1:3
        mask = mask & (im(:,:,j) == 0);
    end

    im = imoverlay(im, mask, [0 0 1]);
end



tstr = h.fn;
tstr = strrep(tstr, '_', ' ');
set(h.vis, 'CData', im);
t = get(h.ax, 'Title');
set(t, 'String', tstr)
drawnow;
