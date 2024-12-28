function getLipoParams(fig)

h = guidata(fig);

h.lipo.bg = zeros(h.Nchans, 1);
h.lipo.amp = zeros(h.Nchans, 1);


crossCoord = h.lipoCrosshair.Position;
x = round(crossCoord(1));
y = round(crossCoord(2));

bgROI = createMask(h.bgROI);

for i = 1:h.Nchans
    tmp = h.stack(:,:, i);
    h.lipo.bg(i) = median(tmp(bgROI));
    h.lipo.amp(i) = h.stack(y, x, i)-h.lipo.bg(i);
    
end
h.lipo.amp(h.lipo.amp<0) = 0;
guidata(fig, h);