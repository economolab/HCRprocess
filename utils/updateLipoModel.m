function updateLipoModel(~, ~, fig)

getLipoParams(fig);
h = guidata(fig);

hold (h.lipoAx, 'Off');
plot(h.lipoAx, h.lipo.bg);
hold (h.lipoAx, 'On');
plot(h.lipoAx, h.lipo.amp);