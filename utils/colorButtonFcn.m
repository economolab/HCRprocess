function colorButtonFcn(~, ~, fig, chan)

h = guidata(fig);

newclr = uisetcolor(h.color_button(chan));
set(h.color_button(chan), 'BackgroundColor', newclr);
updateImage([], [], h.fig);