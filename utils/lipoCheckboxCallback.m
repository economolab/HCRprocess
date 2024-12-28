function lipoCheckboxCallback(~, ~, fig)

h = guidata(fig);
enable = get(h.crosshairs_checkbox  , 'Value');
if enable
    h.bgROI.Visible = 'On';
    h.lipoCrosshair.Visible = 'On';
else
    h.bgROI.Visible = 'Off';
    h.lipoCrosshair.Visible = 'Off';
end