function outsideBounds(fig)
    h = guidata(fig);
    
    cp = get(h.fig, 'CurrentPoint');
    x = round(cp(1,1));
    y = round(cp(1,2));
    
    h.ax.Units = 'pixels';
    axPos = h.ax.Position;
    h.ax.Units = 'normalized';
    bottom = axPos(1,1);
    left = axPos(1,2);
    width = axPos(1,3);
    height = axPos(1,4);
    
    if x<left || x>left+width || y<bottom || y>bottom+height
        fig.WindowButtonMotionFcn = @(f,~)updatePixelIntensity(f);
    end