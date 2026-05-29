function outsideBounds(fig)
    h = guidata(fig);
    %h.ax, 'Units', 'normalized', 'Position', [0.0272 0.0576 0.55 0.853]);
    % ax = get(fig, 'CurrentAxes');
    % if isempty(ax), return; end
    % 
    % cp = get(ax, 'CurrentPoint');
    % x = round(cp(1,1));
    % y = round(cp(1,2));
    % 
    % [rows, cols, ~] = size(h.stack);
    % 
    % if x < 1 || x > cols || y < 1 || y > rows
    %     fig.WindowButtonMotionFcn = @(f,~)updatePixelIntensity(f);
    % end
    
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