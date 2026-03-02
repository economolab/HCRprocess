function grid = makeGrid(im,rects,grid_line_thickness)

    arguments
        im
        rects
        grid_line_thickness = 3
    end
    
    half_thickness = round(grid_line_thickness / 2);
    im_height = size(im,1);
    im_width = size(im,2);
    grid = zeros(im_height,im_width,"logical");
    rect_width = im_width / rects;
    rect_height = im_height / rects;

    for i = 1:rects-1
        x1 = uint16(i*rect_height - half_thickness);
        x2 = uint16(i*rect_height + half_thickness);
        y1 = uint16(i*rect_width - half_thickness);
        y2 = uint16(i*rect_width + half_thickness);
        grid(x1:x2,:) = 1;
        grid(:,y1:y2) = 1;
    end

end