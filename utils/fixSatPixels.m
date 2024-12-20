% fix saturated pixels after registration

function fixSatPixels(fname)

    [V, info] = read_tiff(fname);
    
    sz = size(V);
    for i=1:sz(3)
        for j=1:sz(4)
            section = V(:,:,i,j);
            mask = (section > 60000);
            section(mask) = 0;
            V(:,:,i,j) = section;
        end
    end

    write_tiff(fname, V, info)

end