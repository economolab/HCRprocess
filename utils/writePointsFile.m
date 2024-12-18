function writePointsFile(fn,points)
% Open file
fid = fopen(fn, 'w');

if fid == -1
    error('Cannot open file: %s', filename);
end

try
    % Write header
    fprintf(fid, 'point\n');
    fprintf(fid, '%d\n', size(points, 1));

    % Write data
    for ii = 1:size(points, 1)
        fprintf(fid, '%f %f %f\n', points(ii, 1), points(ii, 2), points(ii, 3));
    end

catch ME
    % Ensure the file is closed even if an error occurred
    fclose(fid);
    rethrow(ME);
end

% Close file
fclose(fid);