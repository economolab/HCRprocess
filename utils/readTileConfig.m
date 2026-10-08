function [names, series, coordStr, dim] = readTileConfig(path)

    lines = readlines(path);
    lines = strtrim(lines);
    lines = lines(lines ~= "" & ~startsWith(lines, "#"));

    dim = 3;
    isDim = startsWith(lines, "dim");
    if any(isDim)
        dim = str2double(extractAfter(lines(find(isDim,1)), "="));
        lines = lines(~isDim);
    end

    tok = regexp(lines, '^(.*?);\s*(\d*)\s*;\s*\(([^)]*)\)\s*$', 'tokens', 'once');
    n = numel(tok);
    names    = strings(n,1);
    series   = strings(n,1);   % kept as text too; empty string if absent
    coordStr = strings(n,1);   % raw "x, y, z" text, untouched
    
    for i = 1:n
        t = tok{i};
        names(i)    = strtrim(t{1});
        series(i)   = strtrim(t{2});
        coordStr(i) = strtrim(t{3});
    end

end