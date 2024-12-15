% input = cell array of strings and string array of keywords
% filters out strings in the input cell array that contain any of the
% input keywords

function strs_filt = filt_kw(strs, kw)

    % strings to ignore
    mask = contains(strs,kw,'IgnoreCase',true);
    
    % keep the strings that don't have any of the keywords
    strs_filt = strs(~mask);
   
end