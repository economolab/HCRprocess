% strs = cell array of strings to be filtered 
% kw = string array of keywords to look for in strs
% mode = 'keep' or 'remove' depending on desired filtering
% removes or keeps strings in the input cell array that contain any of the
% input keywords

function strs_filt = filt_kw(strs, kw, mode)

    % which strings contain the keyword
    mask = contains(strs,kw,'IgnoreCase',true);
    
    % filter appropriately depending on desired filtering mode 
    switch mode
        case 'keep'
            strs_filt = strs(mask);
        case 'remove'
            strs_filt = strs(~mask);
    end
   
end