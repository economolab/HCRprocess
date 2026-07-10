% return parameters

function out = params(in)
    
    d = dictionary();
    
    d('abcDir') = 'D:\allen_brain_atlas';
    d('fijiScriptsDir') = 'C:\Users\ECONOM~1\DOCUME~1\FIJI-W~1\Fiji.app\scripts';
    d('HCRprocessEnv') = 'HCRprocess';

    out = d(in);
   
end