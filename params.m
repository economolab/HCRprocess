% return parameters

function out = params(in)
    
    d = dictionary();
    
    d('abcDir') = 'C:\Users\jpv88\Documents\allen-brain-cell-atlas';
    d('fijiScriptsDir') = 'C:\Users\jpv88\Fiji.app\scripts';
    d('HCRprocessEnv') = 'HCRprocess';

    out = d(in);
   
end