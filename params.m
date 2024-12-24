% return parameters

function out = params(in)
    
    d = dictionary();

    d('fijiScriptsDir') = 'C:\Users\jpv88\Fiji.app\scripts';
    d('HCRprocessEnv') = 'HCRprocess';

    out = d(in);
   
end