function pctval = qprctile(v, p)

pctval = zeros(length(p),1);
p=p./100;
sortv = sort(v, 'ascend');

for i = 1:length(p)
    len = length(v);

    val = len*p(i);
    ind1 = floor(len*p(i));
    ind2 = ind1+1;

    if (ind1 == len)
        pctval(i) = double(sortv(end));
    elseif (ind1<1)
        pctval(i) = double(sortv(1));
    else
        pctval(i) = double(sortv(ind1)) + (double(sortv(ind2))-double(sortv(ind1))).*(val-ind1);
    end
end