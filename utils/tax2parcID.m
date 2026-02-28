% take an allen ccf taxonomy label (presumed at the structure level, but 
% can be set to something different) and convert it to its associated
% parcellation index. Works with single or multiple taxonomy labels
function parcID = tax2parcID(tax,tax_level)

    arguments
        tax
        tax_level = 'structure'
    end

    abcDir = params('abcDir');
    colorDir = fullfile(abcDir,'metadata','Allen-CCF-2020','20230630','views');
    ccfTerms = fullfile(colorDir,'parcellation_to_parcellation_term_membership_acronym.csv');
    Tterms = readtable(ccfTerms);
    
    parcID = zeros(length(tax),1);
    for i=1:length(tax)
        mask = strcmp(Tterms.(tax_level),tax{i});
        parcID(i) = Tterms.parcellation_index(mask); 
    end

end