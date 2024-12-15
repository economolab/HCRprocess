% correctly read a table that includes row names

function T = read_table_with_row_names(path)

    T = readtable(path,'VariableNamingRule','preserve');
    T.Properties.RowNames = T{:,'Row'};
    T = removevars(T,'Row');
    
end