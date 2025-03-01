% correctly read a table that includes row names

function T = read_table_with_row_names(path)

    T = readtable(path,'VariableNamingRule','preserve','Delimiter',',');
    T.Properties.RowNames = T{:,'Row'};
    T = removevars(T,'Row');
    T.Properties.DimensionNames{1} = 'Row';
    
end