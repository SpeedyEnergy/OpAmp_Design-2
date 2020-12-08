function result = isNtype(table)
%ISNTYPE returns 1 if the given TABLE corresponds to an n-type transistor
%

result = strcmpi(tableType(table), 'n');
