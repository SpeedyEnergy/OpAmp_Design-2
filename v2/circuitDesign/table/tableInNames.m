function result = tableInNames(table)
%TABLEINNAMES input variables for the table of a MOS transistor
%
%   CELLARRAY = tableInNames(TABLE) returns an array of strings that 
%   are the names of the input variables of a given TABLE of intrinsic 
%   operating point parameters of a MOS
%   transistor. The order of these strings also corresponds to the order of the
%   indices in the table values. E.g. when 'vgs' is the second element in the
%   resulting CELLARRAY, then this also means that the second index of the
%   table corresponds to vgs.
%
%   EXAMPLE:
%
%     stringArray = tableInNames(N);
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

result = table.Info.inputs;

