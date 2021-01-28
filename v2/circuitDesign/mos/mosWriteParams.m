function parStruct = mosWriteParams(mos, parStruct)
%MOSWRITEPARAMS writes parameters that have been determined in Matlab to a
%  the fields of a structure parStruct
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


for j = 1:length(mos.paramFields)
  parStruct.(mos.name).(mos.paramPrintNames{j}) = mos.(mos.paramFields{j}); 
end
