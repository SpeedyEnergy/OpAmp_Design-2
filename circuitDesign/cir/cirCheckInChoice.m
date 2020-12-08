function circuit = cirCheckInChoice(circuit, choice)
%CIRCHECKINCHOICE defines every field of a structure of design choices as a ...
%    variable in the base workspace and copies this structure to a field of 
%    the given circuit
%
%   CIRCUIT = CIRCHECKINCHOICE(CIRCUIT, CHOICE) adds the field CHOICE to the
%   given CIRCUIT. These fields have the same contents as in the variable
%   CHOICE. Further, each field of CHOICE is assigned to new variables in the
%   base workspace. The name of each new variable is the name of the
%   corresponding field.
%
%   Examples:
%
%   Assume there is a structure "choice" in the base workspace, with
%     choice.Vdd.val = 1;
%     choice.Mn1.vgs = 0.5;
%   while there exists a circuit "ota".
%   Then a call like
%
%   ota = cirCheckInChoice(ota, choice)
%
%   will have the following effects:
%   1. the variable "ota" has a new field "choice" and
%      ota.choice.Vdd.val = 1;
%      ota.choice.Mn1.vgs = 0.5;
%   2. new variables are created, namely
%      Vdd.val 
%      Mn1.vgs
%      Their value is 1 and 0.5, respectively.
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

debug = 0;
if debug
  fprintf(1, '2nd argument of cirCheckInChoice is %s\n', inputname(2));
end

circuit.choice = choice;
checkInStructure(inputname(2), choice, '');


% inline function checkInStructure:
% 
function checkInStructure(main, structure, structureString)
%Assignment of all fields of a structure into the base workspace
%
%  CHECKINSTRUCTURE(STRUCTURE) assigns each field of the given STRUCTURE to new
%  variables in the base workspace. The name of each new variable is the name
%  of the corresponding field. 
%
%  Examples:
%
%  Assume there is a structure "choice" in the base workspace, with
%  choice.a = 7;
%  choice.b = 'string';
%
%  Then the call in the base workspace
%
%  checkInstrStructure(choice)
%
% will define the variables "a" and "b" and "a" has the value 7, while "b" has
% the value 'string';
%

debug = 0;


fieldNames = fieldnames(structure);
for i = 1:length(fieldNames)
  if debug
    fprintf(1, 'fieldNames{%d} is %s\n', i, fieldNames{i});
    fprintf(1, '1st argument of checkInStructure is\n');
    disp(main);
    fprintf(1, '\n2nd argument of checkInStructure is\n');
    disp(structure);
    fprintf(1, '\n3rd argument of checkInStructure is\n');
    disp(structureString);
  end
  if isstruct(structure.(fieldNames{i}))
    if structureString
      checkInStructure(main, structure.(fieldNames{i}), strcat(structureString, '.', ...
	  fieldNames{i}));
    else
      checkInStructure(main, structure.(fieldNames{i}), fieldNames{i});
    end 
  else
    if structureString
      newVarName = strcat(structureString, '.', fieldNames{i});
    else
      newVarName = fieldNames{i};
    end
    if debug
      fprintf(1, 'newVarName = %s\n', newVarName);
      fprintf(1, 'command to execute is\n%s\n', sprintf('%s = %s', newVarName, ...
	  strcat(main, '.', newVarName)));
    end
    evalin('base', sprintf('%s = %s;', newVarName, strcat(main, '.', newVarName)));
  end
end

