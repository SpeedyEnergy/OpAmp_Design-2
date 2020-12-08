function mosPrintCaps(fid, circuit)
% MOSPRINTCAPS prints parasitic capacitors of MOS transistors in a circuit and
% its subcircuits.
%
%   mosPrintCaps(FID, CIRCUIT) prints the parasitic capacitors of the
%   transistors in the given CIRCUIT and its subcircuits to the given file
%   identifier FID. Two file identifiers are automatically available and need
%   not be opened.  They are FID=1 (standard output) and FID=2 (standard
%   error). 
%
%  See also fopen.
%
%
%  (c) IMEC, 2005
%  IMEC confidential 
%

DEBUG = 0;

if (circuit.nnmos > 0) | (circuit.npmos > 0) 
  fprintf(fid, '\n**************************************\n');
  fprintf(fid, 'parasitic capacitors in circuit %s\n', circuit.name); 
  fprintf(fid, '****************************************\n\n');
  
  % print transistor sizes:
  capList = {'cgs', 'cgd', 'cgb', 'csb', 'cdb'};
  mos = {'nmos', 'pmos'};
  for k = 1:length(mos)
    for i = 1:circuit.(strcat('n', mos{k}))
      mosList = strcat(mos{k}, 'List');
      fprintf(fid, '%s: ', circuit.(circuit.(mosList){i}).name)
      fprintf(fid, 'W=%sm L=%sm ', ...
	  eng(circuit.(circuit.(mosList){i}).w, 3), ...
	  eng(circuit.(circuit.(mosList){i}).lg, 3));
      for j = 1:length(capList)
	if DEBUG
	  fprintf(1, 'capList{%d} is %s\n', j, capList{j});
	end
	if isfield(circuit.(circuit.(mosList){i}), capList{j})
	  string.(capList{j}) = ...
	      eng(circuit.(circuit.(mosList){i}).(capList{j}), 3);
	else
	  string.(capList{j}) = 'xxx';
	end
	fprintf(fid, '%s=%sF ', capList{j}, string.(capList{j}));
      end
      fprintf(fid, '\n');
    end
  end
end

  
for i = 1:circuit.nsubckt
  mosPrintCaps(fid, circuit.(circuit.subcktList{i}));
end
  
