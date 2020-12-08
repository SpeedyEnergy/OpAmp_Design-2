function mosPrintSizesAndOpInfo(fid, circuit)
% MOSPRINTSIZESANDOPINFO print sizes and operating point information of MOS
% transistors in a circuit and its subcircuits.
%
%   mosPrintSizesAndOpInfo(FID, CIRCUIT) prints the sizes and operating point
%   information of the transistors in the given CIRCUIT and its subcircuits to
%   the given file identifier FID. Two file identifiers are automatically
%   available and need not be opened.  They are FID=1 (standard output) and
%   FID=2 (standard error).
%
%  See also fopen.
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


if (circuit.nnmos > 0) | (circuit.npmos > 0) 
  fprintf(fid, '\n***********************\n');
  fprintf(fid, 'circuit %s\n', circuit.name); 
  fprintf(fid, '***********************\n\n');
  
  % print transistor sizes:
  fprintf(fid, 'Transistor sizes:\n');
  for i = 1:circuit.nnmos
    fprintf(fid, '%s: ', circuit.(circuit.nmosList{i}).name)
    fprintf(fid, 'W=%sm\tL=%sm\tnFingers=%d\tfingerWidth=%s\n', ...
	eng(circuit.(circuit.nmosList{i}).w, [], 2), ...
	eng(circuit.(circuit.nmosList{i}).lg), ...
	circuit.(circuit.nmosList{i}).nFingers, ...
	eng(circuit.(circuit.nmosList{i}).w / ...
	circuit.(circuit.nmosList{i}).nFingers, [], 2)); 
  end
  for i = 1:circuit.npmos
    fprintf(fid, '%s: W=%sm\tL=%sm\tnFingers=%d\tfingerWidth=%s\n', ...
	circuit.(circuit.pmosList{i}).name, ...
	eng(circuit.(circuit.pmosList{i}).w, [], 2), ...
	eng(circuit.(circuit.pmosList{i}).lg), ...
	circuit.(circuit.pmosList{i}).nFingers, ...
	eng(circuit.(circuit.pmosList{i}).w / ...
	circuit.(circuit.pmosList{i}).nFingers, [], 2)); 
  end
  
  % Print DC parameters:
  fprintf(fid, '\nOperating point parameters:\n');

  % first the n-MOS transistors:
  for i = 1:circuit.nnmos
    if isfield(circuit.(circuit.nmosList{i}), 'gm')
      gmString = eng(circuit.(circuit.nmosList{i}).gm,[],2);
    else
      gmString = 'xxx';
    end
    if isfield(circuit.(circuit.nmosList{i}), 'gds')
      gdsString = eng(circuit.(circuit.nmosList{i}).gds,[],2);
    else
      gdsString = 'xxx';
    end
    if isfield(circuit.(circuit.nmosList{i}), 'ids')
      idsString = eng(circuit.(circuit.nmosList{i}).ids,[],2);
    else
      idsString = 'xxx';
    end
    if isfield(circuit.(circuit.nmosList{i}), 'vgs') & ...
	  isfield(circuit.(circuit.nmosList{i}), 'vth') 
      vovString = eng(circuit.(circuit.nmosList{i}).vgs - ...
	circuit.(circuit.nmosList{i}).vth,[],2); 
    else
      vovString = 'xxx';
    end
    if isfield(circuit.(circuit.nmosList{i}), 'vdsat')
      vdsatString = eng(circuit.(circuit.nmosList{i}).vdsat,[],2);
    else
      vdsatString = 'xxx';
    end
    if isfield(circuit.(circuit.nmosList{i}), 'vds')
      vdsString = eng(circuit.(circuit.nmosList{i}).vds,[],2);
    else
      vdsString = 'xxx';
    end
    
    fprintf(fid, '%s: gm=%sS  gds=%sS ids=%sA  vov=%sV  vdsat=%sV  vds=%sV\n', ...
	circuit.(circuit.nmosList{i}).name, gmString, gdsString, idsString, vovString, ...
	vdsatString, vdsString); 
  end

  % same for p-MOS transistors:
  for i = 1:circuit.npmos
    if isfield(circuit.(circuit.pmosList{i}), 'gm')
      gmString = eng(circuit.(circuit.pmosList{i}).gm,[],2);
    else
      gmString = 'xxx';
    end
    if isfield(circuit.(circuit.pmosList{i}), 'gds')
      gdsString = eng(circuit.(circuit.pmosList{i}).gds,[],2);
    else
      gdsString = 'xxx';
    end
    if isfield(circuit.(circuit.pmosList{i}), 'ids')
      idsString = eng(circuit.(circuit.pmosList{i}).ids,[],2);
    else
      idsString = 'xxx';
    end
    if isfield(circuit.(circuit.pmosList{i}), 'vgs') & ...
	  isfield(circuit.(circuit.pmosList{i}), 'vth') 
      vovString = eng(circuit.(circuit.pmosList{i}).vgs - ...
	circuit.(circuit.pmosList{i}).vth,[],2); 
    else
      vovString = 'xxx';
    end
    if isfield(circuit.(circuit.pmosList{i}), 'vdsat')
      vdsatString = eng(circuit.(circuit.pmosList{i}).vdsat,[],2);
    else
      vdsatString = 'xxx';
    end
    if isfield(circuit.(circuit.pmosList{i}), 'vds')
      vdsString = eng(circuit.(circuit.pmosList{i}).vds,[],2);
    else
      vdsString = 'xxx';
    end
    
    fprintf(fid, '%s: gm=%sS  gds=%sS ids=%sA  vov=%sV  vdsat=%sV  vds=%sV\n', ...
	circuit.(circuit.pmosList{i}).name, gmString, gdsString, idsString, vovString, ...
	vdsatString, vdsString); 
  end
end
for i = 1:circuit.nsubckt
  mosPrintSizesAndOpInfo(fid, circuit.(circuit.subcktList{i}));
end
  
