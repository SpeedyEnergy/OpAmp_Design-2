function cirMosCheckSaturation(circuit, fid)
% checking all MOS transistors for saturation
%
% CIRMOSCHECKSATURATION(CIRCUIT, FID) checks all n-MOS and p-MOS transistors of
% a circuit whether or not they are saturated. The function writes a warning to
% the file identifier FID for every transistor that is not saturated
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


for i = 1:circuit.nnmos
  mosCheckSaturation(circuit.(circuit.nmosList{i}));
end
for i = 1:circuit.npmos
  mosCheckSaturation(circuit.(circuit.pmosList{i}));
end
