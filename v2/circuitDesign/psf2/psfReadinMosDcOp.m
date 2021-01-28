function circuit = psfReadinMosDcOp(circuit, designkit, psf, ...
    oppointAnalysisName)
%PSFREADINMOSDCOP reading DC operating point information from Spectre
% into the datastructure of the circuit.
%
% CIRCUIT = PSFREADINMOSDCOP(CIRCUIT, DESIGNKIT, PSF, OPPOINTANALYSISNAME)
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%



for i = 1:circuit.nnmos
  mos = circuit.nmosList{i};
  circuit = psfReadinMosDcOpAux(circuit, ...
      designkit.spectre.paramsOfInterestNames.nmos, mos, psf, ...
    oppointAnalysisName, designkit.psf.DcExtraString.nmos);
end
for i = 1:circuit.npmos
  mos = circuit.pmosList{i};
  circuit = psfReadinMosDcOpAux(circuit, ...
      designkit.spectre.paramsOfInterestNames.pmos, mos, psf, ...
    oppointAnalysisName, designkit.psf.DcExtraString.pmos);
end


function circuit = psfReadinMosDcOpAux(circuit, paramsOfInterestNames, ...
    mos, psf, oppointAnalysisName, extraString) 
for j = 1:length(paramsOfInterestNames)
  circuit.spectre.(mos).(paramsOfInterestNames{j}) = psfdcop(psf, ...
      oppointAnalysisName, strcat(mos, extraString), paramsOfInterestNames{j});
end

