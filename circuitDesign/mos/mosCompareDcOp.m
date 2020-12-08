function result = mosCompareDcOp(circuit, designkit, simulator, ...
    verboseLevel, fid) 
%MOSCOMPAREDCOP comparison of operating point information computed in matlab 
% and with a circuit simulator
%
%  RESULT = mosCompareDcOp(CIRCUIT, DESIGNKIT, SIMULATOR, VERBOSELEVEL,
%  FID) returns a cell array in which each element is a cell array with three
%  elements, the first element being the name of a parameter of interest (see
%  description of the fields of a designkit, see cirDesignkitInit), the element
%  being the average and the maximum of the (unsigned) relative errors between
%  the value of the 
%  parameter under consideration computed in Matlab and with the given SIMULATOR,
%  averaged over all n-MOS and p-MOS transistors in the given CIRCUIT. 
%  For the time being only spectre is supported.
%  The VERBOSELEVEL is either 0, 1 or 2. For values of 1 and 2, some comparison
%  data is written into a given file identifier FID.
%
%
%  EXAMPLE:
%
%  a = mosCompareDcOp(circuit, designkit, 'spectre', 1, 1);
%
% Then a{1} can be an array that looks like
% a{1}
%
% ans = 
%
%  Columns 1 through 5
%
%    'idsnmos'    'vgsnmos'    'vdsnmos'    'vsbnmos'    'gmnmos'
%
%  Columns 6 through 10
%
%    'gmbsnmos'    'gdsnmos'    'vthnmos'    'vdsatnmos'    'ftnmos'
%
%  Columns 11 through 15
%
%    'idspmos'    'vgspmos'    'vdspmos'    'vsbpmos'    'gmpmos'
%
%  Columns 16 through 20
%
%    'gmbspmos'    'gdspmos'    'vthpmos'    'vdsatpmos'    'ftpmos'
%
% Further, a{2} is a cell array with the average rel. error values (averaged
% over all n-MOS transistors and all p-MOS transistors) and a{3} with the
% maximum value of the relative error.  
%
%  (c) IMEC, 2004
%  IMEC confidential 
%



debug = 0;

if not(ismember(verboseLevel, [0 1 2]))
  error('wrong value of the verbose level');
end

if not(strcmp(simulator, 'spectre'))
  error('simulator %s not supported for this function', simulator);
end

fprintf(fid, 'Comparison of operating point parameters of MOS transistors\n');
fprintf(fid, 'computed with matlab and %s for circuit %s\n', simulator, circuit.name); 
mosTypeArray = {'nmos', 'pmos'};
% result will only contain entries for parameters of interest that
% are different from an empty string. Therefore we define a running index,
% that indicates where we are in the filling of result:
index = 0;
for k = 1:length(mosTypeArray) % begin it. over mosTypeArray (index k)
  mosType = mosTypeArray{k};
  nMos = circuit.(strcat('n', mosType));
  paramsOfInterest = designkit.(simulator).paramsOfInterest.(mosType);
  nParamsOfInterest = length(paramsOfInterest);
  paramsOfInterestNames = designkit.(simulator).paramsOfInterestNames.(mosType);
  nParamsOfInterestNames = length(paramsOfInterest);
  relErrArray = zeros(nParamsOfInterest, nMos);
  for i = 1:nMos % begin it. over nMos (index i)
    mos = circuit.(strcat(mosType, 'List')){i};
    if verboseLevel > 1
      fprintf(fid, 'transistor %s:\n', circuit.(mos).name);
    end

    for j = 1:nParamsOfInterest % begin it. over nParamsOfInt. (index j)
      if designkit.(simulator).paramsOfInterest.(mosType){j} % if paramsOfInt.
	nameInMatlab = paramsOfInterest{j};
	nameInSimulator = designkit.(simulator).paramsOfInterestNames.(mosType){j};
	simulatorValue = circuit.(simulator).(mos).(nameInSimulator);
	matlabValue = circuit.(mos).(nameInMatlab);
	if simulatorValue == 0 | matlabValue == 0 % if (non)zero values
	  relErrArray(j, i) = NaN;
	  if verboseLevel > 1 % if verboseLevel > 1
	    fprintf(fid, '%s in %s: %s, in matlab: %s\n', nameInMatlab, ...
		simulator, eng(abs(simulatorValue)), eng(abs(matlabValue)));
	  end % end if verboseLevel > 1
	else
	  relErrArray(j, i) = (abs(matlabValue) - abs(simulatorValue)) / ...
	      abs(simulatorValue);
	  if verboseLevel > 1 % if verboseLevel > 1
	    fprintf(fid, '%s in %s: %s, in matlab: %s, rel. error = %g%%\n', ...
		nameInMatlab, simulator, eng(abs(simulatorValue)), ...
		eng(abs(matlabValue)) , 100*relErrArray(j, i)); 
	  end % end if verboseLevel > 1
	end % end if (non)zero values
      end % end if paramsOfInt.
    end % end it. over nParamsOfInt. (index j)
  end % end it. over nMos (index i)
  % construction of the return value result:
  for m = 1:nParamsOfInterest % begin it. over nParamsOfInterest (index m) 
    if paramsOfInterest{m} % if paramsOfInterest
      index = index+1;
      if debug
	fprintf(fid, 'm = %d, index = %d\n', m, index);
	fprintf(fid, 'paramsOfInterest{m} is %s\n', paramsOfInterest{m});
      end
      sum = 0;
      nTerms = 0;
      maxErr = 0;
      if debug
	fprintf(fid, 'nMos is %g, mosType is %s\n', nMos, mosType);
      end
      for n = 1:nMos % begin it. over nMos (index n)
	if debug
	  fprintf(fid, 'relErrArray(%d, %d) = %g\n', m, n, relErrArray(m,n));
	end
	if not(isnan(relErrArray(m,n)))
	  sum = sum + abs(relErrArray(m,n));
	  nTerms = nTerms + 1;
	  maxErr = max(abs(relErrArray(m,n)), maxErr);
	end
      end  % end it. over nMos (index n)    
      result{1}{index} = strcat(paramsOfInterest{m}, mosType);
      result{2}{index} = sum/nTerms;
      result{3}{index} = maxErr;
      if verboseLevel > 0 % if verboseLevel > 0
	fprintf(fid, 'error on %s (%s): average = %g%%, maximum = %g%%\n', ...
	    paramsOfInterest{m}, mosType, 100*result{2}{index}, ...
	    100*result{3}{index});
      end % end if verboseLevel > 0
    end % end if paramsOfInterest
  end % end it. over nParamsOfInterest (index m)
end % end it. over mosTypeArray (index k)
if verboseLevel > 0
  fprintf(fid, '**********************************\n');
end

