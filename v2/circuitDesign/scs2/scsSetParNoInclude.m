function  [OrigPAR, status] = scsSetParNoInclude(scsName,PAR, processIncludeFiles)
% [OrigPAR,status] = SCSSETPAR(scsName,PAR)
%     Script that
%       - reads in the spectre input file scsName.scs
%       - updates the parameter statements with the number present in the
%         structure PAR
%       - PAR can contain:
%           - values for parameter statements (numbers or strings)
%               (use "no change" if value may not be changed)
%               (use "remove" if parameter must be deleted from netlist)
%               (no parameters will be added, so all fields must already be
%                 defined in original input file)
%           - device parameters
%               (here device parameters can be added, might be tricky but
%                 it's also usefull.)
%     Outputs:
%       - OrigPAR = original parameter strings
%       - status = number of replaced parameters
%            -1 -> error has occured
%
% JCx, 06-01-2004

% JCx, 10-10-2003
%   First public version
% JCx, 05-11-2003
%   Added the possiblity to change the temperature also.  'temp' cannot be
%   a parameter, it must be specified as an option.
% JCx, 16-12-2003
%   Bug corrections
% JCx. 23-12-2003
%   Added technology dependence with global variable SCS_TECHNOLOGY
% JCx, 06-01-2004
%   Bug fix
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%




% initialize
debug = 0;
OrigPAR = [];
status = 0;
if nargin < 2
    PAR = [];
end
newline = sprintf('\n');
nextline1 = sprintf('\n+');
nextline2 = sprintf('\\\n');
tab = sprintf('\t');
blank = ' ';
% . . technology setting
global SCS_TECHNOLOGY
SCS_TECHNOLOGY_LIST = {'st_bicmos6g','dummy'};  % only these are implemented...
if isempty(SCS_TECHNOLOGY)
    SCS_TECHNOLOGY = 'st_bicmos6g';             % default value
end
if ~ischar(SCS_TECHNOLOGY)
    fprintf('Error: technology is not a string\n');
    disp(SCS_TECHNOLOGY)
    return
end
SCS_TECHNOLOGY = lower(SCS_TECHNOLOGY);
if isempty(strmatch(SCS_TECHNOLOGY, SCS_TECHNOLOGY_LIST, 'exact'))
    fprintf('Warning: unknown technology "%s"\n', SCS_TECHNOLOGY)
    fprintf('   using "dummy" technology...\n')
    SCS_TECHNOLOGY = 'dummy';
end
% . . valid device parameters per technology
validParameters = [];
% . . . . only use lower case names!!
switch SCS_TECHNOLOGY
    case 'st_bicmos6g'
        mos = {'w','l','nfing','number','screfirst','mismatch'};
        validParameters = setfield(validParameters, 'nmos', mos);
        validParameters = setfield(validParameters, 'pmos', mos);
        validParameters = setfield(validParameters, 'enmm9ju', mos);
        validParameters = setfield(validParameters, 'epmm9ju', mos);
        res = {'r','w','nc','number','mismatch'};
        validParameters = setfield(validParameters, 'res', res);
        validParameters = setfield(validParameters, 'rpo1', res);
        validParameters = setfield(validParameters, 'rpo1m1', res);
        validParameters = setfield(validParameters, 'rpo1m2', res);
        cap = {'c','w','nviatop','nviabot','number','mismatch'};
        validParameters = setfield(validParameters, 'cmim', cap);
        validParameters = setfield(validParameters, 'cm4m4b', cap);
        validParameters = setfield(validParameters, 'cm4m4bm1', cap);
        validParameters = setfield(validParameters, 'cm4m4bm2', cap);
    case 'dummy'
        % nothing
    otherwise
        % error
        disp('Error: in setting valid device parameters')
        return
end
% . . technology files
switch SCS_TECHNOLOGY
    case 'st_bicmos6g'
        modelFiles = {'rpolyh','rpolyn','rpolyp','rpolys','rndiff','rnwell', ...
                'cm1m2','cmim','cmim4b','cnwpo1','capane','capape','cspo1', ...
                'nmos','nmos_bs3','nrfmos', ...
                'pmos','pmos_bs3','prfmos', ...
                'npnvlv','npnvhv', ...
                'diode','varactor'};
        modelFiles = modelFiles';
        for i = 1:length(modelFiles)
            modelFiles{i} = [modelFiles{i} '.scs'];
        end
    case 'dummy'
        modelFiles = [];
    otherwise
        % error
        disp('Error: in setting model files')
        return
end


if debug
    fprintf('Technology = %s\n', SCS_TECHNOLOGY)
    struclist(validParameters, 'indent  ')
    modelFiles
    disp('Start of function scsSetPar')
    scsName
end
if ~isstruct(PAR)
    disp(PAR)
    fprintf('Warning : no parameters available -> no substitutions made\n')
    return
end


% open base file
if strcmp(scsName(length(scsName)-3:length(scsName)), '.scs')
    scsName = scsName(1:length(scsName)-4);
end
fid = sprintf('%s.scs', scsName);
fid = fopen(fid,'r');
if fid < 0
    fprintf('Error opening file %s.scs', scsName)
    status = -1;
    return
end
scsString = char(fread(fid)');
fclose(fid);

if processIncludeFiles
% find include statements
f1 = scsFindStatement_1(scsString,'include');
% substitute parameters also in these included files
for i = 1:length(f1)
    % find include file name
    % . . search begin -- ignore leading blanks, tabs or nextlines
    pos = f1(i) + 8;
    go = 1;
    while go
        if strncmp(scsString(pos:length(scsString)), blank, 1)
            pos=pos+1;
        elseif strncmp(scsString(pos:length(scsString)), tab, 1)
            pos=pos+1;
        elseif strncmp(scsString(pos:length(scsString)), nextline1, 2)
            pos=pos+2;
        elseif strncmp(scsString(pos:length(scsString)), nextline2, 2)
            pos=pos+2;
        else
            go = 0;
        end
    end
    % . . search end -- find first blank, tab or newline
    f2 = min([findstr(scsString(pos:length(scsString)), blank)]);
    f2 = min([findstr(scsString(pos:length(scsString)), tab) f2]);
    f2 = min([findstr(scsString(pos:length(scsString)), newline) f2]);
    if isempty(f2) f2 = length(scsString)+1; end
    % . . found
    incName = scsString(pos:pos-2+f2);
    % . . strip quotes
    if strcmp(incName(1),'"') & strcmp(incName(length(incName)),'"')
        incName = incName(2:length(incName)-1);
    end
    % . . interprete ~
    if strcmp(incName(1),'~')
        tmp = pwd;
        incName = [tmp(1:2) incName(2:length(incName))];
    end
    % . . skip model files
    model = 0;
    for j=1:length(modelFiles)
        if strncmp(fliplr(incName),fliplr(modelFiles{j}),length(modelFiles{j}))
            model = 1;
        end
    end
    % open and process this included file
    if ~model
        % this doesn't work if there is only a subckt definition in this file...
        [OrigPARINC, incstatus] = scsSetPar(incName,PAR);
        if incstatus < 0
            return
        else
            status = status + incstatus;
            % adjust OrigPAR structure
            if isstruct(OrigPARINC)
                parIncNames = fieldnames(OrigPARINC);
                for j = 1:length(parIncNames)
                    parIncName = parIncNames{j};
                    OrigPAR = setfield(OrigPAR, parIncName, getfield(OrigPARINC, parIncName));
                end
            end
        end
    end
end
end

% replace parameters in string
[newString, OrigPARSTR, strstatus] = scsSetParString(scsString,PAR, validParameters);
if strstatus < 0
    return
else
    status = status + strstatus;
    % adjust OrigPAR structure
    if isstruct(OrigPARSTR)
        parOrigNames = fieldnames(OrigPARSTR);
        for j = 1:length(parOrigNames)
            parOrigName = parOrigNames{j};
            OrigPAR = setfield(OrigPAR, parOrigName, getfield(OrigPARSTR, parOrigName));
        end
    end
end


% check number of replaced parameters
% . . only in last function call
[ST,index] = dbstack;
check = 0;
if length(ST) == 1
    check = 1;
else
    ST2 = ST(2);
    if isempty(strfind(ST2.name,'scsSetPar.m'))
        check = 1;
    end
end
% . . do the check
if check
    fout = 0;
    parString = struclist(PAR, 'sort');
    orgString = struclist(OrigPAR, 'sort');
    parNewLines = strfind(parString, newline);
    orgNewLines = strfind(orgString, newline);
    parStart = min(strfind(parString, '.'));
    orgStart = min(strfind(orgString, '.'));
    if ~isstruct(OrigPAR)
        fout = 1;
    elseif length(parNewLines) ~= length(orgNewLines)
        fout = 1;
    else
        parNewLines = [0 parNewLines];
        orgNewLines = [0 orgNewLines];
        for i = 2:length(orgNewLines)
            parLine = parString(parNewLines(i-1)+parStart:parNewLines(i));
            orgLine = orgString(orgNewLines(i-1)+orgStart:orgNewLines(i));
            poseq = max(strfind(parLine, '='),strfind(parLine, '='));
            if ~strncmp(parLine, orgLine, poseq)
                fout = 1;
            end
        end
    end
    if fout
        fprintf('Warning : not all given parameters were properly substituted\n');
        status = -1;
        INPUT = PAR;
        struclist(INPUT, 'prefix      ', 'indent  ', 'sort');
        ORIGINAL = OrigPAR;
        struclist(ORIGINAL, 'prefix   ', 'indent  ', 'sort');
    end
end


% dump new string into file
if strcmp(scsString,newString)
    if debug
        fprintf('no subs in %s\n', scsName)
    end
else
    % make back-up copy
    fid = sprintf('.scsSetPar.%s.scs',scsName);
    fid = fopen(fid,'w');
    fprintf(fid,'%s',scsString);
    fclose(fid);
    % write new file
    fid = sprintf('%s.scs',scsName);
    fid = fopen(fid,'w');
    fprintf(fid,'%s',newString);
    fclose(fid);
end



% =========================================================================
% sub-function definition
% =========================================================================

function  [newString, OrigPAR, status] = scsSetParString(scsString,PAR, validParameters)

% updates the parameters in a character string

debug = 0;
if debug
    disp('Start of function scsSetParString')
    scsString
end

[ST,index] = dbstack;
if debug
    for i = 1:length(ST)
        STi = ST(i);
        fprintf('ST %i = %s\n', i, STi.name);
    end
end

% initialize
OrigPAR = [];
status = 0;
parNames = fieldnames(PAR);
newline = sprintf('\n');
nextline1 = sprintf('\n+');
nextline2 = sprintf('\\\n');
tab = sprintf('\t');
blank = ' ';


% find subcircuit definitions
% . . subckt statement
f1 = scsFindStatement_1(scsString, 'subckt');
% . . ends statement
[tmp, f2] = scsFindStatement_1(scsString, 'ends');
% . . check
if length(f1) ~= length(f2)
    fprintf('Error: incorrect subcircuit definitions\n')
    fprintf('  %i subckt statements <-> %i ends statements\n', length(f1), length(f2))
    f1
    f2
    newString = scsString;
    status = -1;
    return
end
% . . (f1 f2) pair should be removed if they are for a subckt definition that is
% . . the complete string
% . . but that is not true is this function is called from the function scsSetPar,
% . . because then the complete file is a subckt!
if length(ST) > 1
    ST2 = ST(2);
    if strfind(ST2.name,'scsSetPar.m (scsSetParString)') & length(f1) > 0
        if (f1(1)==1 & f2(length(f2))==length(scsString))
            f1 = f1(2:length(f1));
            f2 = f2(1:length(f2)-1);
        end
    end
end
f12 = sort([f1 f2]);
f1new = [];
f2new = [];
i = 1;
j = 1;
while i <= length(f1)
    % . . take a subckt statement
    f1new(j) = f1(i);
    % . . find position of this subckt statement in f12 vector
    p12 = find(f12-f1(i)==0);
    % . . scan the rest of f12 to find the ends statement of the correct level
    level = 1;
    while (level ~= 0) & (p12 < length(f12))
        p12 = p12+1;
        if find(f2-f12(p12)==0)
            level = level - 1;
        else
            level = level + 1;
            % skip next subckt statement
            i = i + 1;
        end
    end
    % . . add this position to list of valid ends statements
    f2new(j) = f12(p12);
    i = i + 1;
    j = j + 1;
end
f1 = f1new;
f2 = f2new;


% Treat subckt parameters definitions
for i = 1:length(f1)
    subString = scsString(f1(i):f2(i));
    % error check
    if ~strncmp(subString,'subckt',6)
        fprintf('Error: no subckt definition:\n')
        fprintf('%s', subString)
        fprintf('### end-of-error ###\n')
        return
    end
    % retrieve subckt name
    % . . begin -- ignore leading blanks, tabs or nextlines
    pos = 7;
    go = 1;
    while go
        if strncmp(subString(pos:length(subString)), blank, 1)
            pos=pos+1;
        elseif strncmp(subString(pos:length(subString)), tab, 1)
            pos=pos+1;
        elseif strncmp(subString(pos:length(subString)), nextline1, 2)
            pos=pos+2;
        elseif strncmp(subString(pos:length(subString)), nextline2, 2)
            pos=pos+2;
        else
            go = 0;
        end
    end
    % . . end -- find first blank, tab or nextline
    f3 = min([findstr(subString(pos:length(subString)), blank)]);
    f3 = min([findstr(subString(pos:length(subString)), tab) f3]);
    f3 = min([findstr(subString(pos:length(subString)), nextline1) f3]);
    f3 = min([findstr(subString(pos:length(subString)), nextline2) f3]);
    % . . found
    subName = subString(pos:pos-2+f3);
    if debug
        fprintf('subckt %i of %i : %s\n', i, length(f1), subName)
        subString
    end
    % find associated parameter structure
    % and do the replace operation
    if isfield(PAR, subName)
        PARSUB = getfield(PAR, subName);
        if debug
            struclist(PARSUB);
        end
        % replace parameters
        [newSubString, OrigPARSUB, substatus] = scsSetParString(subString, PARSUB, validParameters);
        if debug
            substatus
            struclist(OrigPARSUB);
        end
        if substatus >= 0
            status = status + 1;
        else
            status = -1
            return
        end
        % update structure with origian parameters
        if isstruct(OrigPARSUB)
            OrigPAR = setfield(OrigPAR, subName, OrigPARSUB);
        end
        % mark this part of the structure as done
        k = strmatch(subName,parNames);
        if k > 0
            parNames{k} = '.';
        end
        % plug new subString into scsString and update positions
        scsString = [scsString(1:f1(i)-1) newSubString scsString(f2(i)+1:length(scsString))];
        f1(find(f1>=f2(i))) = f1(find(f1>=f2(i))) + length(newSubString) - length(subString);
        f2(find(f2>=f2(i))) = f2(find(f2>=f2(i))) + length(newSubString) - length(subString);
    end
end


% Check for special parameters
if strmatch('temp', parNames)
    k = strmatch('temp', parNames);
    % Find options statements
    [f3, f4] = scsFindStatement_2(scsString, 'options');
    % Browse in reverse direction through options statements
    for i = length(f3):-1:1
        statement = scsString(f3(i):f4(i));
        if debug
            statement
        end
        % Find 'par = value' constructions
        [pos_par_beg, pos_par_end, pos_val_beg, pos_val_end] = scsFindParEqVal(statement);
        % Match with parameter names
        for j = length(pos_par_beg):-1:1
            % name
            parName = char(parNames(k));
            if strcmp(statement(pos_par_beg(j):pos_par_end(j)), parName)
                % matching parameter name
                if debug
                    par_eq_val = statement(pos_par_beg(j):pos_val_end(j))
                end
                % . . mark as matched
                parNames{k} = '.';
                % . . store present value in original string
                OrigPAR = setfield(OrigPAR, parName, statement(pos_val_beg(j):pos_val_end(j)));
                % . . get new value
                newParValue = getfield(PAR, parName);
                if isnumeric(newParValue)
                    newParValue = eng(newParValue);
                    if strcmp(newParValue(end), blank)
                        newParValue = newParValue(1:end-1);
                    end
                end
                if ischar(newParValue)
                    if strcmp(lower(newParValue), 'remove')
                        % . . remove parameter
                        stat1 = statement(1:pos_par_beg(j)-1);
                        if pos_val_end(j) < length(statement)
                            stat2 = statement(pos_val_end(j)+1:length(statement));
                        else
                            stat2 = '';
                        end
                        statement = [stat1 stat2];
                        status = status + 1;
                    elseif strcmp(lower(newParValue), 'no change')
                        % . . don't change
                    else
                        % . . replace value
                        stat1 = statement(1:pos_val_beg(j)-1);
                        stat2 = statement(pos_val_end(j)+1:length(statement));
                        statement = [stat1 newParValue stat2];
                        status = status + 1;
                    end
                else
                    disp('Error: value not useable')
                    status=-1;
                    break
                end
            end
        end
        % Replace complete parameter statement
        scsString = [scsString(1:f3(i)-1) statement scsString(f4(i)+1:length(scsString))];
    end
end


% Search for parameter statement
[f3, f4] = scsFindStatement_1(scsString, 'parameters');
% . . skip all occurences that occur within subckt's
for i = 1:length(f1)
    ok = find(~(f1(i) < f3 & f3 < f2(i)));
    f3 = f3(ok);
    f4 = f4(ok);
end

% browse in reverse direction through the parameter definitions
for i = length(f3):-1:1
    statement = scsString(f3(i):f4(i));
    if debug
        statement
    end
    % Find 'par = value' constructions
    [pos_par_beg, pos_par_end, pos_val_beg, pos_val_end] = scsFindParEqVal(statement);
    % Match with parameter names
    for j = length(pos_par_beg):-1:1
        for k = 1:length(parNames)
            % name
            parName = char(parNames(k));
            if strcmp(statement(pos_par_beg(j):pos_par_end(j)), parName)
                % matching parameter name
                if debug
                    par_eq_val = statement(pos_par_beg(j):pos_val_end(j))
                end
                % . . mark as matched
                parNames{k} = '.';
                % . . store present value in original string
                OrigPAR = setfield(OrigPAR, parName, statement(pos_val_beg(j):pos_val_end(j)));
                % . . get new value
                newParValue = getfield(PAR, parName);
                if isnumeric(newParValue)
                    newParValue = eng(newParValue);
                    if strcmp(newParValue(end), blank)
                        newParValue = newParValue(1:end-1);
                    end
                end
                if ischar(newParValue)
                    if strcmp(lower(newParValue), 'remove')
                        % . . remove parameter
                        stat1 = statement(1:pos_par_beg(j)-1);
                        if pos_val_end(j) < length(statement)
                            stat2 = statement(pos_val_end(j)+1:length(statement));
                        else
                            stat2 = '';
                        end
                        statement = [stat1 stat2];
                        status = status + 1;
                    elseif strcmp(lower(newParValue), 'no change')
                        % . . don't change
                    else
                        % . . replace value
                        stat1 = statement(1:pos_val_beg(j)-1);
                        stat2 = statement(pos_val_end(j)+1:length(statement));
                        statement = [stat1 newParValue stat2];
                        status = status + 1;
                    end
                else
                    disp('Error: value not useable')
                    status=-1;
                    break
                end
            end
        end
    end
    % Replace complete parameter statement
    scsString = [scsString(1:f3(i)-1) statement scsString(f4(i)+1:length(scsString))];
end


% Now start replacing device parameters
if debug
    fprintf('Remaining parameters for device substitutions:');
    parNames
end
for i=1:length(parNames)
    devName = char(parNames(i));
    if ~strcmp(devName,'.')
        % search device statement
        [f3, f4] = scsFindStatement_1(scsString, devName);
        if length(f3) > 1
            fprintf('Error: multiple definitions for device "%s" found ', devName);
            fprintf('in string:\n%s\n', scsString);
            status = -1
            return
        end
        if ~isempty(f3)
            statement = scsString(f3:f4);
            DEVPAR = getfield(PAR,devName);
            OrigDEVPAR = [];
            if ~isstruct(DEVPAR)
                % TBD
            else
                devParNames = fieldnames(DEVPAR);
                % find 'par = value' constructions
                [pos_par_beg, pos_par_end, pos_val_beg, pos_val_end, devType] = ...
                    scsFindParEqVal(statement);
                if debug
                    statement
                    DEVPAR
                    devType
                end
                % find and replace matching device parameters
                for j = length(pos_par_beg):-1:1
                    origDevParName = statement(pos_par_beg(j):pos_par_end(j));
                    if isfield(DEVPAR, origDevParName)
                        % mark field as done
                        devParNames{strmatch(origDevParName, devParNames, 'exact')} = '.';
                        % store present value in original string
                        OrigDEVPAR = setfield(OrigDEVPAR, origDevParName, statement(pos_val_beg(j):pos_val_end(j)));
                        % replace parameter
                        newDevParValue = getfield(DEVPAR, origDevParName);
                        if isnumeric(newDevParValue)
                            newDevParValue = eng(newDevParValue);
                            if strcmp(newDevParValue(end), blank)
                                newDevParValue = newDevParValue(1:end-1);
                            end
                        end
                        if ischar(newDevParValue)
                            switch lower(newDevParValue)
                                case {'new parameter','remove','invalid'}
                                    % . . remove parameter
                                    stat1 = statement(1:pos_par_beg(j)-1);
                                    if pos_val_end(j) < length(statement)
                                        stat2 = statement(pos_val_end(j)+1:length(statement));
                                    else
                                        stat2 = '';
                                    end
                                    statement = [stat1 stat2];
                                    status = status + 1;
                                case 'no change'
                                    % . . don't change
                                otherwise
                                    % . . replace value
                                    stat1 = statement(1:pos_val_beg(j)-1);
                                    if pos_val_end(j) < length(statement)
                                        stat2 = statement(pos_val_end(j)+1:length(statement));
                                    else
                                        stat2 = '';
                                    end
                                    statement = [stat1 newDevParValue stat2];
                                    status = status + 1;
                            end
                        else
                            disp('Error: value not useable')
                            status=-1;
                            break
                        end
                    end
                end
            end
            % find new positions of 'par = value' constructions
            [pos_par_beg, pos_par_end, pos_val_beg, pos_val_end, devType] = ...
                scsFindParEqVal(statement);
            % add parameters that are not changed yet
            % only for the list of valid parameters for this device
            for j = 1:length(devParNames)
                devParName = devParNames{j};
                if ~strcmp(devParName, '.')
                    % mark field as done
                    devParNames{j} = '.';
                    % new parameter value to be added
                    newDevParValue = getfield(DEVPAR, devParName);
                    % check if valid parameter
                    devType = lower(devType);
                    if isfield(validParameters, devType)
                        % . . there is a list of valid parameters
                        if isempty(strmatch(lower(devParName), getfield(validParameters, devType), 'exact'))
                            newDevParValue = 'invalid';
                        end
                    else
                        % . . there is no list of valid parameters
                    end
                    if isnumeric(newDevParValue)
                        newDevParValue = eng(newDevParValue);
                        if strcmp(newDevParValue(end), blank)
                            newDevParValue = newDevParValue(1:end-1);
                        end
                    end
                    if ischar(newDevParValue)
                        switch lower(newDevParValue)
                            case 'new parameter'
                                % . . probably something went wrong...
                                OrigDEVPAR = setfield(OrigDEVPAR, devParName, 'no change');
                            case 'remove'
                                % . . probably something went wrong...
                                OrigDEVPAR = setfield(OrigDEVPAR, devParName, 'no change');
                            case 'no change'
                                % . . probably something went wrong...
                                OrigDEVPAR = setfield(OrigDEVPAR, devParName, 'no change');
                            case 'invalid'
                                % . . invalid parameter, don't add
                                OrigDEVPAR = setfield(OrigDEVPAR, devParName, 'invalid');
                            otherwise
                                % . . add value at end
                                stat1 = statement(1:max(pos_val_end));
                                if max(pos_val_end)<length(statement)
                                    stat2 = statement(max(pos_val_end)+1:length(statement));
                                else
                                    stat2 = '';
                                end
                                stat3 = ['  ' devParName '=' newDevParValue ' '];
                                statement = [stat1 stat3 stat2];
                                status = status + 1;
                                % . . add 'new parameter' to original structure
                                OrigDEVPAR = setfield(OrigDEVPAR, devParName, 'new parameter');
                        end
                    else
                        disp('Error: value not useable')
                        status=-1;
                        break
                    end
                end
            end
            % replace statement
            scsString = [scsString(1:f3-1) statement scsString(f4+1:length(scsString))];
            % update original structure
            if isstruct(OrigDEVPAR)
                OrigPAR = setfield(OrigPAR, devName, OrigDEVPAR);
            end
        end         
    end
end

newString = scsString;



% =========================================================================

function [begin, einde] = scsFindStatement_1(scsString, statement)

% [begin, einde] = scsFindStatement_1(scsString,statement)
%
% searches through the spectre-syntax string "scsString" looking for the
% spectre statement "statement" and returns the positions of correct
% matches
% "statement" must be the FIRST word of a line
%
%   begin = vector with begin position
%   einde = vector with end position (newline character)


% initialise
debug = 0;
newline = sprintf('\n');
nextline1 = sprintf('\n+');
nextline2 = sprintf('\\\n');
tab = sprintf('\t');
blank = ' ';

if debug
    fprintf('search statement "%s"\n', statement)
end

% replace tabs and nextlines
scsString = strrep(scsString, tab, char(blank*ones(size(tab))));
scsString = strrep(scsString, nextline1, char(blank*ones(size(nextline1))));
scsString = strrep(scsString, nextline2, char(blank*ones(size(nextline2))));


% first search
% . . statement followed by space
s1 = [statement blank];
f1 = findstr(scsString,s1);
% . . statement followed by newline
s1 = [statement newline];
f1 = sort([f1 findstr(scsString,s1)]);
if debug > 1
    f1
end

% validate search results: must be preceded by a newline, or be the first
%     character of the string
for i = 1:length(f1)
    % skip leading blanks
    while f1(i)>1 & strcmp(scsString(f1(i)-1), blank)
        f1(i) = f1(i)-1;
    end
    % search newline
    if f1(i) == 1
        % . . first char, ok
    else
        % . . check for newline
        if ~strcmp(scsString(f1(i)-1), newline)
            f1(i) = 0;
        end
    end
end

% resulting correct matches
begin = f1(find(f1~=0));


% Search end-of-line statements
% . . search end-of-line
f2 = findstr(scsString, newline);
% . . add end-of-string
f2 = [f2 length(scsString)];

% construct einde-vector
einde = [];
for i=1:length(begin)
    einde(i) = f2(min(find(f2 - begin(i) > 0)));
    if debug
        fprintf('statement found:\n%s\n',(scsString(begin(i):einde(i))))
    end
end



% =========================================================================

function [begin, einde] = scsFindStatement_2(scsString, statement)

% [begin, einde] = scsFindStatement_2(scsString,statement)
%
% searches through the spectre-syntax string "scsString" looking for the
% spectre statement "statement" and returns the positions of correct
% matches
% "statement" must be the SECOND word of a line
%
%   begin = vector with begin position
%   einde = vector with end position (newline character)


% initialise
debug = 0;
newline = sprintf('\n');
nextline1 = sprintf('\n+');
nextline2 = sprintf('\\\n');
tab = sprintf('\t');
blank = ' ';

if debug
    fprintf('search statement "%s"\n', statement)
end

% replace tabs and nextlines
scsString = strrep(scsString, tab, char(blank*ones(size(tab))));
scsString = strrep(scsString, nextline1, char(blank*ones(size(nextline1))));
scsString = strrep(scsString, nextline2, char(blank*ones(size(nextline2))));


% first search
% . . statement preceded and followed by blank
s1 = [blank statement blank];
f1 = findstr(scsString,s1);
% . . statement preceded by blank and followed by newline
s1 = [blank statement newline];
f1 = sort([f1 findstr(scsString,s1)]);
if debug > 1
    f1
end

% validate search results: the statement must be the second word of a new line.
for i = 1:length(f1)
    pos = f1(i);
    % . . skip leading blanks and tabs
    go = 1;
    while go
        if pos == 1
            prev = scsString(pos);
            go = 0;
        else
            pos = pos-1;
            prev = scsString(pos);
            go = strcmp(prev,blank) | strcmp(prev,tab);
        end
    end
    % . . find preceding newline
    f2 = max([1 findstr(scsString(1:pos), newline)+1]);
    % . . there can't be any blanks between f2 and pos
    if isempty(findstr(scsString(f2:pos), blank))
        % . . ok
        f1(i) = f2;
    else
        % . . not ok
        f1(i) = 0;
    end
end
if debug > 1
    f1
end

% resulting correct matches
begin = f1(find(f1~=0));


% Search end-of-line statements
% . . search end-of-line
f2 = findstr(scsString, newline);
% . . add end-of-string
f2 = [f2 length(scsString)];

% construct einde-vector
einde = [];
for i=1:length(begin)
    einde(i) = f2(min(find(f2 - begin(i) > 0)));
    if debug
        fprintf('statement found:\n%s\n',(scsString(begin(i):einde(i))))
    end
end



% =========================================================================

function [pos_par_beg, pos_par_end, pos_val_beg, pos_val_end, devType] = scsFindParEqVal(statement)

% Find 'par = value' constructions in the spectre statement "statement"
% and
% Return the name ofthe type of the device (model or subckt name)


% initialise
newline = sprintf('\n');
nextline1 = sprintf('\n+');
nextline2 = sprintf('\\\n');
tab = sprintf('\t');
blank = ' ';

pos_par_beg = [];
pos_par_end = [];
pos_val_beg = [];
pos_val_end = [];
devType = '';

% replace tabs and nextlines
statement = strrep(statement, tab, char(blank*ones(size(tab))));
statement = strrep(statement, nextline1, char(blank*ones(size(nextline1))));
statement = strrep(statement, nextline2, char(blank*ones(size(nextline2))));



% Find '=' signs
pos_eq = findstr(statement, '=');

% Browse through them in reverse direction
for j = length(pos_eq):-1:1
    % find parameter name, start from '='-sign
    pos = pos_eq(j) - 1;
    % . . end: skip blanks and tabs
    while strcmp(statement(pos), blank) | strcmp(statement(pos), tab)
        pos = pos - 1;
    end
    pos_par_end(j) = pos;
    % . . begin: find leading blank, tab, '+', or newline
    while isempty(strmatch(statement(pos),{blank, tab, '+', newline}))
        pos = pos - 1;
    end
    pos_par_beg(j) = pos + 1;
    % find value definition
    pos = pos_eq(j) + 1;
    % . . begin: skip blanks and tabs
    while strcmp(statement(pos), blank) | strcmp(statement(pos), tab)
        pos = pos + 1;
    end
    pos_val_beg(j) = pos;
    % . . end: start from begin of next parameter name
    if j==length(pos_eq)
        pos = length(statement);
    else
        pos = pos_par_beg(j+1)-1;
    end
    % . . end: skip blanks, tabs, '+', '\', '\n'
    while strmatch(statement(pos),{blank, tab, '+', '\', newline})
        pos = pos - 1;
    end
    pos_val_end(j) = pos;
end

% Find device type
if nargout > 4
    % now look for the first word before pos_par_beg(1)
    if ~isempty(pos_par_beg)
        tmp = statement(1:pos_par_beg(1)-1);
    else
        tmp = statement;
    end
    % remove trailing blanks
    tmp = deblank(tmp);
    % find blanks
    pos_blank = [0 findstr(tmp, blank)];
    devType = tmp(max(pos_blank)+1:end);
end