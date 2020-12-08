function varargout = struclist(struc, varargin)

% string = STRUCLIST(struc, [options])
%     displays all the contents of the structure 'struc', optionally sorted
%     alphabetically
%     Field values in the structure can be either strings, single numbers or
%     vectors
%     Options: 'sort' -> sorted alphabetically
%              'full' -> full field names instead of leading dots
%              'indent char' -> use 'char' as indent character (default: '.')
%              'sep char' -> use 'char' as hierarchy separator (default: '.')
%              'prefix str' -> use 'str' as prefix (default: none)
%              'eq str' -> use 'str' as equals sign (default: ' = ')
%
% JCx, 14-10-2003


persistent indent indentstr

% parse input arguments
sorteer = 0;
sortstr = 'none';
fulldot = 0;
fullstr = 'none';
indentchar = '.';
indentarg = 'none';
sepchar = '.';
separg = 'none';
prefixstr = '';
prefixarg = 'none';
equalstr = ' = ';
equalarg = 'none';
for i = 1:nargin-1
    arg = varargin{i};
    if ischar(arg)
        if strcmp(arg, 'sort')
            sorteer = 1;
            sortstr = 'sort';
        elseif strcmp(arg, 'full')
            fulldot = 1;
            fullstr = 'full';
        elseif strncmp(arg, 'indent ', 7)
            if length(arg) > 7
                indentchar = arg(8);
                indentarg = ['indent ' indentchar];
                fulldot = 0;
                fullstr = 'none';
            end
        elseif strncmp(arg, 'sep ', 4)
            if length(arg) > 4
                sepchar = arg(5:length(arg));
            else
                sepchar = '';
            end
            separg = ['sep ' sepchar];
        elseif strncmp(arg, 'prefix ', 7)
            if length(arg) > 7
                prefixstr = arg(8:length(arg));
                prefixarg = arg;
            end
        elseif strncmp(arg, 'eq ', 3)
            if length(arg) > 3
                equalstr = arg(4:length(arg));
                equalarg = arg;
            end
        elseif ~strcmp(arg,'none')
            fprintf('Error: incorrect argument "%s"\n', arg)
            return
        end
    else
        fprintf('Error: incorrect argument "%s"\n', inputname(i+1))
        return
    end
end

% initialize
string = '';
if isempty(indent)
    indent = 0;
    indentstr = '';
end
if length(inputname(1)>0)
    % first call, initialize indent
    backspace = sprintf('\b');
    if findstr(prefixstr, backspace)
        string = sprintf('\n%s', string);
    end
    string = [string prefixstr sprintf('%s', inputname(1))];
    indent = length(inputname(1));
    indentstr = [prefixstr inputname(1)];
end

% process
if isempty(struc)
    string = [string sprintf(' is empty\n')];
elseif isstruct(struc)
    list = fieldnames(struc);
    if sorteer
        list = sort(list);
    end
    for i=1:length(list)
        list(i);
        if i ~=1
            if ~fulldot
                string = [string prefixstr char(ones(1,indent)*double(indentchar))];
            else
                string = [string indentstr];
            end
        end
        string = [string sepchar sprintf('%s',list{i})];
        if isstruct(getfield(struc,list{i}))
            indent = indent+length(sepchar)+length(num2str(list{i}));
            indentstr = [indentstr sepchar list{i}];
            tmp = struclist(getfield(struc,list{i}), sortstr, fullstr, indentarg, separg, prefixarg, equalarg);
            string = [string tmp];
            indent = indent-length(sepchar)-length(num2str(list{i}));
            indentstr = indentstr(1:length(prefixstr)+indent);
        else
            string = [string equalstr];
            val = getfield(struc,list{i});
            if isempty(val)
                string = [string '[]'];
            elseif ischar(val)
                string = [string sprintf('%s', val)];
            elseif isnumeric(val)
                V = whos('val');
                [nr, nc] = size(val);
                if nr <= 1 & nc <= 1
                    string = [string sprintf('%s', eng(val))];
                elseif (nr == 1 & nc <= 8) | (nc == 1 & nr <= 8)
                    string = [string sprintf('[%s', eng(val(1)))];
                    for j = 2:length(val)
                        string = [string sprintf(' %s', eng(val(j)))];
                    end
                    string = [string ']'];
                    if nr > nc
                        string = [string ''''];
                    end
                else
                    string = [string sprintf('[%ix%i %s]', nc, nr, V.class)];
                end
            elseif iscell(val)
                V = whos('val');
                [nr, nc] = size(val);
                if nr <= 1 & nc <= 1
                    string = [string sprintf('%s', val{1})];
                elseif (nr == 1 & nc <= 8) | (nc == 1 & nr <= 8)
                    string = [string sprintf('{%s', val{1})];
                    for j = 2:length(val)
                        string = [string sprintf(' %s', val{j})];
                    end
                    string = [string '}'];
                    if nr > nc
                        string = [string ''''];
                    end
                else
                    string = [string sprintf('{%ix%i %s}', nc, nr, V.class)];
                end
            else
                string = [string sprintf(' error')];
            end
            string = [string sprintf('\n')];
        end
    end
else
    string = [string sprintf(' is no structure\n')];
end

if nargout == 0
    fprintf(string);
else
    varargout{1} = string;
end
