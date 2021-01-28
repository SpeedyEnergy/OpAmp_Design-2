function tab = mosTable(MOS)
%MOSTABLE returns the table of operating point parameters from which the
%intrinsic operating point of a MOS transistor are derived.
%
%    TABLE = mosTable(MOS) returns the table (this is a structure) from
%    which the intrinsic operating point parameters of the given MOS
%    transistor are derived. 
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


tab = [];

% first try the name given in the MOS structure
if isfield(MOS,'table')
    tab = evalin('base', MOS.table, '[]');
end

% then try the gobal N,P
if isempty(tab)
    global N P;
    switch MOS.eltype
        case 'nmos'
            tab = N;
        case 'pmos'
            tab = P;
        otherwise
            fprintf('Error: unknown table\n');
            tab = [];
    end
end
