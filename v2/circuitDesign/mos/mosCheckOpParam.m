function result = mosCheckOpParam(string)
%MOSCHECKOPPARAM checks whether a given STRING is an operating point
%   parameter that is supported in the table-based circuit sizing approach.
%
%   MOSCHECKOPPARAM(STRING) returns 1 if the given STRING is the name of an 
%   operating point parameter that is supported in the table-base circuit sizing
%   approach. If the STRING does not correspond to the name of a supported
%   parameter, then 0 is returned.
%   A given table may contain only a subset of the supported operating
%   point parameters.
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

switch string
    case {'ids', 'gm', 'gmbs', 'gds', 'vth', 'vdsat', 'cgs', 'cgb', ...
                'cgg', 'cdd', 'cgd', 'csbI', 'cdbI', 'di2_id', 'di2_fn', ...
                'cgsol', 'cgdol', 'igb', 'igd', 'igs'}
        result = 1;
    otherwise
        result = 0;
end
