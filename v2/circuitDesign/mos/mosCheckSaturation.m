function result = mosCheckSaturation(MOS)
%MOSCHECKSATURATION check whether a transistor operates in the saturation
%    region.
%
%    mosCheckSaturation(MOS) gives a warning when transistor MOS is not
%    saturated. This function compares MOS.vds to MOS.vdsat if they exist.
%    When one of these fields does not exist, an error is signaled.
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

result = 1;
if not(isfield(MOS, 'vds')) error('vds is not defined for transistor %s', ...
        inputname(1)); end
if not(isfield(MOS, 'vdsat')) error('vdsat is not defined for transistor %s', ...
        inputname(1)); end
if isNtype(mosTable(MOS))
    if (MOS.vds < MOS.vdsat)
        warning('%s is not saturated: vds = %g < vdsat = %g\n', MOS.name, ...
            MOS.vds, MOS.vdsat);
        result = 0;
    end
else % we have a p-MOS
    if (MOS.vds > MOS.vdsat)
        warning('%s is not saturated: vds = %g > vdsat = %g\n', MOS.name, ...
            MOS.vds, MOS.vdsat);
        result = 0;
    end
end
