function cjTot = mosJuncap(table, voltage, area, perimTot, w)
%MOSJUNCAP computation of junction capacitance of a source or drain region 
%of a MOS transistor.
%
%    cJunction = mosJuncap(TABLE, VDIFFBULK, AREA, PERIMTOT, W) returns
%    the value of the junction capacitance of a source or drain region for a
%    given type of a MOS transistor, specified by the given TABLE. For the
%    junction one must specify the AREA, the total perimeter PERIMTOT of the
%    drain or source region, the gate width W and the voltage over the
%    junction VDIFFBULK. The latter is the voltage difference between the
%    source or drain region and the bulk (should be positive for an n-MOS and 
%    negative for a p-MOS) 
%
%    EXAMPLE :
%
%       for a MOS transistor with one finger, a width of 2 micrometer, a
%       length of 1 micrometer for the diffusion areas, and 1 Volt over the
%       junction between the diffusion area of the bulk and a table name N, the
%       function call is
%       cdb = mosJuncap(N, 1, 2e-12, 6e-6, 2e-6);
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%



cj = tableModelParam('cj', table);
cjsw = tableModelParam('cjsw', table);
if isModelParam('cjswg', table)
  cjswg = tableModelParam('cjswg', table);
else
  cjswg = cjsw;
end
mj = tableModelParam('mj', table);
mjsw = tableModelParam('mjsw', table);
if isModelParam('mjswg', table)
  mjswg = tableModelParam('mjswg', table);
else
  mjswg = mjsw;
end
pb = tableModelParam('pb', table);
pbsw = tableModelParam('pbsw', table);
if isModelParam('pbswg', table)
  pbswg = tableModelParam('pbswg', table);
else
  pbswg = pbsw;
end


if not(isNtype(table))
  % we have a p-MOS
  voltage = voltage * (-1);
end 

if voltage < 0 % junction is forward biased
  if abs(voltage) > tableModelParam('pb', table)/2
    error(sprintf(...
	'junction to bulk of %s-MOS is forward biased : |voltage| = %g',...
	tableType(table), voltage));
  else
    warning(...
	sprintf(...
	'junction to bulk of %s-MOS slightly forward biased : |voltage| = %g',...
	tableType(table), voltage));
  end
end


cArea = cj * (1 + (voltage/pb))^(-mj); % capacitance per unit area
csw = cjsw * (1 + (voltage/pbsw))^(-mjsw); % cap. per unit perimeter
cswg = cjswg * (1 + (voltage/pbswg))^(-mjswg); % cap. per unit length
  


cjTot = area * cArea + (perimTot - w) * csw + w * cswg;
