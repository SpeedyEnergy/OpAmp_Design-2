function techCapDigital(table, capParam, lg, vsb, figNumber)
%TECHCAPS plots intrinsic parasitic (trans)capacitance for |vgs| = VDD.
%
%    techCapDigital(TABLE, CAPPARAM, LG, VSB, FIGNUMBER) plots the intrinsic
%    parasitic capacitance CAPPARAM as a function of VDS.
%    The plot is generated for a MOS transistor characterized by the
%    given TABLE of intrinsic operating point parameters.   
%    The length of the channel and vsb need to be specified as an argument
%    (LG and VSB, respectively). These values should belong to the arrays 
%    tableInArray('lg', TABLE) and tableInArray('vsb', TABLE), respectively.
%
%    See also tableInArray.
%
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


if not(ismember(capParam, fieldnames(table.Table)))
  error('The specified parameter is not a valid capacitance parameter');
end

vds = abs(tableInArray('vds', table));
lgIndex = find(tableInArray('lg', table) == lg);
if isempty(lgIndex)
  error('The specified value of lg, namely %g does not belong to the array of Lg values', ...
  lg);
end

vsbIndex = find(tableInArray('vsb', table) == vsb);
if isempty(vsbIndex)
  error('The specified value of vsb, namely %g does not belong to the array of vsb values', ...
  vsb);
end


eox = 3.45e-11;
tox = tableModelParam('tox', table);
cox = eox/tox;
coxT = cox*tableWref(table)*lg;


cap = squeeze(table.Table.(capParam)(lgIndex, tableInLength('vgs', ...
      table), : , vsbIndex))/coxT; 

axisCap = [0 1.2*vds(length(vds)) -0.1 1.1]; 


% plot of intrinsic capacitance parameter:
figure(figNumber);
clf;
hold off;
plot(vds, abs(cap), 'linewidth', 2);
%axis(axisCap);
title(sprintf(...
    'normalized %s (vsb = %g V, Lg = %g um, %s-MOS)', ...
    capParam, vsb, 1e6*lg, tableType(table)));
xlabel('V_{DS} (V)');
ylabel(sprintf('%s (normalized)', capParam));
grid on;
set(gcf, 'name', capParam);

hold off
