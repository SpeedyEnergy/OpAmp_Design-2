function tableDisplay(table);
%TABLEDISPLAY view the parameters that are present in a table for a MOS transistor
%
%   tableDisplay(TABLE) displays information contained in the table that has
%   been constructed for a given transistor type.
%   The following information is displayed:
%     1) DIMENSIONS: minimum gate length and width, width for which the table
%        has been constructed (so-called reference width), width below which
%        narrow-channel effects become visible (so-called critical width).
%     2) A list of the input variables and their ranges
%     3) A list of operating point parameters that are stored in the table
%         (gm, ids, cgs, vdsat), .... Some parameters scale proportionally to
%         the transistor width (at least above the critical width), other ones 
%         (e.g. vdsat, vth) do not.
%
%   EXAMPLE :
%
%     tableDisplay(N) 
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%


lineLength = 80;
interString = ' -- ';
blanks = 1;

str1 = '';

GeoScale = 1e-6;
scalingString = scaleString(GeoScale);

% minimum L:
str1 = sprintf('%sMinimum L = %sm\n', str1, eng(tableLmin(table))); 
% minimum W:
str1 = sprintf('%sMinimum W = %sm\n', str1, eng(tableWmin(table))); 
% reference W:
str1 = sprintf('%sReference W = %sm\n', str1, eng(tableWref(table)));
% critical W:
str1 = sprintf('%sCritical W = %sm\n', str1, eng(tableWcrit(table)));

% blank lines:
str1 = sprintf('%s\n', str1);

% input variables:
str1 = sprintf('%sInput variables:  ', str1);

inputNames = tableInNames(table);
for i = 1:length(inputNames)
  str1 = sprintf('%s%s  ', str1, char(inputNames(i))); 
end

% blank lines:
str1 = sprintf('%s\n\n', str1);

str1 = sprintf('%s vgs from %g to %g V in steps of %g V\n', str1, ...
    tableInInit('vgs', table), tableInFinal('vgs', table), ...
    tableInStep('vgs', table));

str1 = sprintf('%s vds from %g to %g V in steps of %g V\n', str1, ...
    tableInInit('vds', table), tableInFinal('vds', table), ...
    tableInStep('vds', table));

str1 = sprintf('%s vsb from %g to %g V in steps of %g V\n', str1, ...
    tableInInit('vsb', table), tableInFinal('vsb', table), ...
    tableInStep('vsb', table));

strAux = sprintf('array of L (in %s m):  ', scalingString);
str1 = sprintf('%s\n%s%s', str1, strAux, ...
    printArrayToString(1/GeoScale*tableInArray('lg', table), length(strAux), ...
    lineLength, interString, blanks));

% blank lines:
str1 = sprintf('%s\n\n', str1);

strAux = 'operating point parameters:  ';
str1 = sprintf('%s\n%s%s', str1, strAux, ...
    printArrayToString(fieldnames(table.Table), length(strAux), lineLength, ...
    interString, blanks));

figure(1)
clf;
hold off;
text(0.5, 0.5, str1, 'HorizontalAlignment','center', ...
    'color', [0 0 0], 'FontName', 'Arial');
axis([0 1 0 1]);
axis off;










  
