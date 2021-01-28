function generateMasterDocTex(destDir)
%GENERATEMASTERDOCTEX generation of LaTeX code from the help files of all
%  functions related to circuit design.
%  generateMasterDocTex(DESTDIR) generates LaTeX files in the subdirectory 
%  DESTDIR.
%
%  Example:
%
%  generateMasterDocTex('d:\piet\matlab\circuitDesign\help');
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

generateHelp(strcat(destDir, '\tech.tex'), ...
    'd:\piet\matlab\circuitDesign\tech', ...
'Matlab files related to technology characterization');

generateHelp(strcat(destDir, '\mos.tex'), ...
    'd:\piet\matlab\circuitDesign\mos', ... 
'Matlab files related to the MOS transistor');

generateHelp(strcat(destDir, '\table.tex'), ...
    'd:\piet\matlab\circuitDesign\table', ... 
'Matlab files related to tables of a MOS transistor');

generateHelp(strcat(destDir, '\cir.tex'), ...
    'd:\piet\matlab\circuitDesign\cir', ... 
'Matlab files related to the circuit level');

generateHelp(strcat(destDir, '\noise.tex'), ...
    'd:\piet\matlab\circuitDesign\cir\noise', ... 
'Matlab files related to noise calculations');

