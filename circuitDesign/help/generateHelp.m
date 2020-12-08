function generateHelp(helpfile, directory, title)
%
%
%  (c) IMEC, 2004
%  IMEC confidential 
%

w = what(directory);
fid = fopen(helpfile, 'w');
fprintf(fid, '\\section{%s}\n', title);
for i = 1:length(w.m)
  fileName = char(w.m(i));
  fileNameLength = length(fileName);
  if ((fileNameLength > 2) & strcmp(fileName(fileNameLength - 1), '.') & ...
	strcmpi(fileName(fileNameLength), 'm'))  % i.e. we have an m-file
    fileNameNoExtension = fileName(1:fileNameLength-2);
    helpString = help(fileNameNoExtension);
    if helpString
      fprintf(fid, '\\subsection{%s}\n', fileNameNoExtension);
      fprintf(fid, '\\label{sec:%s}\n', fileNameNoExtension);
      fprintf(fid, '\\index{%s@\\textsf{\\em %s}}\n', fileNameNoExtension, ...
          fileNameNoExtension);
      fprintf(fid, '\\begin{verbatim}\n');
      fprintf(fid, '%s\n', helpString);
      fprintf(fid, '\\end{verbatim}\n\n');
      fprintf(fid, '\\newpage\n');
    end
  end
end
fprintf(fid, '\\newpage\n');
fclose(fid);
