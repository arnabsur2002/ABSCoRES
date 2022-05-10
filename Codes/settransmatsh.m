function settransmatsh(transraw, datei, namefile)
% Auxiliary function to transform transition probability matrices
% Inputs:
%     TRANSRAW  :     two dimensional transmat
%     DATEI     :     date stamp on file
%     NAMEFILE  :     name to save file ouput
% 
% Assumptions
% - transraw is in current format (nt x nj x nj)
% - number of scenarios is always the same
%
% 2022.04.25
% Alberto J. Lamadrid

nt = size(transraw, 3);
sel = 3;                            % column of transition prob matrix selected for period 1
transraw(find(transraw==0))=1e-7;   % change probability for scenarios with zero prob

fid = fopen(sprintf('%s.m', namefile), 'w');
fprintf(fid, 'function transmat = %s\n', namefile);
fprintf(fid, '%% %s\n', datei);
fprintf(fid, '%% Alberto J. Lamadrid\n');

for t = 1:nt
    if t == 1
      fprintf(fid, 'transmat%d = [\n', t-1);
    else
      fprintf(fid, 'transmat{%d} = [\n', t);
    end
    for it = 1:size(transraw, 1)
      fprintf(fid, '%12.8f', transraw(it, :, t));
      fprintf(fid, '\n');
    end
    fprintf(fid, '];\n');
    if t == 1
      fprintf(fid, 'transmat{%d} = transmat0(:, %d);\n', t, sel);
    end
end

fclose(fid);