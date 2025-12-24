function [] = write_Du_2020_option(file,varargin)


fileID = fopen(file,'w');

form = 'Tutte';
alphaRatio = 1.e-6; 
alpha = -1;
ftol_abs = 1.e-8;
ftol_rel = 1.e-8; 
xtol_abs = 1.e-8; 
xtol_rel = 1.e-8; 
gtol_abs = 1.e-8; 
algorithm = 'Projected_Newton'; 
maxeval = 75; % 10000;
stopCode = 'none'; % 'all_good'; 

% 'record'

vert = 1; 
energy = 1; 
minArea = 1; 
gradient = 0; 
gNorm = 0; 
searchDirection = 0;
searchNorm = 1;
stepSize = 0; 
stepNorm = 0;

% 'save'

save_vert = 1;



fprintf(fileID,'form\n%s\n',form);
fprintf(fileID, 'alphaRatio\n%g\n', alphaRatio);
fprintf(fileID, 'alpha\n%g\n', alpha);
fprintf(fileID, 'ftol_abs\n%g\n', ftol_abs);
fprintf(fileID, 'ftol_rel\n%g\n', ftol_rel);
fprintf(fileID, 'xtol_abs\n%g\n', xtol_abs);
fprintf(fileID, 'xtol_rel\n%g\n', xtol_rel);
fprintf(fileID, 'gtol_abs\n%g\n', gtol_abs);
fprintf(fileID, 'algorithm\n%s\n', algorithm);
fprintf(fileID, 'maxeval\n%d\n', maxeval);
fprintf(fileID, 'stopCode\n%s\n', stopCode);

fprintf(fileID, 'record\n');

fprintf(fileID, 'vert\t%d\n', vert);
fprintf(fileID, 'energy\t%d\n', energy);
fprintf(fileID, 'minArea\t%d\n', minArea);
fprintf(fileID, 'gradient\t%d\n', gradient);
fprintf(fileID, 'gNorm\t%d\n', gNorm);
fprintf(fileID, 'searchDirection\t%d\n', searchDirection);
fprintf(fileID, 'searchNorm\t%d\n', searchNorm);
fprintf(fileID, 'stepSize\t%d\n', stepSize);
fprintf(fileID, 'stepNorm\t%d\n', stepNorm);

fprintf(fileID, 'save\n');

fprintf(fileID, 'vert\t%d', save_vert);



fclose(fileID);
