function [D] = mesh_dirac(n,F)

% x1 y2 + x2 y3 + x3 y1 - (x2 y1 + x3 y2 + x1 y3)
Dh = sparse( F(:,1), F(:,2), 1, n, n) ...
   + sparse( F(:,2), F(:,3), 1, n, n) ...
   + sparse( F(:,3), F(:,1), 1, n, n);

D =  0.5*(Dh - Dh');
