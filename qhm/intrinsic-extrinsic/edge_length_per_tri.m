function [EL] = edge_length_per_tri(V,F)

% f = size(F,1);

rnr = @(x) sqrt(sum(x.^2,2));

EL = [rnr(V(F(:,2))- V(F(:,3))), rnr(V(F(:,3))- V(F(:,1))), rnr(V(F(:,1))- V(F(:,2)))];