function [COV,COVD] = mesh_covariance(n,F,u)

% example:
% n = size(V,1);
% u = [V(:,1);V(:,2)];
% COV = mesh_covariance(n,F,u);
%%

II = [F(:,1);F(:,2);F(:,3)];
JJ = [F(:,2);F(:,3);F(:,1)];
%%
PL = sparse(II,JJ,1,n,n);
degree = full(sum(PL+PL',2));
% interior edge is counted twice, boundary edge is counted only once. 
% degree/2 is the number of adjacent triangles per vertex. 
%%
MI = [II;II;JJ;JJ]; % note the last JJ should be added, or it will not be psd.
MJ = [II;JJ;II;JJ]; 
% indices to assemble each 3x3 small convariance matrix. 
%%
covII = [MI;MI;MI+n;MI+n];
covJJ = [MJ;MJ+n;MJ;MJ+n];

covDII = [MI;MI+n];
covDJJ = [MJ;MJ+n];
%%
COV = sparse(covII,covJJ,u(covII).*u(covJJ),2*n,2*n);
%%
COVD = sparse(covDII,covDJJ,u(covDII).*u(covDJJ),2*n,2*n);