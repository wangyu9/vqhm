function [MVL] = signed_mean_value_laplacian(V,F)
%%

n = size(V,1);
f = size(F,1);

mesh = struct();
mesh.V = V;
mesh.F = F;

EL = edge_lengths(V,F); 
mesh.EL = EL;

v12 = V(F(:,2),:) - V(F(:,1),:); 
v13 = V(F(:,3),:) - V(F(:,1),:); 

v23 = V(F(:,3),:) - V(F(:,2),:); 
v21 = V(F(:,1),:) - V(F(:,2),:); 

v31 = V(F(:,1),:) - V(F(:,3),:); 
v32 = V(F(:,2),:) - V(F(:,3),:); 
        
dot23 = sum(v12.*v13, 2);
dot31 = sum(v23.*v21, 2);
dot12 = sum(v31.*v32, 2);

mesh.dot23 = dot23;
mesh.dot31 = dot31;
mesh.dot12 = dot12;

% 2D cross
cross23 = v12(:,1) .* v13(:,2) - v12(:,2) .* v13(:,1); 
cross31 = v23(:,1) .* v21(:,2) - v23(:,2) .* v21(:,1); 
cross12 = v31(:,1) .* v32(:,2) - v31(:,2) .* v32(:,1); 

mesh.cross23 = cross23;
mesh.cross31 = cross31;
mesh.cross12 = cross12;

% https://en.wikipedia.org/wiki/Tangent_half-angle_formula
% the *signed* tan angle 213 (right hand rule for triangle 123 for positive sign)
% tan(theta/2) = sin(theta) / (1+cos(theta))
stan1 = cross23 ./ (dot23 + EL(:,2).*EL(:,3) ); 
stan2 = cross31 ./ (dot31 + EL(:,3).*EL(:,1) );
stan3 = cross12 ./ (dot12 + EL(:,1).*EL(:,2) );

% mesh.stan1 = stan1;
% mesh.stan2 = stan2;
% mesh.stan3 = stan3;

stan = [stan1,stan2,stan3];
mesh.stan = stan;


VV = [stan./EL(:,[3,1,2]), stan./EL(:,[2,3,1])];
II = [F,F];
JJ = [F(:,[2,3,1]),F(:,[3,1,2])];

% MVL = sparse(II,JJ,VV,n,n);
% MVL = MVL - diag(sum(MVL,2));

% equivalently:
MVL = sparse(II,JJ,VV,n,n) - sparse(II,II,VV,n,n);

% sum(sum(abs( MVL-mean_value_laplacian(V,F) ))) 
% gives a number around 0.  
end
