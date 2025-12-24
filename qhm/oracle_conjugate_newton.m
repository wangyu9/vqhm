function [E, g, Hess, out] = oracle_conjugate_newton(mesh, da, BC)
%%
sparse_diag = @(ddd) sparse(1:numel(ddd), 1:numel(ddd), ddd, numel(ddd), numel(ddd));
V = mesh.V;
F = mesh.F;

%%
BE = outline(F);
B = BE(:,1);

n = size(V,1);
nb = numel(B);

% eq_lhs = [  sparse(1:nb, B, 1, nb, n*2);...
%             sparse(1:nb, B+n, 1, nb, n*2)];
% eq_rhs = [ VT(B,1); VT(B,2)];
% 
% if false
% %
% [V2,F2,~,~] = subdivide_with_constraint(V,F,eq_lhs,eq_rhs,1);
% [VT2,~,~,~] = subdivide_with_constraint(VT,F,eq_lhs,eq_rhs,1);
% 
% V = V2;
% F = F2;
% VT = VT2;
% end

BE = outline(F);
B = BE(:,1);

n = size(V,1);
nb = numel(B);

para = struct();
% coefficients for constraints
para.CC = 10;
% coefficients for energy
para.CE = 1; 
% coefficients for regularizer.
para.CR = 0; %

%
%mesh = triangle_mesh(V,F);
%mesh2 = triangle_mesh(VT,F);

d01 = mesh.DEC.d01;

G = mesh.G;

n = mesh.n;
ne = mesh.ne;
f = mesh.f;
ft = ones(1,f);

known = mesh.B;
unknown = mesh.UB;
row_slice_matrix = @(indices) sparse(indices,1:numel(indices),1,n,numel(indices));
R = row_slice_matrix(known);
S = row_slice_matrix(unknown);


Area = doublearea(V,F)/2;

%%


DH = [sparse_diag(zeros(f,1)), -sparse_diag(ones(f,1)); ...
    sparse_diag(ones(f,1)), sparse_diag(zeros(f,1))];

u = V(:,1);
v = V(:,2);
u(known,:) = BC(:,1);
v(known,:) = BC(:,2);
%%
%ua = da .* [Area;Area];
%va = 1./da .* [Area;Area];

A_u = (G)' * sparse_diag(da) * G;
A_v = (DH * G)' * sparse_diag(1./da) * DH * G;

u(unknown,:) = - A_u(unknown,unknown) \ (A_u(unknown,known) * u(known,:));
v(unknown,:) = - A_v(unknown,unknown) \ (A_v(unknown,known) * v(known,:));


g_u = 0.5 * sparse_diag( G * u) * G * u -  (sparse_diag(G*u) * G * S) * ((S'*A_u*S)\(S'*A_u*u));

%
g_v = - sparse_diag(1./da.^2) * ( ...
        0.5 * sparse_diag(DH*G * v) * DH*G * v ...
        -  sparse_diag(DH*G*v) *DH * G * S * ((S'*A_v*S)\(S'*A_v*v)) ...
        );

%
Hess_v = sparse_diag( (DH * G * v).^2 ./ da.^3 ) ...
    + sparse_diag(1./da.^2) * sparse_diag(DH * G * v) * ...
    DH * G * S * inv( S' * (DH * G)' * sparse_diag(1./da) * DH * G * S ) * S' * G' * DH' * ...
     sparse_diag(DH * G * v) * sparse_diag(1./da.^2);

Hess_u = - sparse_diag(G * u) * ...
     G * S * inv( S' * (G)' * sparse_diag(da) * G * S ) * S' * G' * ...
     sparse_diag( G * u);

E_u = 0.5 * u' * A_u * u;
E_v = 0.5 * v' * A_v * v;

E = E_u + E_v;

%
% render_mesh2([u,v],F,'EdgeColor',[0,0,0]);
% axis equal;

% dda = 1e-7 * normrnd(0,1,size(da));
% da = da_old + dda;
%%

% da = da - 0.001 *  (0.001*speye(numel(da)) + Hess_v) \ (g_u + g_v);
% da = da - 0.01 * (g_u + g_v);

%%
g = g_u + g_v;
Hess = Hess_u + Hess_v;


out = struct();
out.u = u;
out.v = v;
%%

