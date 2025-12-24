function [E, g, Hess, out] = oracle_anisotropic_laplace_joint_grad(mesh, at, uv,BC, tp)
%%
sparse_diag = @(ddd) sparse(1:numel(ddd), 1:numel(ddd), ddd, numel(ddd), numel(ddd));
V = mesh.V;
F = mesh.F;

dim = 2;
%%
BE = outline(F);
B = BE(:,1);

n = size(V,1);
nb = numel(B);


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


% assert(tp.conj_vmap==false);

if tp.conj_vmap

    DH = [sparse_diag(zeros(f,1)), -sparse_diag(ones(f,1)); ...
        sparse_diag(ones(f,1)), sparse_diag(zeros(f,1))];

else

    DH = speye(2*f);
    
end
    

u = uv(:,1);
v = uv(:,2);

%% setup parameterizations of A

% [s_at2au,s_at2au_mul,s_pdapdt_lmul] = para_fun_old(Area,dim,'diag');

% [s_at2au_v,s_at2au_mul_v,s_pdapdt_lmul_v] = para_fun_old(Area,dim,'invdiag');

if true

para_type = 'diag'; % 'diag-no-mass';

%[s_at2au,s_pdapdt_lmul,s_at2au_v,s_pdapdt_lmul_v] = ...
%    para_fun(Area,dim,para_type);

s_at2au = tp.s_at2au;
s_pdapdt_lmul = tp.s_pdapdt_lmul;
s_at2au_v = tp.s_at2au_v;
s_pdapdt_lmul_v = tp.s_pdapdt_lmul_v;


end

%%
%ua = da .* [Area;Area];
%va = 1./da .* [Area;Area];

au = s_at2au(reshape(at,[f,3]));

AA = zeros(f,2,2);
IA = zeros(f,2,2);

a11 = au(1:f);
a12 = au(f+1:2*f);
a22 = au(2*f+1:end);

AA(:,1,1) = a11(:);
AA(:,1,2) = a12(:);
AA(:,2,1) = a12(:);
AA(:,2,2) = a22(:);

av = s_at2au_v(reshape(at,[f,3]));

av11 = av(1:f);
av12 = av(f+1:2*f);
av22 = av(2*f+1:end);

IA(:,1,1) = av11(:);
IA(:,1,2) = av12(:);
IA(:,2,1) = av12(:);
IA(:,2,2) = av22(:);

A_u = (G)' * [sparse_diag(AA(:,1,1)), sparse_diag(AA(:,1,2));...
              sparse_diag(AA(:,2,1)), sparse_diag(AA(:,2,2));] * G;
A_v = (DH * G)' * [sparse_diag(IA(:,1,1)), sparse_diag(IA(:,1,2));...
                   sparse_diag(IA(:,2,1)), sparse_diag(IA(:,2,2));] * DH * G;


% the span op
d = 2;
span = @(xx) symmetric_tensor_span(xx,d);

g_u = A_u * u;

g_v = A_v * v;

g_at_u = s_pdapdt_lmul( reshape(at,[f,3]), ...
    0.5 * span( G * u)' * G * u  ...
    );

%
g_at_v = s_pdapdt_lmul_v ( reshape(at,[f,3]), ...
        0.5 * span(DH*G * v)' * DH*G * v ...
        );

g_at = g_at_u + g_at_v + tp.reg_grad(at);
    
%
% Hess_v = sparse_diag( (DH * G * v).^2 ./ da.^3 ) ...
%     + sparse_diag(1./da.^2) * sparse_diag(DH * G * v) * ...
%     DH * G * S * inv( S' * (DH * G)' * sparse_diag(1./da) * DH * G * S ) * S' * G' * DH' * ...
%      sparse_diag(DH * G * v) * sparse_diag(1./da.^2);
% 
% Hess_u = - sparse_diag(G * u) * ...
%      G * S * inv( S' * (G)' * sparse_diag(da) * G * S ) * S' * G' * ...
%      sparse_diag( G * u);

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
% g = [zeros(size(g_at)); g_u; g_v;];
g = [g_at; g_u; g_v;];

g = [g;tp.append_grad];

% Hess = Hess_u + Hess_v;
Hess = [];


out = struct();
out.u = u;
out.v = v;

end
%%