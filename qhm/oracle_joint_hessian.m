function [E, g, Hess, out] = oracle_joint_hessian(mesh, ddaa, BCBN, tp)

sparse_diag = @(ddd) sparse(1:numel(ddd), 1:numel(ddd), ddd, numel(ddd), numel(ddd));
V = mesh.V;
F = mesh.F;

n = size(V,1);
dim = 2;

BC = BCBN(:,1:2);
BN = BCBN(:,3:4);

known = mesh.IKB;
unknown = mesh.IUB;

nu = numel(unknown);

u = zeros(n,1);
v = zeros(n,1);

u(known) = BC(:,1);
v(known) = BC(:,2);


u(unknown) = ddaa(1:nu);
v(unknown) = ddaa(nu+1:2*nu);
at = ddaa(2*nu+1:end);


f = size(F,1);

% if size(at,1)~=f
%    at = reshape(at,[f,numel(at)/f]);
% end

assert(numel(at)==f*2);
at = reshape(at,[f,2]);
at = [at, zeros(f,1)];


out = struct();

%%

G = mesh.GI;

known = mesh.IKB;
unknown = mesh.IUB;

row_slice_matrix = @(indices) sparse(indices,1:numel(indices),1,n,numel(indices));
R = row_slice_matrix(known);
S = row_slice_matrix(unknown);


% the span op
d = 2;
% this is slow
span = @(xx) symmetric_tensor_span(xx,d);
span_dot = @(xx,yy) symmetric_tensor_span_dot(xx,yy,d); 

%%
old_imp = false; 

s_at2au = tp.s_at2au;
s_pdapdt_lmul = tp.s_pdapdt_lmul;


au = s_at2au(reshape(at,[f,numel(at)/f]));

AA = zeros(f,2,2);

a11 = au(1:f);
a12 = au(f+1:2*f);
a22 = au(2*f+1:end);

AA(:,1,1) = a11(:);
AA(:,1,2) = a12(:);
AA(:,2,1) = a12(:);
AA(:,2,2) = a22(:);


DiagA = [sparse_diag(AA(:,1,1)), sparse_diag(AA(:,1,2));...
              sparse_diag(AA(:,2,1)), sparse_diag(AA(:,2,2));];

A = (G)' * DiagA * G;

Auu = A(unknown,unknown);

Gx = G(1:f,:);
Gy = G(f+1:2*f,:);
%%

E_u = 0.5 * u' * A * u - 0.5 * BC(:,1)' * BN(:,1);
E_v = 0.5 * v' * A * v - 0.5 * BC(:,2)' * BN(:,2);

E = E_u + E_v;

%%

Gu = G*u;
Gv = G*v;

Gxu = Gx * u;
Gxv = Gx * v;

Gyu = Gy * u;
Gyv = Gy * v;

%%%%%%%%%%

h11 = 0.5 * (Gxu.^2+Gxv.^2) .* tp.h_dA1_dP11(at) ...
    +   (Gxu.*Gyu+Gxv.*Gyv) .* tp.h_dA2_dP11(at) ...
    + 0.5 * (Gyu.^2+Gyv.^2) .* tp.h_dA3_dP11(at);

h12 = 0.5 * (Gxu.^2+Gxv.^2) .* tp.h_dA1_dP12(at) ...
    +   (Gxu.*Gyu+Gxv.*Gyv) .* tp.h_dA2_dP12(at) ...
    + 0.5 * (Gyu.^2+Gyv.^2) .* tp.h_dA3_dP12(at);

h22 = 0.5 * (Gxu.^2+Gxv.^2) .* tp.h_dA1_dP22(at) ...
    +   (Gxu.*Gyu+Gxv.*Gyv) .* tp.h_dA2_dP22(at) ...
    + 0.5 * (Gyu.^2+Gyv.^2) .* tp.h_dA3_dP22(at);

H = [sparse_diag(h11), sparse_diag(h12); sparse_diag(h12), sparse_diag(h22) ];

SS = S;
% or set: SS = speye(n);

Zss = sparse([],[],[],size(SS,2),size(SS,2));

% pdapdt is (2f x 3f)
pdapdt = tp.s_pdapdt(at)';

Fu = SS' * G' * span(Gu) * pdapdt';
Fv = SS' * G' * span(Gv) * pdapdt';

g = [SS'*A*u; SS'*A*v; pdapdt*0.5*(span(Gu)'*(Gu)+span(Gv)'*(Gv))];

Hess = [SS'*A*SS, Zss, Fu; ...
        Zss, SS'*A*SS, Fv; ...
        Fu', Fv', H];

out.u = u;
out.v = v;

%%
flipped = doublearea([u,v], F) < 0;
% number of flipped triangles. 
num_flipped = numel(find(flipped));
fprintf('********** flipps %04d\n', num_flipped); 