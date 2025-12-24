% assume input: folder, args.

%

[V,F] = readOBJ([folder+"input.obj"]);
%[V,F] = readOBJ('bunny_flatten.obj');

% V = V(:,1:2);
[VT, ~] = readOBJ([folder+"result.obj"]);
% VT = V;
VT = VT(:,1:2);

for ii=1:args.sub_div_level
F_old = F;
[V,F] = upsample(V,F_old);
[VT,F2] = upsample(VT,F_old);
end


mesh = triangle_mesh(V,F);
mesh2 = triangle_mesh(VT,F);

% input: mesh, VT
NVT = - cotmatrix(VT, mesh.F)*VT; % cotmatrix is negative definite.

TVB = VT(mesh.B,:);
% target boundary vertices. 

%%
% input: V,F, TVB, NVT
%%
out = struct();
out.finish_with_local_global = false;
% mprint = args.mprint;

sparse_diag = @(ddd) sparse(1:numel(ddd), 1:numel(ddd), ddd, numel(ddd), numel(ddd));
% Example 2


tt = args.rotate; % pi/4 * 0.0;
% VT = VT * [cos(tt),-sin(tt);sin(tt),cos(tt)]';


BE = outline(F);
B = BE(:,1);

n = size(V,1);
nb = numel(B);

% 
% eq_lhs = [  sparse(1:nb, B, 1, nb, n*2);...
%             sparse(1:nb, B+n, 1, nb, n*2)];
% eq_rhs = [ VT(B,1); VT(B,2)];

if false
%
%[V2,F2,~,~] = subdivide_with_constraint(V,F,eq_lhs,eq_rhs,1);
%[VT2,~,~,~] = subdivide_with_constraint(VT,F,eq_lhs,eq_rhs,1);

%V = V2;
%F = F2;
%VT = VT2;

end



known = mesh.B;
unknown = mesh.UB;
row_slice_matrix = @(indices) sparse(indices,1:numel(indices),1,n,numel(indices));
R = row_slice_matrix(known);
S = row_slice_matrix(unknown);

BC = TVB;
f = size(F,1);

Area = doublearea(V,F)/2;
mesh.FA = []; % to avoid any code to use mesh.FA directly. 

if args.graph_tutte
    fprintf('Init with uniform graph Laplacian\n');
    mesh.GI = intrinsic_grad_equilateral(n,F);
    Area = ones(size(Area));  
end

% an nonlinear energy


%
% args.energy_type =  "none"; 
% args.energy_type = "ARAP"; % "area", "mass-spring"; %ARAP
% %"symmetric-Dirichlet"; 
% args.energy_type = "symmetric-Dirichlet"; 

if ~strcmp(args.energy_type,'none')

    [energy_ef,energy_gf, energy_hf] = intrinsic_grad_hessian(args.energy_type);

    % s_value_hess_aa = @(FW) IntrinsicHessianClass.ProjectedHessian(FW, F, V, energy_gf, energy_hf);

    e_value_grad = @(FW) deal(...
        IntrinsicHessianClass.Value(FW(:,1:2), F, V(:,1:2), energy_ef), ...
        IntrinsicHessianClass.Grad(FW(:,1:2), F, V(:,1:2), energy_gf));

end
%% complex-plane-det1: equivalent
da = [0*ones(f,1);0*ones(f,1);0*ones(f,1)];
rng(1,'philox'); % set random seed. 
da = (rand([f,3])-0.5) * 0.001;
LB = [];
UB = [];
NONLCON = [];

tp = tensor_para(Area, 2, 'complex-plane-det1'); 
tp.conj_vmap = false;

tp.alpha = args.alpha;
tp.beta = args.beta;
tp.gamma = args.gamma;

if ~strcmp(args.energy_type,'none')
    tp.energy_uv = e_value_grad;
end
%%
value_grad_fun = @(ddaa) oracle_conjugate_newton_symmetric(mesh, ddaa, [BC,R'*NVT], tp);
%% LBFGS
history = {};
num_flipped_old = f; 
tic;
%%
for ii=1:args.max_iter

if true
    
    attempt_local_global;
    
    flipped = doublearea([u,v], F) < 0;
    % number of flipped triangles. 
    num_flipped = numel(find(flipped));
    if num_flipped==0
        out.finish_with_local_global = true;
        break;
    end
    
end
%
stopping_criteria = (num_flipped==0) || ((num_flipped-num_flipped_old<=0)&&(ii>=100)); 
if stopping_criteria
    break;
end

num_flipped_old = num_flipped; 
end
t_end = toc;
fprintf('***************** Iter %04d, flipps %04d, time: %g********************\n', ii, num_flipped, t_end); 
%
%
%
%
%
%