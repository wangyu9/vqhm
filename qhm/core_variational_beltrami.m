% function [] = core_variational_beltrami(V,F,VT,args)
% assume input: folder, args.

% args = default_args();
% args.vis = true;
% folder = "../Locally-Injective-Mappings-Benchmark/Simple/square_rot180";
% core_variational_beltrami;

switch args.load_format
    
    case 'du2020'

        % for [Du et al. 2020]

        [V,F,VT] = readOBJ([folder+"/input.obj"]);

        VT = VT(:,1:2);

    case 'du2020-old'

        % for [Du et al. 2020]

        [V,F] = readOBJ([folder+"/input.obj"]);
        %[V,F] = readOBJ('bunny_flatten.obj');

        % V = V(:,1:2);
        [VT, ~] = readOBJ([folder+"/result.obj"]);
        % VT = V;
        VT = VT(:,1:2);
    
    case 'none'    
        
        
    case 'du2020-proj'

        % for [Du et al. 2020]

        [V,F] = readOBJ([folder+"/input.obj"]);
        %[V,F] = readOBJ('bunny_flatten.obj');

        VT = V(:,1:2);
        
    otherwise 
        
        error('unknown format!');

end


%%
% input: V,F, VT, for VT only the values on the boundary matter. 
% VT: target boundary vertices. 

assert(size(V,1)==size(VT,1));
assert(size(VT,2)==2);

% normalize target area. 
if args.load_normalize_target 
    [VT] = normalize_mesh(V,F,VT);
end

tt = args.rotate; % pi/4 * 0.0;
VT = VT * [cos(tt),-sin(tt);sin(tt),cos(tt)]';

for ii=1:args.sub_div_level
    F_old = F;
    [V,F] = upsample(V,F_old);
    [VT,F2] = upsample(VT,F_old);
end
%%
out = struct();
out.finish_with_local_global = false;
% mprint = args.mprint;

sparse_diag = @(ddd) sparse(1:numel(ddd), 1:numel(ddd), ddd, numel(ddd), numel(ddd));

mesh = triangle_mesh_basic(V,F,args.indEC);
% mesh = triangle_mesh(V,F);

known = mesh.IKB;
unknown = mesh.IUB;

TVB = VT(known,:);

n = size(V,1);
f = size(F,1);

row_slice_matrix = @(indices) sparse(indices,1:numel(indices),1,n,numel(indices));
R = row_slice_matrix(known);
S = row_slice_matrix(unknown);

NVT = mesh.D*(R* [TVB(:,2),-TVB(:,1)] ); 
% the old code that uses target Lap.
% do not use any more: NVT = - cotmatrix(VT, mesh.F)*VT; % cotmatrix is negative definite.


BC = TVB;
BN = R'*NVT;

% the total area is 0.5 * ( BC(:,1)' * R'*NVT(:,1) + BC(:,2)' * R'*NVT(:,2)  )
% converge to which  0.5 * sum(doublearea([u,v],F)) will, upon convergence

Area = mesh.FA;
mesh.FA = []; % to avoid any code to use mesh.FA directly. 

if ~isempty( args.epsilon_gradarea_angle)
    [mesh.GI, Area, mesh.GIS] = intrinsic_grad_area(V,F,args.epsilon_gradarea_angle);
end

if args.graph_tutte
    fprintf('Init with uniform graph Laplacian\n');
    [mesh.GI, mesh.GIS] = intrinsic_grad_equilateral(n,F);
    Area = ones(size(Area));  
end

mesh.AI = Area;

% an nonlinear energy


%
% args.energy_type =  "none"; 
% args.energy_type = "ARAP"; % "area", "mass-spring"; %ARAP
% %"symmetric-Dirichlet"; 
% args.energy_type = "symmetric-Dirichlet"; 

int_support_type = {'area', 'mass-spring', 'symmetric-Dirichlet', 'ARAP', 'symmetric-Dirichlet-capped', 'area-change'};

if ~strcmp(args.energy_type,'none')


    if any(contains(int_support_type, args.energy_type))
    
        % assert(size(V,2)==2||norm(V(:,3))==0);
        
        [energy_ef,energy_gf, energy_hf] = intrinsic_grad_hessian(args.energy_type);

        % s_value_hess_aa = @(FW) IntrinsicHessianClass.ProjectedHessian(FW, F, V, energy_gf, energy_hf);

        if isempty(args.energy_Wf)

            e_value_grad = @(FW) deal(...
                IntrinsicHessianClass.Value(FW(:,1:2), F, V, energy_ef), ...
                IntrinsicHessianClass.Grad(FW(:,1:2), F, V, energy_gf));

        else

            e_value_grad = @(FW) deal(...
                IntrinsicHessianClass.ValueWithWeights(FW(:,1:2), F, V, energy_ef, args.energy_Wf), ...
                IntrinsicHessianClass.GradWithWeights(FW(:,1:2),  F, V, energy_gf, args.energy_Wf));

        end
    elseif true
        
        assert(size(V,2)==2||norm(V(:,3))==0);
        
        energy_obj = OptimProblemArap(V(:,1:2), F, [], [], V(:,1:2));
        
        e_value_grad = @(FW) energy_obj.evaluateFunctional(FW(:), true, true, false);
    
    end

end
%% complex-plane-det1: equivalent
da = [0*ones(f,1);0*ones(f,1);0*ones(f,1)];
rng(1,'philox'); % set random seed. 
da = (rand([f,3])-0.5) * 0.001;
LB = [];
UB = [];
NONLCON = [];
%%

if isempty(args.L2_reg)

    tp = tensor_para(Area, 2, args.tp_type);  % 'complex-plane-det1'
    tp.reg_value = @(at) args.reg_value(at,Area); 
    tp.reg_grad  = @(at) args.reg_grad(at,Area); 

else
    
    tp = tensor_para(Area, 2, args.tp_type, 'RegCoeff', 0); % not using it here. 'Experimental', 

    tp.reg_grad  = @(at) args.L2_reg * [Area.*at(:,1); Area.*at(:,2); 0*at(:,3);];
    tp.reg_value = @(at) args.L2_reg * 0.5 * sum( Area .* ((at(:,1).^2) + (at(:,2).^2)) );

end




tp.conj_vmap = false; % true;

tp.alpha = args.alpha;
tp.beta = args.beta;
tp.gamma = args.gamma;



if ~strcmp(args.energy_type,'none')
    tp.energy_uv = e_value_grad;
end

reuse = reuse_struct();
reuse.fine_recorder = args.fine_recorder;
%%
value_grad_fun = @(ddaa) oracle_conjugate_newton_symmetric(mesh, ddaa, [BC,BN], tp, reuse);
%% LBFGS
history = {};
num_flipped_old = f; 
tic;
%%
for ii=1:args.max_iter
%% keep this block same as core_optimize_block.m

core_optimize_block;

%%
if (num_flipped>0) && (num_flipped<=args.F) && (ii>=args.min_iter_attempt_term) && (num_flipped-num_flipped_old<=2)
    
    attempt_local_global; % it computes num_flipped
    
    if num_flipped==0 && args.stop_when_no_flip 
        out.finish_with_local_global = true;
        break;
    end
    
end
%%
% stopping_criteria = ((num_flipped<=args.F)&&(num_flipped==0)) || ((num_flipped-num_flipped_old<=0)&&(ii>=100));
stopping_criteria = ((num_flipped<=args.F)&&(num_flipped==0));
if args.stop_when_no_flip  && stopping_criteria
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

%%
if false
%%
au = reshape(tp.s_at2au(da), [f,3]);
a11 = au(:,1);
a12 = au(:,2);
a22 = au(:,3);
mu = (a22 - a11 - 2i*a12) ./ ((1+a11).*(1+a22)-a12.^2);
uf =  abs(mu); % angle(mu); % abs(mu); % imag(mu);
render_mesh3(V,F,'EdgeColor',[0,0,0],'FaceScaleColor',uf,'ColorMap','default');
colorbar;
end
%% Other possible solvers.
if false    
%% adam
tic
sOpt = optimset('fmin_adam');
% sOpt.GradObj = 'off';
sOpt.MaxFunEvals = 20;
sOpt.MaxIter = 10;
sOpt.Display = 'iter'; % 'none'; 

[da] = fmin_adam(value_grad_fun, da, 0.01, [], [], [], 1, sOpt);
toc

[~,~,~,out_data] = value_grad_fun(da);


u = out_data.u;
v = out_data.v; 

[t,l,h] = render_mesh2([u,v],F,'EdgeColor',[0,0,0],'FaceColor',[1,1,1]);
axis off;
axis equal;
% alpha(t,0)
%%
for ii=1:20
% figure;
%% vector adam
tic
sOpt = optimset('fmin_adam');
sOpt.MaxIter = 50;

[da] = fmin_vector_adam_simple(value_grad_fun, da, 0.01, 0.9, 0.999, sOpt, numel(da)/f);
toc

[~,~,~,out_data] = value_grad_fun(da);


u = out_data.u; 
v = out_data.v; 

flipped = doublearea([u,v], F) < 0;
% number of flipped triangles. 
numel(find(flipped))

[t,l,h] = render_mesh2([u,v],F,'EdgeColor',[0,0,0],'FaceColor',[1,1,1]);
axis off;
axis equal;
% alpha(t,0)g
end

%% A local global solver: Initialization.

au = tp.s_at2au(reshape(da,[f,3]));
au = reshape(au,[f,3]);
a11 = au(:,1);
a12 = au(:,2);
a22 = au(:,3);

GI = mesh.GI;
% Area = Area;
%%
for jj=1:20
% global step:

Lw = GI'*[sparse_diag(a11),sparse_diag(a12);sparse_diag(a12),sparse_diag(a22)]*GI;
lhs = S'*Lw*S;
rhs = - S'*Lw*R * TVB;

U = lhs \ rhs;

W = S * U + R * TVB;

u = W(:,1);
v = W(:,2);

trace( W' * Lw * W / 2 )

% render meshes.

render_mesh2(W,F,'EdgeColor',[0,0,0]);
axis equal;

% local step:

if false
    % the naive slow implementation:
    for j=1:f
        Jj = (GI([j,j+f],:) * W)';
        Aj = abs(det(Jj)) * inv(Jj' * Jj) * Area(j);
        a11(j) = Aj(1,1);
        a12(j) = Aj(1,2);
        a22(j) = Aj(2,2);
    end
else
    % the fast implementation. 
GW = GI * W;
Gxu = GW(1:f,1);
Gxv = GW(1:f,2);
Gyu = GW(f+1:2*f,1);
Gyv = GW(f+1:2*f,2);
% Jj = [Gxu, Gxv; Gyu, Gyv]';

adetJ = abs( Gxu.*Gyv - Gxv.*Gyu );

a22 = (Gxu .* Gxu + Gxv .* Gxv) ./ adetJ .* Area;
a12 = - (Gxu .* Gyu + Gxv .* Gyv) ./ adetJ .* Area; 
a11 = (Gyu .* Gyu + Gyv .* Gyv) ./ adetJ .* Area;
end

flipped = doublearea([u,v], F) < 0;
% number of flipped triangles. 
numel(find(flipped))

end
%% double check
[~,g1] = oracle_conjugate_newton_symmetric(mesh, da+0.1, BC);
[~,g2] = oracle_conjugate_newton(mesh, da([1:f,2*f+1:3*f])+0.1, BC);
norm(g1([1:f,2*f+1:3*f]) -g2)
%% double check intrinsic Laplacian. 
GG = intrinsic_grad(V,F);
LL = GG' * sparse_diag([Area;Area]) * GG;
sum(sum(abs(LL + cotmatrix(V,F)))) / sum(sum(abs(LL)))
%%
end


%%
if false
%%
render_mesh3([u,v],F,'EdgeColor',[0,0,0],'FaceScaleColor', flipped, 'ColorMap', 'jet');
axis equal;
end


%
%
%
%
%

%%
% rng('shuffle'); 
% ridx = randi(2^52);
% % this is to save env to struct, slow but works.
% rname = ['this-is-a-temp-file-' , num2str(ridx) , '.mat'];

rname = args.save_file;
if false
save(rname);
end
% mat = load(rname);
% delete(rname);

% out.workspace = mat;
out.u = u;
out.v = v;
out.da = da;
out.V = V;
out.F = F;
out.num_flipped = num_flipped;
if ~isempty(args.fine_recorder)
out.history = history;
out.reuse = reuse;
end
out.tp = tp;
if false
out.mesh = mesh;
end
out.BC = BC;
out.BN = R'*NVT;

out.time = t_end;

if(out.finish_with_local_global)
    out.a11 = a11;
    out.a12 = a12;
    out.a22 = a22;
end
