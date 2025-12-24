
%% 

suitesparse_path = '~/workspace/tools/SuiteSparse-7.0.1/CHOLMOD/MATLAB/';

if ~exist(suitesparse_path, 'dir')
    error('set suitesparse_path to your ');
end

addpath(suitesparse_path); 

addpath(genpath('./')); 
% addpath(genpath('../../gptoolbox/'));  % now includes a subset of gptoolbox in the repo

if false
    path = './test_cases/Simple/';
    meshname = 'WeberZorin14_fig19'; % 'cross', 'Lshape', 'WeberZorin14_fig19', 'square_transRot90', 'square_rot180', 'WeberZorin14_fig20';
else
    path = './test_cases/Letters/';
    meshname = 'david_o_A'; % 'david_o_A', 'bunny_i_H', 'lucy_o_G'
end

folder = [path,meshname];

%% set hyperparameters and options.

args = default_args();

% set variable parameterization
args.tp_type = 'complex-plane-det1-tanh'; % 'complex-plane-det1-tanh' or 'complex-plane-det1' (sigmod or tanh). 
args.solver = 'vadam'; % 'lbfgs' or 'vadam' (vector-adam)

% lr schedule for adam: 
args.solver_adam_lr = [0.1*ones(1,5), 0.01*ones(1,5)];
args.max_iter = numel(args.solver_adam_lr); 

args.vis = false;
args.load_format = 'du2020';
args.alpha = 1; % 100; 10000;
args.F = -1;
args.energy_type = 'none'; % 'none', 'symmetric-Dirichlet-capped', 'symmetric-Dirichlet', 'ARAP', 'area-change', 'mass-spring';

% initial map: 
args.graph_tutte = true;  % true: more robust numerically but less conformal. false: less robust but more conformal.

%% run the solver

core_variational_beltrami;

%%
ut = u; vt = v;

figure('position',[0 0 1000 1000]);
[t,l,h] = render_mesh2([ut,vt],F,'EdgeColor',[0,0,0],'FaceColor',[1,1,1]);
set(t,'LineWidth',1);
axis off;
axis equal;
