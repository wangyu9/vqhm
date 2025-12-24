function [args] = default_args()

args = struct();
args.solver = 'lbfgs'; % other choice is 'vadam'. 

args.solver_adam_lr = 0.01; % the default lr if using adam or vadam 
args.solver_lbfgs_m = 10;

args.solver_adam_args = struct(); 
args.solver_adam_args.break_with_stop = false; 
% weather to stop when out_data.stop_sign is true. 

args.sub_div_level = 0;
args.energy_type = "none";
args.energy_Wf = [];

args.vis = false;
args.max_iter = 20;
args.rotate = 0;

args.alpha = 1;
args.beta = 0;
args.gamma = 0;

args.F = 10; 
% less equal than the number invokes a faster routine.

args.load_format = 'du2020'; % 'none'; 
args.load_normalize_target = true;

args.graph_tutte = false;

args.save_file = [];


args.reg_value = @(at,Area) 0;
args.reg_grad  = @(at,Area) 0;

args.tp_type = 'complex-plane-det1';

args.min_iter_attempt_term = 16;

args.fine_recorder = [];

args.L2_reg = [];

args.stop_when_no_flip = true; % whether to finish iterations when there are no more flip. 

args.epsilon_gradarea_angle = [];



args.indEC = []; 