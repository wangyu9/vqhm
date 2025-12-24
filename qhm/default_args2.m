function [args] = default_args2()

args = default_args(); 

args.sub_div_level = 0; 
args.energy_type = 'none'; 
args.vis = true; 

args.F = 0; % do not use faster routine.



args.min_iter_attempt_term = -1;

% args.energy_type = 'area';

%
args.solver = 'vadam';
% args.solver = 'lbfgs';

args.solver_adam_lr = [0.1*ones(1,15), 0.01*ones(1,8), 0.001*ones(1,2)];
% args.solver_adam_lr = 10.^linspace(-1,-3,10);
% args.solver_adam_lr = 0.1*ones(1,5); 

args.max_iter = numel(args.solver_adam_lr); 

args.tp_type = 'complex-plane-det1-tanh';
