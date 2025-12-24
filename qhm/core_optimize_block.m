switch args.solver
    
case "lbfgs"
      
    options = optimoptions('fmincon','display','iter','specifyobjectivegradient',true,'hessianapproximation',{'lbfgs',args.solver_lbfgs_m},...
             'maxiterations',50,'OptimalityTolerance',1e-18,'StepTolerance',0e-18); 

    [da] = fmincon(value_grad_fun,da,[],[],[],[],LB,UB,NONLCON,options); 
    
case "vadam"
      
    sOpt = optimset('fmin_adam');
    sOpt.MaxIter = 50;

    if (numel(args.solver_adam_lr)==1)
        adam_lr = args.solver_adam_lr;
    else
        assert(numel(args.solver_adam_lr)==args.max_iter);
        adam_lr = args.solver_adam_lr(ii);
    end

    fprintf('\t adam_lr=%g',adam_lr);

    [da] = fmin_vector_adam_simple(value_grad_fun, da, adam_lr, 0.9, 0.999, sOpt, numel(da)/f, args.solver_adam_args);
   
otherwise
    
    error('undefined solver type!\n');
    
end
    

%
[~,~,~,out_data] = value_grad_fun(da); 


u = out_data.u; 
v = out_data.v; 
%%
% number of flipped triangles. 
% flipped = doublearea([u,v], F) < 0;
% num_flipped = numel(find(flipped));
num_flipped = out_data.num_flipped; 
fprintf('\nIter %04d, flipps %04d\n', ii, num_flipped); 

record = struct();
record.u = u;
record.v = v;
record.num_flipped = num_flipped;
record.da = da;

history{ii} = record;

if args.vis
    
    % figure;
    [t,l,h] = render_mesh2([u,v],F,'EdgeColor',[0,0,0],'FaceColor',[1,1,1]); 
    axis off; 
    axis equal; 

% alpha(t,0)g
end
