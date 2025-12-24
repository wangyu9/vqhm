function [out] = example_variational_beltrami(folder, args)

if false
%% example of use

name = "Simple/square_rot180"; % "square_transRot90" "Simple/square_rot180"; "Simple/WeberZorin14_fig19" "Simple/square_rot90" % "Simple/cross" 
%name = "Letters/gargoyle_i_A"; %"Letters/hand_3_i_A"; % "Letters/gargoyle_i_A"; %	

folder = "../Locally-Injective-Mappings-Benchmark/" + name + "/";
% folder = "D:/WorkSpace/fastsol/Locally-Injective-Mappings-Benchmark/Simple/" + name + "/";

addpath(genpath('./'));
args = default_args();

%"Letters/david_o_A"; % "Letters/lucy_o_G"; % "Letters/bunny_i_G";

% args.mprint = @(fid, varargin) 

%log_path = folder + '/200x50-lbfgs-'
%diary(log_path+'log.txt')
%%
[out] = example_variational_beltrami(folder, args); 

%save(log_path+'out.mat' );

u = out.u;
v = out.v;
F = out.F;
render_mesh2([u,v],F,'EdgeColor',[0,0,0],'FaceColor',[1,1,1]);
axis off;
axis equal;

%diary off
end
%%
core_variational_beltrami;