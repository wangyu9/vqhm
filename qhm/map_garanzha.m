function [] = map_garanzha(V,F,VT,folder,args)

% see also. total_lifted_content

folder = char(folder); 

% mkdir(folder);


input_mesh = [folder,'/input.obj']; 

writeOBJ(input_mesh, V, F, VT(:,1:2));

% linux 
% myPath = getenv("PATH");
myPath = '/home/wangyu/workspace/tools/ws_moveit2/install/moveit_core/bin:/opt/ros/humble/bin:/home/wangyu/workspace/tools/anaconda3/bin:/home/wangyu/workspace/tools/anaconda3/condabin:/home/wangyu/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/snap/bin'; 
%  system(['export PATH=' myPath ' ; which gcc']);

cpath = ['export PATH=' myPath ' ; ~/workspace/projects/qhm/invertible-maps/cpp/build/untangle2d ../input.obj --save_iter 1 ' char(args)]; 

cpath = ['~/workspace/projects/qhm/invertible-maps/cpp/build/untangle2d ../input.obj --save_iter 1 ' char(args)]; 

command1 = [cpath ' '  ];


fprintf('%s\n',command1); 
% [status, result] = system( command1 );
% 
% if(status~=0)
%     error('');
% end
% disp(result);
