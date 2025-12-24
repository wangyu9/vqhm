function [VR,history] = total_lifted_content(V,F,VT,folder)

folder = char(folder); 

mkdir(folder);

folder = [pwd,'/',folder,'/'];

% folder = [folder,'/'];

IO_path = '../lifting_simplices_to_find_injectivity/IO/'; 

if ismac
    python = '/Users/yu/opt/anaconda3/bin/python'; 
else
    python = '/home/wangyu/anaconda3/bin/python';
end

convert_input_2D = [python ' ' IO_path 'convert_input_2D.py'];

extract_boundary_vert = [python ' ' IO_path 'extract_boundary_vert.py'];

get_result_mesh = [python ' ' IO_path 'get_result_mesh.py'];

input_mesh = [folder,'/input.obj']; % /Users/yu/workspace/matlab/qhm/AQP/temp/current

handle = [folder,'/handles.txt']; 

% %  % handle = '../Locally-Injective-Mappings-Benchmark/Simple/square_rot180/handles.txt';

% 'temp/current/';

DuInputFormat = [folder,'DuInputFormat'];

DuOutput = [folder,'DuOutput'];

result_mesh = [folder,'result.obj'];

% writeOBJ_safe(input_mesh, V, F, VT(:,1:2));
writeOBJ(input_mesh, V, F, VT(:,1:2));

% if isfile([folder,'/input-meshlab.obj'])
%     input_mesh = [folder,'/input-meshlab.obj'];
% end

% command0 = [extract_boundary_vert ' ' input_mesh ' ' handle];
% 
% fprintf('%s\n',command0); 
% [status, result] = system( command0 );
% 

[BE, ~, ~] = circular_laplacian(size(V,1),F);
B = sort(BE(:,1));

dlmwrite(handle,B-1,'delimiter','\n', 'precision', '%d'); 
% dlmwrite will have an extra '\n' at the end, which seems not to be a problem


command1 = [convert_input_2D ' ' input_mesh ' '  handle ' ' DuInputFormat];


fprintf('%s\n',command1); 
[status, result] = system( command1 );

if(status~=0)
    error('');
end

disp(result);

% ./extract_boundary_vert.py [inputMeshFile] [outputHandleFile] 
% ./convert_input_2D.py [inputObjFile] [handleFile] [outFile]

% ../lifting_simplices_to_find_injectivity/IO/extract_boundary_vert.py $WD/input.obj temp/temp_handle.txt


% export fjPN=../Total-Lifted-Content-PN/build/findInjective_PN

% $fjPN  $WD/DuInputFormat my_PN_solver_options0 $WD/Du-PN-result

fjPN = '../Total-Lifted-Content-PN/build/findInjective_PN'; 

config = [folder,'/solver_option'];
% config = '../scripts/my_PN_solver_options0'; 

write_Du_2020_option(config);


command2 = [fjPN ' ' DuInputFormat ' ' config ' ' DuOutput]; 


fprintf('%s\n',command2); 
[status, result] = system( command2 );

if(status~=0)
    error('');
end

disp(result);

command3 = [ get_result_mesh ' ' DuInputFormat ' ' DuOutput ' ' result_mesh];


fprintf('%s\n',command3); 
[status, result] = system( command3 );

if(status~=0)
    error('');
end

disp(result);
% ./get_result_mesh.py [inputFile] [resultFile] [outFile]



VR = readOBJ(result_mesh);

n = size(V,1);
history = [];

for ii=0:100000
    iter_file = [folder,'/DuOutput_vert_iter', num2str(ii)];

    if isfile(iter_file)
        VDD = readDuMat(iter_file, n);

        history(ii+1).UV = VDD;

    else
        break;
    end
end
% 
%
% figure('position',[0 0 1000 1000]);
% draw_uv_with_flips2(VDD(:,1),VDD(:,2),F,[]);


