clear;
num_elem = 30;
stiffness = 300;
elem_leng = 40;
boundary = [[-30, 0]; [30, 0]];
%% position
x = sym("x", [num_elem-1, 2]);
eng = cal_energy(x, num_elem, stiffness, elem_leng, boundary);
energy_fun = matlabFunction(eng, 'var', {x});
% x_g = reshape(gradient(eng)', [num_elem-1, 2]);
x_g = reshape(jacobian(eng, reshape(x, [], 1)), size(x, 1), size(x, 2));
grad_fun = matlabFunction(x_g, 'var', {x});
%% length
l = sym("l", [num_elem, 1]);
eng_len = cal_energy(l, num_elem, stiffness, elem_leng, boundary, true);
energy_len_fun = matlabFunction(eng_len, 'var', {l});
% x_g = reshape(gradient(eng)', [num_elem-1, 2]);
x_len_g = jacobian(eng_len, l);
grad_len_fun = matlabFunction(x_len_g, 'var', {l});
%% multiple exprs
expr_num = 10;
results_x_pos = zeros(num_elem-1, 2*expr_num);
results_len = zeros(num_elem, expr_num);
for expr = 1:expr_num
    %% init
    % rod_init = [25, 40];
    rod_init = (rand(num_elem-1, 2) - 0.5 ) * 100;
    % rod_init(:, 1) = (num_elem / 2 - 1) * (-1) : 1 : (num_elem / 2 - 1);
    length_init = rand(num_elem, 1) * 80;
    
    %% optimize
    rod_pos = [boundary(1, :); rod_init; boundary(2, :)];
    
    options = optimoptions('fmincon', ...
                            'GradObj', 'on', ...
                            'Display', 'iter');
    [xopt,fval,exitflag,output,lambda,grad,hessian] = ...
                fmincon({energy_fun, grad_fun}, rod_init, [], [], [], [], -100.0, 100.0, [], options);
    
    [xopt_len,fval_len,~,~,~,~,~] = ...
                fmincon({energy_len_fun, grad_len_fun}, length_init, [], [], [], [], 0.0, 100.0, [], options);
    
    %% save results
    results_x_pos(:, (expr-1)*2+1:expr*2) = xopt;
    results_len(:, expr) = xopt_len;
    %% plot
    if num_elem == 2
        figure;
        [X,Y] = meshgrid(-50:1:50,-50:1:50);
        all_energy = zeros(size(X, 1), size(X, 1));
        for i = 1:size(X, 1)
            for j = 1:size(X, 2)
    %             disp([X(i, j), Y(i, j)]);
                all_energy(i, j) = energy_fun([X(i, j), Y(i, j)]); %, num_elem, stiffness, elem_leng, boundary);
            end
        end
        surf(X, Y, all_energy);
        hold on;
        % plot3(round(rod_pos(:, 1))', round(rod_pos(:, 2))', all_energy(round(rod_pos(:, 1)) + 50, round(rod_pos(:, 2)) + 50), '-ok');
        % hold on;
        plot3(round(rod_init(1)), round(rod_init(2)), all_energy(round(rod_init(2)) + 50, round(rod_init(1)) + 50), '.', 'Color', 'b', 'MarkerSize', 20);
        hold on;
        grad_init = grad_fun(rod_init);
        end_grad = rod_init + grad_init;
        end_grad = end_grad / norm(end_grad);
        plot3([rod_init(1) end_grad(1)], ...
            [rod_init(2) end_grad(2)], ...
            [all_energy(round(rod_init(2)) + 50, round(rod_init(1)) + 50) ...
            all_energy(round(end_grad(2)) + 50, round(end_grad(1)) + 50)], 'r-', 'LineWidth', 2);
        hold on;
        plot3(round(xopt(1)), round(xopt(2)), all_energy(round(xopt(2)) + 50, round(xopt(1)) + 50), '.', 'Color', 'r', 'MarkerSize', 20);
    end
end
% figure;
% rod_pos_final = [boundary(1, :); xopt; boundary(2, :)];
% plot(rod_pos_final(:, 1), rod_pos_final(:, 2), '-ok');
%% functions
function l2_len =  cal_length(x1, x2)
    l2_len = sqrt((x1(1) - x2(1))^2 + (x1(2) - x2(2))^2);
end
function energy = cal_energy(x_input, num_elem, stiffness, elem_leng, boundary, using_len, varargin)
    if nargin < 6
        using_len = false;
        rod_pos_tmp = [boundary(1, :); x_input; boundary(2, :)];
    end
    energy = 0;
    for i = 1:num_elem
        if using_len
            energy = energy + 1/2 * stiffness * (x_input(i) - elem_leng)^2;
        else
%             disp(cal_length(rod_pos_tmp(i, :), rod_pos_tmp(i+1, :)));
            energy = energy + 1/2 * stiffness * (cal_length(rod_pos_tmp(i, :), rod_pos_tmp(i+1, :)) - elem_leng)^2;
        end
    end
end