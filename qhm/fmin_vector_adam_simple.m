function [x] = fmin_vector_adam_simple(fun, x0, stepSize, beta1, beta2, options, d, args)

% https://arxiv.org/pdf/2205.13599.pdf

epsilon = sqrt(eps);

MaxIter = options.MaxIter;

x0 = reshape(x0, [numel(x0)/d,d]);

x = x0;

n = size(x,1);

m = zeros(n,d);
v = zeros(n,1);

for iter=1:MaxIter
    
    [value, grad, ~, out_data] = fun(x(:));
    
    fprintf('vadam: iter=%04d, f=%g\n', iter, value);

    if out_data.stop_sign && args.break_with_stop 
        fprintf('VAdam: stop_sign set.\n'); 
        break;
    end
    
    grad = reshape(grad, [numel(grad)/d,d]);
    
    m = beta1 * m + (1-beta1) * grad;

    v = beta2 * v + (1-beta2) * sum(grad.^2, 2);

    mt = m / (1-beta1^iter);
    vt = v / (1-beta2^iter);
    
    x = x - stepSize * mt ./ (sqrt(vt) + epsilon);

end

x = x(:);

