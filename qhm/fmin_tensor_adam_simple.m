function [x] = fmin_tensor_adam_simple(fun, x0, stepSize, beta1, beta2, options, d)

% https://arxiv.org/pdf/2205.13599.pdf

epsilon = sqrt(eps);

MaxIter = options.MaxIter;

x0 = reshape(x0, [numel(x0)/d,d]);

x = x0;

n = size(x,1);

m = zeros(n,d);
v = zeros(n,1);

sum_fun = @(gg) gg(:,1).^2 + 2 * gg(:,2).^2 + gg(:,3).^2;

for iter=1:MaxIter
    
    [value, grad] = fun(x(:));
    
    grad = reshape(grad, [numel(grad)/d,d]);
    
    m = beta1 * m + (1-beta1) * grad;

    v = beta2 * v + (1-beta2) * sum_fun(grad);

    mt = m / (1-beta1^iter);
    vt = v / (1-beta2^iter);
    
    x = x - stepSize * mt ./ (sqrt(vt) + epsilon);

end

x = x(:);
