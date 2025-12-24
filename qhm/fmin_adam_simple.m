function [x] = fmin_adam_simple(fun, x0, stepSize, beta1, beta2, options)

epsilon = sqrt(eps);

MaxIter = options.MaxIter;

x = x0;

n = size(x,1);

m = zeros(n,1);
v = zeros(n,1);

for iter=1:MaxIter
    
    [value, grad] = fun(x);
    
    m = beta1 * m + (1-beta1) * grad;

    v = beta2 * v + (1-beta2) * grad.^2;

    mt = m / (1-beta1^iter);
    vt = v / (1-beta2^iter);
    
    x = x - stepSize * mt ./ (sqrt(vt) + epsilon);

end
