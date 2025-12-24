function [reg_value, reg_grad] = reg_value_grad(reg_type, weight)

    %weight = 1e-3;
    %weight = 1e-7;
    energy_radical = @(at) weight * 4 * ( cosh(sqrt( at(:,1).^2 + at(:,2).^2  )).^2 - 1 ); 
       
    reg_value = @(at,AI) sum( AI .* energy_radical(at)  );
    
    reg_grad  = @(at,AI) weight * reshape( [grad_radical(at) .* AI, 0*at(:,3)] ,[numel(at),1]);

end


function [g] = grad_radical(at)
 
    % only at(:,1:2) will be used. 

    r = sqrt( at(:,1).^2 + at(:,2).^2  ); 

    g = (4 * 2) * cosh(r) .* sinh(r) .* [at(:,1), at(:,2)] ./ r; 

end
