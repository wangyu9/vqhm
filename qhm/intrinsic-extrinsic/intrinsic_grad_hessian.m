function [e,g,h] = intrinsic_grad_hessian(energy) 
%%
syms aa1;
syms bb1;
syms cc1;

syms aa0;
syms bb0;
syms cc0;

assume(aa1,'real')
assume(bb1,'real')
assume(cc1,'real')

assume(aa0,'real')
assume(bb0,'real')
assume(cc0,'real')
%

R = sqrt( 2*aa0.*bb0+2*aa0.*cc0+2*bb0.*cc0 - aa0.^2 - bb0.^2 - cc0.^2  );

A = ( (-aa0 + bb0 + cc0).* aa1 + (aa0 - bb0 + cc0).* bb1 + (aa0 + bb0 - cc0).* cc1) ./ (R.^2);

ysq = ((bb0.*cc0).*aa1.^2 + (aa0.*cc0).*bb1.^2 + (aa0.*bb0).*cc1.^2 ...
    + (-aa0 - bb0 + cc0).* cc0 .* bb1 .* aa1 ...
    + (-aa0 + bb0 - cc0).* bb0 .* aa1 .* cc1 ...
    + ( aa0 - bb0 - cc0).* aa0 .* cc1 .* bb1 ...
 ) ./ (R.^4);

S = [ sqrt(A + 2*sqrt(ysq)), sqrt(A - 2*sqrt(ysq))];
% singular values.


switch energy
    % https://en.wikipedia.org/wiki/Heron%27s_formula
    case "area"
        e = 1/4*sqrt(2*(aa1.*bb1+bb1.*cc1+cc1.*aa1)-(aa1.^2+bb1.^2+cc1.^2)) ./ ...
            (1/4*sqrt(2*(aa0.*bb0+bb0.*cc0+cc0.*aa0)-(aa0.^2+bb0.^2+cc0.^2)));
    case "area-r"
        epsilon = 1e-18;
        e = 1/4*sqrt( epsilon + 2*(aa1.*bb1+bb1.*cc1+cc1.*aa1)-(aa1.^2+bb1.^2+cc1.^2)) ./ ...
            (1/4*sqrt(2*(aa0.*bb0+bb0.*cc0+cc0.*aa0)-(aa0.^2+bb0.^2+cc0.^2)));
    case "mass-spring"
        e = 1e-3 * ( (sqrt(aa1)-sqrt(aa0)).^2 + (sqrt(bb1)-sqrt(bb0)).^2 +  (sqrt(cc1)-sqrt(cc0)).^2  )...
            ./ (aa0 + bb0 + cc0);
    case "symmetric-Dirichlet"
        e = (A * 2 + A * 2 ./ (A.^2 - 4*ysq))/2;
    case "symmetric-Dirichlet-capped"
        epsilon = 1e-6;
        e = 1e-6*(A * 2 + A * 2 ./ (A.^2 - 4*ysq + epsilon))/2;
    case "symmetric-Dirichlet-equiv"
        % sym Dirichlet in an equivalent form. 
        e = (S(1).^2+S(2).^2 + 1./(S(1).^2)+1./(S(2).^2))/2; 
    case "ARAP"
        e = 1e-3*(2*A + 2 - 2*(S(1)+S(2))); % ARAP
    case "ARAP-equiv"
        e = 1e-3*( (S(1)-1).^2+(S(2)-1).^2 ); % ARAP
    case "area-change"
        e = 1e-3 * (1 - sqrt(2*(aa1.*bb1+bb1.*cc1+cc1.*aa1)-(aa1.^2+bb1.^2+cc1.^2)) ./ ...
            (sqrt(2*(aa0.*bb0+bb0.*cc0+cc0.*aa0)-(aa0.^2+bb0.^2+cc0.^2))) ).^2;
    otherwise
        error(['unsupport energy']);
end


H2 = ...
[diff(e,aa1,aa1), diff(e,aa1,bb1), diff(e,aa1,cc1);
 diff(e,bb1,aa1), diff(e,bb1,bb1), diff(e,bb1,cc1);
 diff(e,cc1,aa1), diff(e,cc1,bb1), diff(e,cc1,cc1);];

g = matlabFunction([diff(e,aa1), diff(e,bb1), diff(e,cc1)]');

h = matlabFunction(H2);

e = matlabFunction(e);