function [S] = singular_values_from_edge_lengths(E,E0)


a0 = E0(:,1).^2;
b0 = E0(:,2).^2;
c0 = E0(:,3).^2;

a1 = E(:,1).^2;
b1 = E(:,2).^2;
c1 = E(:,3).^2;

R = sqrt( 2*a0.*b0+2*a0.*c0+2*b0.*c0 - a0.^2 - b0.^2 - c0.^2  );

A = ( (-a0 + b0 + c0).* a1 + (a0 - b0 + c0).* b1 + (a0 + b0 - c0).* c1) ./ (R.^2);

ysq = ((b0.*c0).*a1.^2 + (a0.*c0).*b1.^2 + (a0.*b0).*c1.^2 ...
    + (-a0 - b0 + c0).* c0 .* b1 .* a1 ...
    + (-a0 + b0 - c0).* b0 .* a1 .* c1 ...
    + ( a0 - b0 - c0).* a0 .* c1 .* b1 ...
 ) ./ (R.^4);

S = [A + 2*sqrt(ysq), A - 2*sqrt(ysq)];