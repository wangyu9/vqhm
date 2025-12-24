function [V2,F2,eq_lhs2,eq_rhs2] = subdivide_with_constraint(V,F,eq_lhs,eq_rhs,num_iters)



V2 = V;
F2 = F;
eq_lhs2 = eq_lhs;
eq_rhs2 = eq_rhs;

for outloop=1:num_iters 

%
V_old = V2;
F_old = F2;
eq_lhs_old = eq_lhs2;
eq_rhs_old = eq_rhs2;


[VVB,FF] = upsample([V_old speye(size(V_old,1))],F_old); % 'OnlySelected',invTri
% [VVB,FF] = upsample([V_old speye(size(V_old,1))],F_old);

VV = full(VVB(:,1:size(V_old,2)));
M = VVB(:,size(V_old,2)+1:end);
max(max(abs(VV - M*V_old)));

MM = blkdiag(M,M);

fc = find(sum(MM * (eq_lhs_old)',2)==1);

eq_lhs_new = sparse(1:numel(fc), fc, 1, numel(fc), 2*size(VV,1));
eq_rhs_new = eq_lhs_new * MM * eq_lhs_old' * eq_rhs_old;

%
V2 = VV;
F2 = FF;
eq_lhs2 = eq_lhs_new;
eq_rhs2 = eq_rhs_new;


%
end