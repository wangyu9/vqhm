function [e, g, FW] = val_grad_qhm_edge(xx,pre,s_at2au,s_value_grad,option)

B = pre.known;

assert(pre.dim==2);

ne = size(pre.star,1);

if option.at_only
aatt = xx(:);
BC = option.BC;
else
aatt = xx(1:ne);
BC = xx(1+ne:end);
end

BC = reshape(BC, [size(B,1), numel(BC)/size(B,1)]);

[e, g_at, g_BC, FW] = FEMvalue_ew(s_at2au(aatt),pre,aatt, BC, s_value_grad, option.CR);

if ~isempty(option.PBC)
   g_BC = option.PBC(g_BC); 
end

% g_BC(1:size(option.fixed,1),:) = 0;
g_BC = g_BC(:);

if ~isempty(option.PA)
   g_at = option.PA(g_at); 
end

if option.at_only
g = g_at;
else    
g = [g_at;g_BC];
end

