function [e, g, FW] = val_grad_qhm(xx,pre,s_at2au,s_value_grad,fixed)

B = pre.known;

assert(pre.dim==2);
mcdim = 3;

aatt = xx(1:pre.f*mcdim);
BC = xx(1+pre.f*mcdim:end);

aatt = reshape(aatt,[pre.f,numel(aatt)/pre.f]);
BC = reshape(BC, [size(B,1), numel(BC)/size(B,1)]);

[e, g_at, g_BC, FW] = QHWvalue(s_at2au(aatt),pre,aatt, BC, s_value_grad);

% g_BC(1:size(fixed,1),:) = 0;

g = [g_at;g_BC(:)];

