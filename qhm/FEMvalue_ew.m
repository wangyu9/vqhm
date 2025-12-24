function [e,g,g_BC,W,extra] = FEMvalue_ew(au,pre,at, BC,s_val_grad, CR)

%%
verbose = 1;
extra = struct();

%%

f = pre.f;
dim = pre.dim;
n = pre.n;
solver = 'cholmod';
m = size(BC,2);
W = zeros(n,m);
W(pre.known,:) = BC(:,:);
pre.iter = pre.iter + 1;
fprintf('Iter: %d',pre.iter);
nu = length(pre.unknown);
%%
V = pre.V;
F = pre.F;

%
X = V;
if size(X,2)==2
   X = [X,zeros(size(X,1),1)]; 
end
mesh = getMeshData(X,F,2,'none');
DEC = discreteExteriorCalculus(mesh);
Lfz = DEC.d01;
%%
ne = size(Lfz,1);
%%

%A = symmetric_tensor_assemble(au,dim); 
A = spdiags(au,0,ne,ne);

Lw = Lfz' * A * Lfz;

extra.Lw = Lw;

% pre.old_Wu = - linear_solver( Lw(pre.unknown,pre.unknown), (Lw(pre.unknown,pre.known)*BC(:,:)), solver);
pre.old_Wu = -  Lw(pre.unknown,pre.unknown) \ (Lw(pre.unknown,pre.known)*BC(:,:));


extra.LA = Lw(pre.unknown,pre.unknown);
extra.LB = - Lw(pre.unknown,pre.known);

W(pre.unknown,:) = pre.old_Wu;




%%

[e, dW] = s_val_grad(W);

dW = reshape(dW, size(W));

% dW = dW + CR * 1/2 * Lw * W;


% e = trace((W'*pre.L)*(pre.invMass*(pre.L*W)));
% fprintf('Bilaplacian energy: %g, quasi Dirichlet energy: %g, with min=%g.\n',e,trace(W'*Lw*W),min(min(W)));
% Res = 2 * pre.L*(pre.invMass*(pre.L*W));

if verbose>=1
    fprintf('energy=%f\n',e);
end

%%

% pre.old_BR = linear_solver(Lw(pre.unknown,pre.unknown),dW(pre.unknown,:),solver);
pre.old_BR = Lw(pre.unknown,pre.unknown) \ dW(pre.unknown,:);

BR = pre.old_BR;
PS = Lfz(:,pre.unknown) * BR;

%
% GW = pre.G(:,pre.unknown)*W(pre.unknown,:) + pre.G(:,pre.known)*W(pre.known,:);
GW = Lfz * W;
for j=1:m
    if j==1
        dEda = zeros(ne,1);
    end

    if true % this is much slower. 
        spXj = spdiags(GW(:,j),0,ne,ne);
        dEda = dEda - spXj' * PS(:,j) + CR * 1/2 * spXj' * GW(:,j);
        if j==1
        extra.Tx = spXj' * Lfz(:,pre.unknown);
        end
        if j==2
        extra.Ty = spXj' * Lfz(:,pre.unknown);
        end
        if j==3
        extra.Tz = spXj' * Lfz(:,pre.unknown);
        end
    else
        dEda = dEda - symmetric_tensor_span_dot(GW(:,j),PS(:,j),dim);
    end
end

g = dEda;

extra.g_au = g;

% % g = pre.s_pdapdt(at)' * g; % too slow
g = pre.s_pdapdt_lmul(at,g);
% g = g; % no reparameterization. 

g_BC = dW(pre.known,:) - Lw(pre.unknown,pre.known)' * (Lw(pre.unknown,pre.unknown)\dW(pre.unknown,:));


%H = speye(size(g,1));
H = [];

%if ~isempty(pre.handler)
%    pre.handler(pre, at, g);
%end

if false
    folder_path = 'examples/rect-200by40/W_qhe_lbfgs/';
    save(sprintf([folder_path,'%5d.mat'],pre.iter), 'at', 'g', 'au'); 
end


%%
if false
%
render_mesh2(W,pre.F,'EdgeColor',[0,0,0]);
hold on; 
quiver(W(:,1),W(:,2),-dW(:,1),-dW(:,2));
%
render_mesh2(W,pre.F,'EdgeColor',[0,0,0]);
hold on;
quiver(W(pre.known,1),W(pre.known,2),-g_BC,-g_BC);
end