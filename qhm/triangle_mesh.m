function [mesh] = triangle_mesh(V,F,varargin)
%%
mesh = struct();

n = size(V,1);
f = size(F,1);

mesh.n = n;
mesh.f = f;

dim = size(F,2) - 1;
assert(dim==2||dim==3);
mesh.dim = dim;

mesh.V = V;
mesh.F = F;



mesh.D = mesh_dirac(n,F);

%
mesh.FA = doublearea(V,F)/2;

EL = edge_lengths(V,F);
E = edges(F);
ne = size(E,1);

mesh.EL = EL;
mesh.E = E;
mesh.ne = ne;

assert(norm(sort(E')'-E)==0);

% BE = outline(F);
[BE, cind, circular_ops] = circular_laplacian(n,F);
mesh.cycle = circular_ops;
mesh.cind = cind;

B = BE(:,1);
SBE = sort(BE,2); % sort will sort columns, this does the trick for rows. 
nbe = size(BE,1);
nb = nbe; % asserting simply connected domain. 

UB = find(~sparse(1,B,true,1,n));
UB = UB(:);

mesh.B = B;
mesh.UB = UB;
mesh.BE = BE;
[~,~,be2e] = intersect(SBE,E,'stable','rows'); 
assert(norm(SBE-E(be2e,:))==0);

mesh.be2e = be2e;

try

% BT: list of boundary triangles:
ind_B = zeros(n,1);
ind_B(B) = 1;
BT = find(ind_B(F(:,1))+ind_B(F(:,2))+ind_B(F(:,3))>=2); % note triangle with only one boundary vertex has no boundary edge.

assert(max(ind_B(F(:,1))+ind_B(F(:,2))+ind_B(F(:,3)))<3); % assuming each boundary triangle has at most one boundary edge. 

% BTibe: such that for the k-th triangle in BT, the boundary edge is the
% BTibe(k)-th edge in that triangle. 
[~,BTibe] = min( [ind_B(F(BT,1)),ind_B(F(BT,2)),ind_B(F(BT,3))], [], 2);
BTvbe = zeros(numel(BT),2);
BTcbe = zeros(numel(BT),2);
for i=1:numel(BT)
    BTvbe(i,:) = ...
        - V( F(BT(i), mod(BTibe(i)-1-1,3)+1  ), 1:2 ) ...
        + V( F(BT(i), mod(BTibe(i)+1-1,3)+1  ), 1:2 );
    BTcbe(i,:) = ...
         ( V( F(BT(i), mod(BTibe(i)-1-1,3)+1  ), 1:2 ) ...
         + V( F(BT(i), mod(BTibe(i)+1-1,3)+1  ), 1:2 ) ) / 2;
end

mesh.BT = BT;
mesh.BTibe = BTibe;
mesh.BTvbe = BTvbe;
mesh.BTcbe = BTcbe;

catch
    
end

if false
    %%
    figure;
    hold on;
    BBB = BE(:,1);
    qline = @(p1,p2) quiver(p1(1),p1(2),p2(1) - p1(1), p2(2) - p1(2),0);

    for i=1:numel(BT)
        TTT = zeros(2,2);
        TTT(1,:) = V( F(BT(i), mod(BTibe(i)-1-1,3)+1  ), 1:2 );
        TTT(2,:) = V( F(BT(i), mod(BTibe(i)+1-1,3)+1  ), 1:2 );
        % line(TTT(:,1),TTT(:,2));
        qline(TTT(1,:),TTT(2,:));
        pause(0.1);
    end
end

if false
    %%
    figure;
    hold on;
    BBB = BE(:,1);
    qline = @(p1,d) quiver(p1(1),p1(2),d(1), d(2),0);

    for i=1:numel(BT)
        TTT = zeros(2,2);
        TTT(1,:) = V( F(BT(i), mod(BTibe(i)-1-1,3)+1  ), 1:2 );
        TTT(2,:) = V( F(BT(i), mod(BTibe(i)+1-1,3)+1  ), 1:2 );
        % line(TTT(:,1),TTT(:,2));
        qline( TTT(1,:), BTvbe(i,:));
        pause(0.1);
    end
end


%% FEM related

g1 = basis_grad_core(V,F(:,[1,2,3]));
g2 = basis_grad_core(V,F(:,[2,3,1]));
g3 = basis_grad_core(V,F(:,[3,1,2]));

g1 = g1(:,1:2);
g2 = g2(:,1:2);
g3 = g3(:,1:2);

mesh.g1 = g1;
mesh.g2 = g2;
mesh.g3 = g3;

%% Boundary Curve related

JJJ = sparse(be2e,1:nbe,1,ne,nbe);
GB = JJJ * ...
    ( sparse(1:nbe,BE(:,2),1,nbe,n) - sparse(1:nbe,BE(:,1),1,nbe,n) );


% VB = V(B,:);
% dV = V(BE(:,2),:) - V(BE(:,1),:); % these are counter-clockwise. 
% quiver(VB(:,1),VB(:,2),dV(:,1),dV(:,2)); 

% IBE is similar to BE with the same size, but the indices are into B s.t. 
% B(IBE(:))==BE(:)

[~,~,IBE1] = intersect(BE(:,1),B,'stable'); 
[~,~,IBE2] = intersect(BE(:,2),B,'stable'); 
IBE = [IBE1,IBE2];

assert(norm(BE - [B(IBE1),B(IBE2)])==0);

Dp = sparse(IBE(:,1),1:nbe,1,nb,nbe) * JJJ';
Ds = sparse(IBE(:,2),1:nbe,1,nb,nbe) * JJJ';

v12 = V(F(:,2),:) - V(F(:,1),:); 
v13 = V(F(:,3),:) - V(F(:,1),:); 

v23 = V(F(:,3),:) - V(F(:,2),:); 
v21 = V(F(:,1),:) - V(F(:,2),:); 

v31 = V(F(:,1),:) - V(F(:,3),:); 
v32 = V(F(:,2),:) - V(F(:,3),:); 
        
dot23 = sum(v12.*v13, 2);
dot31 = sum(v23.*v21, 2);
dot12 = sum(v31.*v32, 2);

mesh.dot23 = dot23;
mesh.dot31 = dot31;
mesh.dot12 = dot12;


% 2D cross
cross23 = v12(:,1) .* v13(:,2) - v12(:,2) .* v13(:,1); 
cross31 = v23(:,1) .* v21(:,2) - v23(:,2) .* v21(:,1); 
cross12 = v31(:,1) .* v32(:,2) - v31(:,2) .* v32(:,1); 

mesh.cross23 = cross23;
mesh.cross31 = cross31;
mesh.cross12 = cross12;

% https://en.wikipedia.org/wiki/Tangent_half-angle_formula
% the *signed* tan angle 213 (right hand rule for triangle 123 for positive sign)
% tan(theta/2) = sin(theta) / (1+cos(theta))
stan1 = cross23 ./ (dot23 + EL(:,2).*EL(:,3) ); 
stan2 = cross31 ./ (dot31 + EL(:,3).*EL(:,1) );
stan3 = cross12 ./ (dot12 + EL(:,1).*EL(:,2) );

%assert(min(stan1)>=0); 
%assert(min(stan2)>=0);
%assert(min(stan3)>=0); 

% mesh.stan1 = stan1;
% mesh.stan2 = stan2;
% mesh.stan3 = stan3;

stan = [stan1,stan2,stan3];
mesh.stan = stan;


mV = [stan./EL(:,[3,1,2]), stan./EL(:,[2,3,1])];
mI = [F,F];
mJ = [F(:,[2,3,1]),F(:,[3,1,2])];

% MVL = sparse(II,JJ,VV,n,n);
% MVL = MVL - diag(sum(MVL,2));

mesh.WDE = mV;
mesh.mI = mI;
mesh.mJ = mJ;

mesh.assemble_MVL = @(II,JJ,VV) sparse(II,JJ,VV,n,n) - sparse(II,II,VV,n,n);

% equivalently:
MVL = sparse(mI,mJ,mV,n,n) - sparse(mI,mI,mV,n,n);

MVL = mesh.assemble_MVL(mI,mJ,mV);
mesh.MVL = MVL;


% sum(sum(abs( MVL-mean_value_laplacian(V,F) ))) 
% gives a number around 0.  


mesh.wde2vec = @(WDE) WDE(:); % WDE is a fx3 vector storing data for each directed edge. 
mesh.vec2wde = @(aa) reshape(aa, [f,3]);


X = V;
if size(X,2)==2
   X = [X,zeros(size(X,1),1)]; 
end

try
    % some experiments that is not used in the paper. 
    meshdec = getMeshData(X,F,2,'none');
    DEC = discreteExteriorCalculus(meshdec);
    assert(norm(E-meshdec.Elist)==0);
    
    mesh.DEC = DEC;
catch

end


G = grad(V,F);
% assert(size(G,1)==dim*f);
mesh.G = G;

[mesh.GI, mesh.GIS] = intrinsic_grad(V,F);
end


function [r] = basis_grad_core(V,F)

% v10, v20

% this computes gradient of the dirac function centered at vertex 0
% in this triangle. The gradient is extrinsic in the 3d space.


if size(V,2)==2
   V = [V,zeros(size(V,1),1)]; 
end

v20 = V(F(:,3),:) - V(F(:,1),:);
v10 = V(F(:,2),:) - V(F(:,1),:); 
        
dot12 = sum(v20.*v10, 2);

ns12 = sum(cross(v20, v10).^2, 2);


r = bsxfun(@times,v10, sum(v20.*v20, 2)./ns12) ...
  + bsxfun(@times,v20, sum(v10.*v10, 2)./ns12) ...open
  - bsxfun(@times,v10 + v20, dot12./ns12);

end

function [g0] = basis_grad_core3d(V,T)

% this computes gradient of the dirac function centered at vertex 0
% in this tet. The gradient is extrinsic in the 3d space.


assert(size(V,2)==3)

v30 = V(T(:,4),:) - V(T(:,1),:);
v20 = V(T(:,3),:) - V(T(:,1),:);
v10 = V(T(:,2),:) - V(T(:,1),:); 

v31 = v30 - v10;
v21 = v20 - v10;

% 1/6 det is the volume of tet. 
% 1/6 det = 1/3 area(123) * h
% area(123) = 1/2 cross(v31, v21)


% h = det([v10;v20;v30]) ./ cross(v31, v21);
% the code above does not work, need det per row

% n is the vector pendicular to the plane 
% of tri(123)  

n123 = cross(v21, v31);
n123 = bsxfun(@times, n123, 1./sqrt(sum(n123.^2,2)));

g0 = bsxfun(@times, n123, -1./sum(n123.*v10, 2) );

end