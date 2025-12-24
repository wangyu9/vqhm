function [BE, cind, circular_ops] = circular_laplacian(n,F)

% BE: each row is the two indices for each edge. 
% BE are indieces into V, i.e. all vertices. 
% BE(:,1) is the list of boundary vertices. 

% cind is the list of boundary edges that is in *counterclockwise* order. 

% cind is index into rows of BE. 
% cind is also index into B=BE(:,1), the list of boundary vertices. 

% BE(cind,:) is the list of boundary vertices in the same order as boundary
% edges. 

% ---

% BE = outline(F);
% copied from outline: 

% build directed adjacency matrix
A = sparse(F,F(:,[2:end 1]),1);
% Find single occurance edges
[OI,OJ,OV] = find(A-A');
% Maintain direction
BE = [OI(OV>0) OJ(OV>0)];%;OJ(OV<0) OI(OV<0)];

B = BE(:,1);
[~,~,indE1] = intersect(BE(:,2),BE(:,1),'stable');
[~,~,indE2] = intersect(BE(:,1),BE(:,2),'stable');
assert(norm(BE(indE1,1)-BE(:,2))==0);
assert(norm(BE(indE2,2)-BE(:,1))==0);
%
cind = [1];
for i=1:size(indE1,1)-1
   new_ind = indE1(cind(end));
   cind = [cind;new_ind]; 
end

% # boundary edges 
nbe = size(BE,1);
% # boundary vertices
nbv = nbe; % asserting simply connected domain. 

assert(numel(cind)==nbe); % it only works for simply connected meshes. 
% cind is index into rows of BE.
%
scind = [cind,[cind(2:end);cind(1)]];

% for each boundary edge BE(i,:), BEindF(i) is the index of triangle that
% BE belongs to. 

f = size(F,1);
[CC,~,IB] = intersect(BE,[F(:,[1,2]); F(:,[2,3]); F(:,[3,1]) ],'stable','rows');
IF = [(1:f)';(1:f)';(1:f)'];
assert(norm(CC-BE)==0);

BEindF = IF(IB);

circular_ops = struct();

circular_ops.BEindF = BEindF;

L_be = 2* speye(numel(cind)) - sparse(scind,scind(:,[2,1]), 1 );
% ---

succ = [cind(2:end);cind(1)];
prev = [cind(end);cind(1:end-1)];
% like cind, succ, prev are indice into rows of BE and B. 

circular_ops.L_be = L_be;
circular_ops.succ = succ;
circular_ops.prev = prev;

Dp = sparse(prev,1:nbe,1,nbv,nbe);
Ds = sparse(succ,1:nbe,1,nbv,nbe);

circular_ops.Dp = Dp;
circular_ops.Ds = Ds;

% 
v2bv = sparse(1:nbv,BE(:,1),1,nbv,n);
GB = (Ds - Dp)' * v2bv;                                                                                                                                                                                      
% GB calculates gradients along the boundary edges. 

circular_ops.GB = GB;

circular_ops.Dirac = (sparse(B(cind), B(prev), 0.5, n, n) - sparse(B(cind), B(succ), 0.5, n, n) );
%% debugging space 

if false
%% visualization
figure;
hold on;
BBB = BE(:,1);
for i=1:numel(cind)
    vv = V(BBB(cind(i)),:);
    vv = V(BBB(succ(i)),:);
    vv = V(BBB(prev(i)),:);
    vx = vv(1);
    vy = vv(2);
    scatter(vx,vy,10);
    pause(0.1);
end
%
%%
VB = V(B,:);
dV = V(BE(:,2),:) - V(BE(:,1),:);
quiver(VB(:,1),VB(:,2),dV(:,1),dV(:,2));
end

