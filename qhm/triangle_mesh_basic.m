function [mesh] = triangle_mesh_basic(V,F,indEC,varargin)
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

if false
    EL = edge_lengths(V,F);
    E = edges(F);
    ne = size(E,1);
    
    mesh.EL = EL;
    mesh.E = E;
    mesh.ne = ne;
    
    assert(norm(sort(E')'-E)==0);
end

BE = outline(F);

B = BE(:,1);

UB = find(~sparse(1,B,true,1,n));
UB = UB(:);

% mesh.B = B;
% mesh.UB = UB;

mesh.BE = BE;


mesh.IKB = [B;indEC]; 

mesh.IUB = find(~sparse(1,mesh.IKB,true,1,n));
mesh.IUB= mesh.IUB(:); 

assert(numel(mesh.IKB)==numel(unique(mesh.IKB)));


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
