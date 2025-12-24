function [VT_out] = normalize_mesh(V,F,VT)

area = sum(0.5*doublearea(V,F));

mesh = triangle_mesh(V,F);

n = size(V,1);

known = mesh.B;
unknown = mesh.UB;
row_slice_matrix = @(indices) sparse(indices,1:numel(indices),1,n,numel(indices));

R = row_slice_matrix(known);
S = row_slice_matrix(unknown);

BC = VT(mesh.B,:);

NVT = mesh.D*(R* [BC(:,2),-BC(:,1)] ); 

BN = R' * NVT;

area_t = 0.5 * ( BC(:,1)' * BN(:,1) + BC(:,2)' * BN(:,2)  );

assert(area_t>0); 

VT_out = sqrt(area./area_t) * VT;
