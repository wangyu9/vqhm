function [J] = linear_tensor_to_edge_weights(V,F)

g1 = core(V,F(:,[1,2,3]));
g2 = core(V,F(:,[2,3,1]));
g3 = core(V,F(:,[3,1,2]));

g1 = g1(:,1:2);
g2 = g2(:,1:2);
g3 = g3(:,1:2);


f = size(F,1);

J = zeros(f,3,3);
    
for i=1:f
    J(i,:,:) = [...
        g1(i,1)*g2(i,1), g1(i,1)*g2(i,2)+g1(i,2)*g2(i,1), g1(i,2)*g2(i,2);...
        g2(i,1)*g3(i,1), g2(i,1)*g3(i,2)+g2(i,2)*g3(i,1), g2(i,2)*g3(i,2);...
        g3(i,1)*g1(i,1), g3(i,1)*g1(i,2)+g3(i,2)*g1(i,1), g3(i,2)*g1(i,2);...
    ];
end

%% optionally double check 
if false
%%
for i=1:f
   assert(abs(det(squeeze(A2w(i,:,:))))>1e-10); 
end
end

end
    
function [r] = core(V,F)

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

