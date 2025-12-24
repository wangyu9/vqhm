function [G,Area,out] = intrinsic_grad_area(V,F,epsilon)


absorb_mass = true;


out = struct();
out.absorb_mass = absorb_mass;


n = size(V,1);
f = size(F,1);

E = edge_lengths(V,F);



if false
    
E = E ./ max(E,[],2);

% ll = find(E<epsilon); 
% E(ll) = E(ll) + epsilon;

ep1 = 0.6;
ep2 = 0.3;

tt = find(min(E,[],2)<ep1); 
% E(tt,:) = E(tt,:) + epsilon;
Et = E(tt,:);

Et = abs(ep1-min(Et,[],2)) + Et;

% tl = find(Et<epsilon);
% Et(tl) =  Et(tl) + epsilon;
E(tt,:) = Et;

Area = doublearea_intrinsic(E)/2; 

end


% F(:,1) are always put at the origin.
% F(:,2) are put at (e3, 0)
% F(:,3) are stored at (e2 * cos(theta1), h3)
% cos(theta1):
% cos1 = (E(:,2).^2 + E(:,3).^2 - E(:,1).^2) ./ (2*E(:,2).*E(:,3));


% V1 = zeros(f,2);
% V2 = [ E(:,3), zeros(f,1)];
% V3 = [ (E(:,2).^2 + E(:,3).^2 - E(:,1).^2) ./ (2*E(:,3)), 2*Area ./ E(:,3)];



AAA = internalangles_intrinsic(E) / pi; 

if true
% the old version
    AAA(find(AAA<epsilon)) = epsilon;

else
    upper = max(0, 1-2*epsilon); 
    AAA(find(AAA>upper)) = upper;

    AAA(find(AAA<epsilon)) = epsilon;
end

AAA = pi * AAA ./ sum(AAA,2);

    
V1 = zeros(f,2);
V2 = [ cot(AAA(:,1)) + cot(AAA(:,2)) , zeros(f,1)];
V3 = [ cot(AAA(:,1)), ones(f,1)];

Area = ( cot(AAA(:,1)) + cot(AAA(:,2)) )/2; 

assert(min(min(Area))>0);

if absorb_mass
    V1 = V1 ./ sqrt(Area); 
    V2 = V2 ./ sqrt(Area); 
    V3 = V3 ./ sqrt(Area); 
    Area = ones(size(Area)); 
end



g1 = core2d(V1,V2,V3);
g2 = core2d(V2,V3,V1);
g3 = core2d(V3,V1,V2);

g1 = g1(:,1:2);
g2 = g2(:,1:2);
g3 = g3(:,1:2);


Gx =  sparse(1:f,F(:,1),g1(:,1),f,n) ...
    + sparse(1:f,F(:,2),g2(:,1),f,n) ...
    + sparse(1:f,F(:,3),g3(:,1),f,n);

Gy =  sparse(1:f,F(:,1),g1(:,2),f,n) ...
    + sparse(1:f,F(:,2),g2(:,2),f,n) ...
    + sparse(1:f,F(:,3),g3(:,2),f,n);

G = [Gx;Gy];

out.g1 = g1;
out.g2 = g2;
out.g3 = g3;

end

function [r] = core2d(V1, V2, V3)

% v10, v20

% this computes gradient of the dirac function centered at vertex 0
% in this triangle. The gradient is extrinsic in the 3d space.


v20 = V3 - V1;
v10 = V2 - V1;
        
dot12 = sum(v20.*v10, 2);

% aa(1)*bb(2) - aa(2)*bb(1);
% cross(v20, v10)

c21 = v20(:,1).*v10(:,2) - v20(:,2).*v10(:,1);

ns12 = c21.^2;

r = bsxfun(@times,v10, sum(v20.*v20, 2)./ns12) ...
  + bsxfun(@times,v20, sum(v10.*v10, 2)./ns12) ...open
  - bsxfun(@times,v10 + v20, dot12./ns12);

end