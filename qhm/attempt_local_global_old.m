%% A local global solver: Initialization.

au = tp.s_at2au(reshape(da,[f,3]));
au = reshape(au,[f,3]);
a11 = au(:,1);
a12 = au(:,2);
a22 = au(:,3);

GI = mesh.GI;
% Area = mesh.FA;
%%
for jj=1:20
% global step:

disp('attempting local-global step\n');

Lw = GI'*[sparse_diag(a11),sparse_diag(a12);sparse_diag(a12),sparse_diag(a22)]*GI;
lhs = S'*Lw*S;
rhs = - S'*Lw*R * TVB;

U = lhs \ rhs;

W = S * U + R * TVB;

u = W(:,1);
v = W(:,2);

trace( W' * Lw * W / 2 );

% render meshes.

% render_mesh2(W,F,'EdgeColor',[0,0,0]);
% axis equal;

% local step:

if false
    % the naive slow implementation:
    for j=1:f
        Jj = (GI([j,j+f],:) * W)';
        Aj = abs(det(Jj)) * inv(Jj' * Jj) * Area(j);
        a11(j) = Aj(1,1);
        a12(j) = Aj(1,2);
        a22(j) = Aj(2,2);
    end
else
    % the fast implementation. 
GW = GI * W;
Gxu = GW(1:f,1);
Gxv = GW(1:f,2);
Gyu = GW(f+1:2*f,1);
Gyv = GW(f+1:2*f,2);
% Jj = [Gxu, Gxv; Gyu, Gyv]';

adetJ = abs( Gxu.*Gyv - Gxv.*Gyu );

a22 = (Gxu .* Gxu + Gxv .* Gxv) ./ adetJ .* Area;
a12 = - (Gxu .* Gyu + Gxv .* Gyv) ./ adetJ .* Area; 
a11 = (Gyu .* Gyu + Gyv .* Gyv) ./ adetJ .* Area;
end

%
newArea = doublearea([u,v], F);
flipped = newArea < 0;
% number of flipped triangles. 
num_flipped = numel(find(flipped));
fprintf('flipps=%d, min_area=%g',num_flipped,min(newArea));

if num_flipped==0
    break;
end
%%
end

