function [V,F] = rectangle_mesh_simple(x0,y0,dx,dy,m,n)
% (x0,y0) is the left down point of rect domain.
% (dx,dy) is the size of each grid.
% (m,n) is the number of grids along two direction.
% note the normal should point to z direction if using the rand hand rule.
% by wangyu

V = [];
for j=1:1:n+1
   x = x0+(0:m)*dx;  
   y = (y0+(j-1)*dy)*ones(1,m+1);
   V = [V;x',y'];
end

% F = [];
% assert(m>=2);
% if(j>0)
%     F = 1
% end



F = delaunay(V(:,1),V(:,2));