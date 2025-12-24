function [V,F] = symmetric_rectangle_mesh2(x0,y0,dx,dy,hm,hn)
% (x0,y0) is the left down point of rect domain.
% (dx,dy) is the size of each grid.
% (hm,hn) is half of the number of grids along two direction.

% by wangyu

m = 2*hm;
n = 2*hn;

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



%F = delaunay(V(:,1),V(:,2)); this does not give me symmetric mesh.

F = [];
for j=1:1:hn
    
   leftdown = 1+(j-1)*(m+1);
   num_x = hm;
   index_delta_x = 1;
   index_delta_y = m+1;   
   
   ld = (leftdown:(leftdown+num_x-1))';
   rd = ld + index_delta_x;
   lu = ld + index_delta_y;
   ru = rd + index_delta_y;
   % \ style
   dF = [ld,rd,ru;ld,ru,lu];    
   F = [F;dF];
   
   
   %
   leftdown = 1+(j-1)*(m+1)+hm;
   num_x = hm;
   index_delta_x = 1;
   index_delta_y = m+1;

   ld = (leftdown:(leftdown+num_x-1))';
   rd = ld + index_delta_x;
   lu = ld + index_delta_y;
   ru = rd + index_delta_y;
   % / style
   dF = [ld,rd,lu;rd,ru,lu];
   F = [F;dF];
end

for j=hn+1:1:2*hn
    
   leftdown = 1+(j-1)*(m+1);
   num_x = hm;
   index_delta_x = 1;
   index_delta_y = m+1;   
   
   ld = (leftdown:(leftdown+num_x-1))';
   rd = ld + index_delta_x;
   lu = ld + index_delta_y;
   ru = rd + index_delta_y;
   % / style
   dF = [ld,rd,lu;rd,ru,lu];
   F = [F;dF];
   
   
   %
   leftdown = 1+(j-1)*(m+1)+hm;
   num_x = hm;
   index_delta_x = 1;
   index_delta_y = m+1;

   ld = (leftdown:(leftdown+num_x-1))';
   rd = ld + index_delta_x;
   lu = ld + index_delta_y;
   ru = rd + index_delta_y;
   % \ style
   dF = [ld,rd,ru;ld,ru,lu];    
   F = [F;dF];
end