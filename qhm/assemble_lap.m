function [A_lower, A_full] = assemble_lap(n,GIS,F,a11,a12,a22)

A_lower = [];
A_full = [];

%

g1 =  GIS.g1;
g2 =  GIS.g2;
g3 =  GIS.g3;


output_full = true;
output_lower = true; % true;

if output_full

    II = [F(:,1);F(:,2);F(:,3);F(:,1);F(:,2);F(:,3);F(:,2);F(:,3);F(:,1);];
    JJ = [F(:,1);F(:,2);F(:,3);F(:,2);F(:,3);F(:,1);F(:,1);F(:,2);F(:,3);];

    % incorrect as it will merge directed edges: [RI, RJ, RO] = find( sparse(II,JJ,1:numel(II),n,n)  );

    % sorted first by columns and then by rows. 
    [BB,RO] = sortrows([JJ,II],'ascend'); % so that BB = [...](RO,:)
    RI = BB(:,2);
    RJ = BB(:,1);

end

%II = [F(:,1);F(:,2);F(:,3);F(:,2);F(:,3);F(:,1);];
%JJ = [F(:,2);F(:,3);F(:,1);F(:,1);F(:,2);F(:,3);];

%II = [F(:,1);F(:,2);F(:,3);];
%JJ = [F(:,2);F(:,3);F(:,1);];

%II2 = [F(:,1);F(:,2);F(:,3);F(:,1);F(:,2);F(:,3);F(:,2);F(:,3);F(:,1);];
%JJ2 = [F(:,1);F(:,2);F(:,3);F(:,2);F(:,3);F(:,1);F(:,1);F(:,2);F(:,3);];

if output_lower

    maxFF12 = max(F(:,1),F(:,2)); 
    minFF12 = min(F(:,1),F(:,2)); 
    
    maxFF23 = max(F(:,2),F(:,3)); 
    minFF23 = min(F(:,2),F(:,3));
    
    maxFF31 = max(F(:,3),F(:,1)); 
    minFF31 = min(F(:,3),F(:,1)); 
    
    % [F(:,1);F(:,2);F(:,3);F(:,1);F(:,2);F(:,3);];
    % [F(:,1);F(:,2);F(:,3);F(:,2);F(:,3);F(:,1);];
    
    II_lower = [F(:,1);F(:,2);F(:,3);maxFF12;maxFF23;maxFF31;];
    JJ_lower = [F(:,1);F(:,2);F(:,3);minFF12;minFF23;minFF31;];

    
    % sorted first by columns and then by rows. 
    [BB_lower,RO_lower] = sortrows([JJ_lower,II_lower],'ascend');
    RI_lower = BB_lower(:,2);
    RJ_lower = BB_lower(:,1);

end

inner_prod = @(ss,tt) ss(:,1) .* a11 .* tt(:,1) +  a12 .* ( ss(:,1) .* tt(:,2) +  ss(:,2) .* tt(:,1)) + ss(:,2) .* a22 .* tt(:,2); 

% v12=v21, v23=v32, v31=v13

v12 = inner_prod(g1, g2); 
v23 = inner_prod(g2, g3); 
v31 = inner_prod(g3, g1); 
% no need to multiple with mesh.AI again. already did so. 

v11 = - (v12 + v31); 
v22 = - (v23 + v12);
v33 = - (v31 + v23); 

if output_full

    VV = [v11;v22;v33;v12;v23;v31;v12;v23;v31];
    
    %VV = [v12;v23;v31];
    
    % old slow code: A_full = sparse2(II,JJ,VV,n,n,numel(VV)); 

    A_full = sparse2(RI,RJ,VV(RO),n,n,numel(VV)); 
end

if output_lower

    VV_lower = [v11;v22;v33;v12;v23;v31;];
    
    % old slow code: A_lower = sparse2(II_lower,JJ_lower,VV_lower,n,n,numel(VV_lower));

    A_lower = sparse2(RI_lower,RJ_lower,VV_lower(RO_lower),n,n,numel(VV_lower)); 

end

% sum(sum(abs(A2-A)))

% https://www.mathworks.com/matlabcentral/answers/509170-update-a-sparse-matrix-efficiently?s_tid=srchtitle
% fastest way to construct a sparse matrix:
% inputs are sorted, first by columns and then by rows. 

% [siA, sjA, svA] = find(A_full);
% tic; sA = sparse(siA, sjA, svA, n, n); toc