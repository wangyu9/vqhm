function [out] = assemble_lap_core(n,F)

out = struct();

out.asb_lower = [];
out.asb_full = [];



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


if output_full

    
    
    %VV = [v12;v23;v31];
    
    % old slow code: A_full = sparse2(II,JJ,VV,n,n,numel(VV)); 

    % out.asb_full = @(GIS,a11,a12,a22) sparse(RI,RJ,get_value(GIS,a11,a12,a22,RO,true),n,n,numel(RI)); 
    out.asb_full = @(GIS,a11,a12,a22) sparse2(RI,RJ,get_value(GIS,a11,a12,a22,RO,true),n,n,numel(RI)); 

end

if output_lower

    
    
    % old slow code: A_lower = sparse2(II_lower,JJ_lower,VV_lower,n,n,numel(VV_lower));

    % out.asb_lower = @(GIS,a11,a12,a22) sparse(RI_lower,RJ_lower, get_value(GIS,a11,a12,a22,RO_lower, false),n,n,numel(RI_lower)); 
    out.asb_lower = @(GIS,a11,a12,a22) sparse2(RI_lower,RJ_lower, get_value(GIS,a11,a12,a22,RO_lower, false),n,n,numel(RI_lower)); 

end

% sum(sum(abs(A2-A)))

% https://www.mathworks.com/matlabcentral/answers/509170-update-a-sparse-matrix-efficiently?s_tid=srchtitle
% fastest way to construct a sparse matrix:
% inputs are sorted, first by columns and then by rows. 

% [siA, sjA, svA] = find(A_full);
% tic; sA = sparse(siA, sjA, svA, n, n); toc


end

function [V_out] = get_value(GIS, a11,a12,a22,perm, is_full)


    %
    
    g1 =  GIS.g1;
    g2 =  GIS.g2;
    g3 =  GIS.g3;


    inner_prod = @(ss,tt) ss(:,1) .* a11 .* tt(:,1) +  a12 .* ( ss(:,1) .* tt(:,2) +  ss(:,2) .* tt(:,1)) + ss(:,2) .* a22 .* tt(:,2); 
    
    % v12=v21, v23=v32, v31=v13
    
%    epsilon = 0; %1e-8; 

    v12 = inner_prod(g1, g2); 
    v23 = inner_prod(g2, g3); 
    v31 = inner_prod(g3, g1); 

%     v12 = v12 + epsilon;
%     v23 = v23 + epsilon;
%     v31 = v31 + epsilon;

    % no need to multiple with mesh.AI again. already did so. 
    
    v11 = - (v12 + v31); 
    v22 = - (v23 + v12);
    v33 = - (v31 + v23); 

    if is_full

        VV = [v11;v22;v33;v12;v23;v31;v12;v23;v31];

        V_out = VV(perm); 

    else

        VV_lower = [v11;v22;v33;v12;v23;v31;];

        V_out = VV_lower(perm); 

    end


end