%% A local global solver: Initialization.

old_imp = false;

if old_imp 

au = tp.s_at2au(reshape(da,[f,3]));
au = reshape(au,[f,3]);
GI = mesh.GI;

else

[au, ~] = tensor_para_faster(reshape(da,[f,3]),tp.para_type);
au = au .* Area;    

end

a11 = au(:,1);
a12 = au(:,2);
a22 = au(:,3);

%%
for jj=1:20
% --- global step ---

disp('attempting local-global step\n');


if old_imp 
    Lw = GI'*[sparse_diag(a11),sparse_diag(a12);sparse_diag(a12),sparse_diag(a22)]*GI;
    lhs = S'*Lw*S;
    rhs = - S'*Lw*R * TVB;
    
    U = lhs \ rhs;
    
    W = S * U + R * TVB;
    
    u = W(:,1);
    v = W(:,2);
    
    trace( W' * Lw * W / 2 );

else

    uvc = V(:,1) + 1i * V(:,2); 
    uvc(known,:) = BC(:,1) + 1i * BC(:,2);

    grad_xy = reuse.grad_xy; 
    div_xy = reuse.div_xy;
    q = reuse.q;
    Q = reuse.Q;

    Auqkuqk_lower = reuse.RL.asb_lower(mesh.GIS, a11,a12,a22); 
    Auu_qq_lower = Auqkuqk_lower(1:length(q),1:length(q));
    Auqk = Auqkuqk_lower(length(q)+1:end,1:length(q))';

    [Lcq,~] =lchol(Auu_qq_lower); 

    Asol = @(bb) Q * (Lcq' \ (Lcq \ ( Q' * bb )));

    Aqq_sol = @(bb) (Lcq' \ (Lcq \ ( bb )));

    
    uvc_uqk = Aqq_sol(- Auqk * (BC(:,1) + 1i * BC(:,2)) );
    uvc(unknown(q),:) = uvc_uqk; 

    u = real(uvc);
    v = imag(uvc);
end

if any(isnan(u)) || any(isnan(v))
    warning('local-global: UV has NaN\n');
    break; 
end


% render meshes.
% render_mesh2(W,F,'EdgeColor',[0,0,0]);
% axis equal;


% -- local step ---


if old_imp 

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

else

    [Gx_uvc, Gy_uvc] = grad_xy(uvc); 

    % the fast implementation. 
    Gxu = real(Gx_uvc); 
    Gxv = imag(Gx_uvc); 
    Gyu = real(Gy_uvc);
    Gyv = imag(Gy_uvc); 

    ccc = Area ./ abs( Gxu.*Gyv - Gxv.*Gyu );
    
    a22 = (Gxu .* Gxu + Gxv .* Gxv) .* ccc;
    a12 = - (Gxu .* Gyu + Gxv .* Gyv) .* ccc; 
    a11 = (Gyu .* Gyu + Gyv .* Gyv) .* ccc;

    newArea = Area .* ( Gxu.*Gyv - Gxv.*Gyu );

    % newArea2 = doublearea([u,v], F)/2; 
    % assert( norm(newArea - newArea2 ) < 1e-8 );

end



flipped = newArea < 0;
% number of flipped triangles. 
num_flipped = numel(find(flipped));
fprintf('flipps=%d, min_area=%g',num_flipped,min(newArea));

if num_flipped==0
    break;
end
%%
end

