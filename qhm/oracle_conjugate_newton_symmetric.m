function [E, g, Hess, out] = oracle_conjugate_newton_symmetric(mesh, at, BCBN, tp, reuse)
%%
sparse_diag = @(ddd) sparse(1:numel(ddd), 1:numel(ddd), ddd, numel(ddd), numel(ddd));
V = mesh.V;
F = mesh.F;

n = size(V,1);
dim = 2;

BC = BCBN(:,1:2);
BN = BCBN(:,3:4);

alpha = tp.alpha;
beta = tp.beta; 
gamma = tp.gamma; 

f = size(F,1);
if size(at,1)~=f
   at = reshape(at,[f,numel(at)/f]);
end

out = struct();
out.stop_sign = false;
%%

para = struct();
% coefficients for constraints
para.CC = 10;
% coefficients for energy
para.CE = 1; 
% coefficients for regularizer.
para.CR = 0; %

%
%mesh = triangle_mesh(V,F);
%mesh2 = triangle_mesh(VT,F);


G = mesh.GI;

% not used.
% Area = mesh.FA;
Area = mesh.AI;

f = mesh.f;

known = mesh.IKB;
unknown = mesh.IUB;
row_slice_matrix = @(indices) sparse(indices,1:numel(indices),1,n,numel(indices));
R = row_slice_matrix(known);
S = row_slice_matrix(unknown);

nb = numel(known);
%%

old_imp = false;


assert(tp.conj_vmap==false);

if old_imp
    u = V(:,1);
    v = V(:,2);
    u(known,:) = BC(:,1);
    v(known,:) = BC(:,2);
else
    uvc = V(:,1) + 1i * V(:,2); 
    uvc(known,:) = BC(:,1) + 1i * BC(:,2);
end
%% setup parameterizations of A

% [s_at2au,s_at2au_mul,s_pdapdt_lmul] = para_fun_old(Area,dim,'diag');

% [s_at2au_v,s_at2au_mul_v,s_pdapdt_lmul_v] = para_fun_old(Area,dim,'invdiag');



para_type = 'diag'; % 'diag-no-mass';



if old_imp

%[s_at2au,s_pdapdt_lmul,s_at2au_v,s_pdapdt_lmul_v] = ...
%    para_fun(Area,dim,para_type);

s_at2au = tp.s_at2au;
s_pdapdt_lmul = tp.s_pdapdt_lmul;
else

end

%%
%ua = da .* [Area;Area];
%va = 1./da .* [Area;Area];

if old_imp

    au = s_at2au(reshape(at,[f,numel(at)/f]));

else
    
    [au, paupat] = tensor_para_faster(at,tp.para_type);
    au = au .* Area;
    
end
    
if old_imp    

    if size(au,2)==1
        au = reshape(au, [f,3]);
    end    
    
    
    DiagA = [sparse_diag(au(:,1)), sparse_diag(au(:,2));...
             sparse_diag(au(:,2)), sparse_diag(au(:,3));];
    
    A = (G)' * DiagA * G;
    
    % if false
    % thres = 100;
    % CA = A; 
    % phind = find(CA>thres);
    % bhind = find(CA<-thres);
    % if numel(phind)>0 || numel(bhind)>0
    %     fprintf("warning: extreme Laplacian entries are capped!");
    % end
    % CA(phind) =  thres; 
    % CA(bhind) = -thres;
    % 
    % %
    % CA = CA - sparse(1:n,1:n,diag(CA),n,n);
    % % diagonal is minus sum of offdiagonal entries
    % CA = CA - sparse(1:n,1:n,sum(CA,2),n,n);
    % 
    % A = CA;
    % end
    
    Auu = A(unknown,unknown);
    Auk = A(unknown,known);

    update = true;
    switch reuse.mode 
        case 'none'
            update = true;
        case 'periodic'
            update = (0==mod(reuse.count,reuse.period));
        otherwise
            error('unknown\n');
    end
    
    
    % reuse.inited == false || isempty(reuse.Lc)
    
    if update 
        % reuse.inited = true;
        % disp('initing reuse\n');
    
        % make solver reusable
    
        [Lc,~,Qc] = chol(Auu,'lower');
    
        reuse.Lc = Lc;
        reuse.Qc = Qc;
    else
        disp('reusing Lc, Qc.\n')
        
        Lc = reuse.Lc;
        Qc = reuse.Qc;
    end
    reuse.count = reuse.count + 1;
    
    Asol = @(bb) Qc * (Lc' \ (Lc \ ( Qc' * bb )));


else
%%


    if isempty(reuse.RL)

        [L_sym] = symbolic_lap(n,F); 
    
        q = analyze(L_sym(unknown,unknown),'sym')'; % should be the same as: q = analyze(Auu,'sym')';
        % so q is a column vector. 
            
        Q = sparse(q,1:numel(q),ones(numel(q),1),numel(q),numel(q)); 
        % such that Q'*x=x(q). 
    
        iq = zeros(size(q));
        iq(q) = (1:length(q))'; % inverse permuting of q. 
    
    
        % so that the following is 0: sum(sum(abs( Auu(q,q) - A(unknown(q),unknown(q)) )))
        
    
        p_all = [unknown(q);known]; 
        
        ip_all = zeros(size(p_all));
        ip_all(p_all) = (1:length(p_all))'; 
    
        perm_F = @(invq, FF) [invq(FF(:,1)),invq(FF(:,2)),invq(FF(:,3))];
    
        FP = perm_F(ip_all, F); 
    
        % [Auu_qq_lower, Auqkuqk] = assemble_lap(n,mesh.GIS,FP,au(:,1),au(:,2),au(:,3));
        
        GGx = G(1:f,:); 
        GGy = G(f+1:2*f,:);  

        grad_xy = @(uu) deal(GGx * uu, GGy * uu);
        div_xy = @(ux,uy) GGx' * ux + GGy' * uy;

        reuse.grad_xy = grad_xy; 
        reuse.div_xy = div_xy;

        reuse.RL = assemble_lap_core( n,FP); 
        reuse.q = q;
        reuse.Q = Q;
    else
        grad_xy = reuse.grad_xy; 
        div_xy = reuse.div_xy;
        q = reuse.q;
        Q = reuse.Q;

    end



    % Auqkuqk = reuse.RL.asb_full( mesh.GIS,au(:,1),au(:,2),au(:,3)); 

    Auqkuqk_lower = reuse.RL.asb_lower(mesh.GIS, au(:,1),au(:,2),au(:,3)); 

    Auu_qq_lower = Auqkuqk_lower(1:length(q),1:length(q));
    % Auu_qq = Auqkuqk(1:length(q),1:length(q));

    % Auqk = Auqkuqk(1:length(q),length(q)+1:end);
    Auqk = Auqkuqk_lower(length(q)+1:end,1:length(q))';

    % [~, A] = assemble_lap(n,mesh.GIS,mesh.F,au(:,1),au(:,2),au(:,3));
    

    if false    
    
        Auu = A(unknown,unknown);
        Auk = A(unknown,known);
        
        Auu_qq = Auu(q,q); 

    end

    % [Lc,~] =lchol(Auu_qq);
    [Lcq,~] =lchol(Auu_qq_lower); 

    Asol = @(bb) Q * (Lcq' \ (Lcq \ ( Q' * bb )));

    Aqq_sol = @(bb) (Lcq' \ (Lcq \ ( bb )));

    % to debug, double check that:  assert(norm( Asol(bbb) - Auu\bbb )<1e-8)
end

ATAmul = @(xxx) [ au(:,1).* xxx(1:f) + au(:,2).* xxx(f+1:2*f);  au(:,2).* xxx(1:f) + au(:,3).* xxx(f+1:2*f) ];


if old_imp
    u(unknown,:) = - Auu \ (Auk * u(known,:));
    v(unknown,:) = - Auu \ (Auk * v(known,:));

% u(unknown,:) = Asol(-Auk * u(known,:));
% v(unknown,:) = Asol(-Auk * v(known,:));

% uvc_uk = Asol(- Auk * (u(known,:)+1i*v(known,:)) );
% u(unknown,:) = real(uvc_uk); 
% v(unknown,:) = imag(uvc_uk); 

% uvc_uqk = Aqq_sol(- Auqk * (u(known,:)+1i*v(known,:)) );
% u(unknown(q),:) = real(uvc_uqk); 
% v(unknown(q),:) = imag(uvc_uqk); 
else
    uvc_uqk = Aqq_sol(- Auqk * uvc(known,:) );
    uvc(unknown(q),:) = uvc_uqk; 
end
% the span op
d = 2;
% this is slow
span = @(xx) symmetric_tensor_span(xx,d);

span_dot = @(xx,yy) symmetric_tensor_span_dot(xx,yy,d); 

% span_dot = @(xx,yy) span(xx)' *  yy; 



%% weighted Dirichlet energy etc. 
% g_u = s_pdapdt_lmul( reshape(at,[f,numel(at)/f]), ...
%     0.5 * span( G * u)' * G * u -  (span(G*u)' * G * S) * (Asol(S'*A*u)) ...
%     );
% 
% %
% g_v = s_pdapdt_lmul( reshape(at,[f,numel(at)/f]), ...
%     0.5 * span( G * v)' * G * v -  (span(G*v)' * G * S) * (Asol(S'*A*v)) ...
%     );

% Gu = G * u;
% g_u = s_pdapdt_lmul( at, ...
%     0.5 * span_dot( Gu, Gu) -  span_dot( Gu,  (G * (S * Asol(S'*(A*u)))) ) ...
%     );
% %
% Gv = G * v;
% g_v = s_pdapdt_lmul( at, ...
%     0.5 * span_dot( Gv, Gv) -  span_dot( Gv,  (G * (S * Asol(S'*(A*v)))) ) ...
%     );


if old_imp
Gu = G * u;
Gv = G * v;
Gux = Gu(1:f,:); Guy = Gu(f+1:2*f,:);  
Gvx = Gv(1:f,:); Gvy = Gv(f+1:2*f,:);

else
% Gx = G(1:f,:); 
% Gy = G(f+1:2*f,:);  

% Guvc = G * uvc;
% Gu = real(Guvc);
% Gv = imag(Guvc);

% Gx_uvc = Gx * uvc;
% Gy_uvc = Gy * uvc;

[Gx_uvc, Gy_uvc] = grad_xy(uvc); 

end

span_dot_xy = @(bbx, bby, ccx, ccy) [bbx.*ccx; bby.*ccx+bbx.*ccy; bby.*ccy]; 

if old_imp

    res_u = Asol(S'*(A*u));
    res_v = Asol(S'*(A*v));
    
    % complex solve saves one linear solve.
    % res_uvc = Asol(S'*(A*(u+1i*v)));  


    % res_uvc = Asol(S'*(G'*(ATAmul(Guvc))));  
    % res_u = real(res_uvc);
    % res_v = imag(res_uvc);
    
    gn_u = 0.5 * span_dot( Gu, Gu) -  span_dot( Gu,  (G * (S * res_u )) ); 
    gn_v = 0.5 * span_dot( Gv, Gv) -  span_dot( Gv,  (G * (S * res_v )) );

else

    if false

        % GtAGuvc = G'*(ATAmul(Guvc)); 
    
        % GtAGuvc = Gx' * (au(:,1) .* Gx_uvc + au(:,2) .* Gy_uvc) ...
        %         + Gy' * (au(:,2) .* Gx_uvc + au(:,3) .* Gy_uvc); 
    
        GtAGuvc = div_xy( au(:,1) .* Gx_uvc + au(:,2) .* Gy_uvc,  ...
                          au(:,2) .* Gx_uvc + au(:,3) .* Gy_uvc ); 
    
        q_res_uvc = Aqq_sol(GtAGuvc(unknown(q),:));  
        
        % q_res_u = real(q_res_uvc);
        % q_res_v = imag(q_res_uvc); 
        % gn_u = 0.5 * span_dot( Gu, Gu) -  span_dot( Gu,  (G *(S*(Q*(q_res_u )))) ); 
        % gn_v = 0.5 * span_dot( Gv, Gv) -  span_dot( Gv,  (G *(S*(Q*(q_res_v )))) );
    
    
        % rrr = (G *(S*(Q*(q_res_uvc )))) ; 
        % gn_u = 0.5 * span_dot( Gu, Gu) -  span_dot( Gu, real(rrr) ); 
        % gn_v = 0.5 * span_dot( Gv, Gv) -  span_dot( Gv, imag(rrr) );
    
        % gn_u = span_dot( Gu, 0.5 * Gu - real(rrr) ); 
        % gn_v = span_dot( Gv, 0.5 * Gv - imag(rrr) );

        
        tt_uvc = (S*(Q*(q_res_uvc )));
        %tx_uvc = Gx * tt_uvc; 
        %ty_uvc = Gy * tt_uvc; 
        [tx_uvc, ty_uvc] = grad_xy(tt_uvc); 
    
        % rux = 0.5 * Gux - real(tx_uvc); 
        % rvx = 0.5 * Gvx - imag(tx_uvc); 
        % ruy = 0.5 * Guy - real(ty_uvc);
        % rvy = 0.5 * Gvy - imag(ty_uvc);
    
        rx = 0.5 * Gx_uvc - tx_uvc;
        ry = 0.5 * Gy_uvc - ty_uvc;
    
        % gn_u = span_dot_xy(Gux,Guy,rux,ruy); 
        % gn_v = span_dot_xy(Gvx,Gvy,rvx,rvy); 
    
        gn_u = span_dot_xy(real(Gx_uvc),real(Gy_uvc),real(rx),real(ry)); 
        gn_v = span_dot_xy(imag(Gx_uvc),imag(Gy_uvc),imag(rx),imag(ry)); 

    else

        Gxu = real(Gx_uvc); 
        Gyu = real(Gy_uvc); 
        Gxv = imag(Gx_uvc); 
        Gyv = imag(Gy_uvc); 

        gn_u = 0.5 * span_dot_xy(Gxu,Gyu,Gxu,Gyu); 
        gn_v = 0.5 * span_dot_xy(Gxv,Gyv,Gxv,Gyv); 

        GtAGuvc = div_xy( au(:,1) .* Gx_uvc + au(:,2) .* Gy_uvc,  ...
                  au(:,2) .* Gx_uvc + au(:,3) .* Gy_uvc ); 

    end
       
end

if old_imp
    out.newArea = doublearea([u,v], F); 
else
    out.newArea = Area .* ( Gxu.*Gyv - Gxv.*Gyu );
end

flipped = out.newArea < 0;
out.num_flipped = numel(find(flipped));

out.stop_sign = out.num_flipped==0;

if old_imp
    g_u = s_pdapdt_lmul( at, gn_u);
    g_v = s_pdapdt_lmul( at, gn_v);
    g_drlt = g_u + g_v;
else
    % g_u = paupat.s_pdapdt( gn_u, Area);
    % g_v = paupat.s_pdapdt( gn_v, Area);
    % g_drlt = g_u + g_v;

    g_drlt = paupat.s_pdapdt( gn_u + gn_v, Area);

end



assert(tp.conj_vmap==false);

if old_imp

    E_u = 0.5 * u' * A * u - 0.5 * BC(:,1)' * BN(:,1);
    E_v = 0.5 * v' * A * v - 0.5 * BC(:,2)' * BN(:,2);

    % E_u = 0.5 * Gu' * ATAmul( Gu ) - 0.5 * BC(:,1)' * BN(:,1);
    % E_v = 0.5 * Gv' * ATAmul( Gv ) - 0.5 * BC(:,2)' * BN(:,2);
    E_drlt = E_u + E_v; 

else
    % E_drlt = 0.5 * Guvc' * ATAmul( Guvc ) - 0.5 * BC(:,1)' * BN(:,1) - 0.5 * BC(:,2)' * BN(:,2); 
    
    % since "'" does the Hermitian conjungate
    E_drlt = 0.5 * uvc' * GtAGuvc - 0.5 * BC(:,1)' * BN(:,1) - 0.5 * BC(:,2)' * BN(:,2); 
    E_drlt = real(E_drlt); % possible small error. 

end

% beta term

g_bn = 0;
E_bn = 0;

if beta~=0

RE = G * (A * [u,v] - R * BN);
RW = RE - G * S * (Asol(S'*G'*(DiagA*RE)));

% g_bn = s_pdapdt_lmul( at, span( G * u)' * RW(:,1)) ... 
%      + s_pdapdt_lmul( at, span( G * v)' * RW(:,2));

g_bn = s_pdapdt_lmul( at, span_dot( G * u, RW(:,1))) ... 
     + s_pdapdt_lmul( at, span_dot( G * v, RW(:,2)));

E_bn = 0.5 * sum(sum(( A * [u,v] - R * BN ).^2)); % todo: fill in it

end

% gamma term 
 
g_wbn = 0;
E_wbn = 0;

if gamma~=0

XX = [A,ones(n,1); ones(1,n),1] \ [R*BN; zeros(1,size(BN,2))];
PW = XX(1:n,:);

% numerical unstable: 
% PW = A\(R*BN); % pseudo inverse, PW is unique up to a constant. 

GPW = G * PW;

% TODO: speed the implementation.

% g_wbn =  s_pdapdt_lmul( at, 0.5 * span( G * u)' * (G * u) -  span(G*u)' * (G * (S * (Asol(S'*(A*u))))) ) ...
%        + s_pdapdt_lmul( at, 0.5 * span( G * v)' * (G * v) -  span(G*v)' * (G * (S * (Asol(S'*(A*v))))) ) ... % assume A_u==A_v here
%        + s_pdapdt_lmul( at, -0.5 * span( GPW(:,1))' * GPW(:,1)) ... 
%        + s_pdapdt_lmul( at, -0.5 * span( GPW(:,2))' * GPW(:,2));

g_wbn =  s_pdapdt_lmul( at, 0.5 * span_dot( G * u, G * u) -  span_dot(G*u, G * (S * (Asol(S'*(A*u))))) ) ...
       + s_pdapdt_lmul( at, 0.5 * span_dot( G * v, G * v) -  span_dot(G*v, G * (S * (Asol(S'*(A*v))))) ) ... % assume A_u==A_v here
       + s_pdapdt_lmul( at, -0.5 * span_dot( GPW(:,1), GPW(:,1)) ) ... 
       + s_pdapdt_lmul( at, -0.5 * span_dot( GPW(:,2), GPW(:,2)) );
   
E_wbn =  0.5 * u' * A * u + 0.5 * v' * A * v + ...
    0.5 * trace(PW(:,1:2)'*A*PW(:,1:2));
end

%%

if old_imp

else
    % can move earlier. 
    u = real(uvc);
    v = imag(uvc); 
end


E = alpha * (E_drlt) + E_bn * beta + E_wbn * gamma;
g = alpha * (g_drlt) + g_bn * beta + g_wbn * gamma;

fprintf("\t\t E_drlt=%g,\t E_bn=%g,\t E_wbn=%g,\t ", E_drlt, E_bn, E_wbn);

if ~isempty(tp.energy_uv)

[er_uv, r_uv] = tp.energy_uv([u,v]);

r_u = r_uv(1:n);
r_v = r_uv(n+1:2*n);

if old_imp

    gr_u = s_pdapdt_lmul( at, ...
        - span_dot( G*u,  (G * (S * Asol(S'*r_u))) ) ...
        );
    
    %
    gr_v = s_pdapdt_lmul( at, ...
        - span_dot( G*v,  (G * (S * Asol(S'*r_v))) ) ...
        );

    gr_uv = gr_u + gr_v;

else
    [trx_uvc, try_uvc] = grad_xy((S * Asol(S'*(r_u+1i*r_v)))); 

    % gr_uv = paupat.s_pdapdt( ...
    %     - span_dot( Gu, real([trx_uvc;try_uvc])  ) ...
    %     - span_dot( Gv, imag([trx_uvc;try_uvc])  ), ...
    %     Area);

    % er_uv = 0; 
    gr_uv = 1.000 * paupat.s_pdapdt( ...
         - span_dot_xy(real(Gx_uvc),real(Gy_uvc),real(trx_uvc),real(try_uvc))...
         - span_dot_xy(imag(Gx_uvc),imag(Gy_uvc),imag(trx_uvc),imag(try_uvc)),...  
        Area);

end


E = E + er_uv;
g = g + gr_uv; 

fprintf("Er_uv=%g,\t \t", er_uv); % \n

else

fprintf("\n");

end


%%

reg = tp.reg_value(at);
% + tp.rau_value(au);

if tp.rau_value(au)~=0
    error('not implemented!\n');
end

fprintf("Reg=%g,\t \t", reg); 

E = E + reg;
g = g + tp.reg_grad(at);

% + s_pdapdt_lmul( at, tp.rau_grad(au)); % not implemented yet.

Hess = [];

out.u = u;
out.v = v;

if ~isempty(reuse.fine_recorder)
    if isa(reuse.fine_recorder,'logical') && reuse.fine_recorder==true
        tt = 1;
        reuse.fine_recorder = []; % it soon becomes not empty. 
    else
        tt = 1 + numel(reuse.fine_recorder);
    end
    reuse.fine_recorder(tt).u = u;
    reuse.fine_recorder(tt).v = v;
end

%%
%
% render_mesh2([u,v],F,'EdgeColor',[0,0,0]);
% axis equal;

% dda = 1e-7 * normrnd(0,1,size(da));
% da = da_old + dda;

%%
% if false
% Aef = edge_triangle_adjacency(F,edges(F));
% Aef = Aef(find(Aef(:,1)>0 & Aef(:,2)>0),:);
% Def = sparse(1:size(Aef,1), Aef(:,1), 1, size(Aef,1),f) - sparse(1:size(Aef,1), Aef(:,2), 1, size(Aef,1),f);
% L_dual = Def' * Def;
% 
% end

%%

% g = g * 0.1 + 0.9 * mean(g);
%%

% da = da - 0.001 *  (0.001*speye(numel(da)) + Hess_v) \ (g_u + g_v);
% da = da - 0.01 * (g_u + g_v);
