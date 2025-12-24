function [E, g, Hess, out] = oracle_conjugate_newton_symmetric(mesh, at, BCBN, tp)
%%
sparse_diag = @(ddd) sparse(1:numel(ddd), 1:numel(ddd), ddd, numel(ddd), numel(ddd));
V = mesh.V;
F = mesh.F;

dim = 2;

BC = BCBN(:,1:2);
BN = BCBN(:,3:4);
%%
BE = outline(F);
B = BE(:,1);

n = size(V,1);
nb = numel(B);

% eq_lhs = [  sparse(1:nb, B, 1, nb, n*2);...
%             sparse(1:nb, B+n, 1, nb, n*2)];
% eq_rhs = [ VT(B,1); VT(B,2)];
% 
% if false
% %
% [V2,F2,~,~] = subdivide_with_constraint(V,F,eq_lhs,eq_rhs,1);
% [VT2,~,~,~] = subdivide_with_constraint(VT,F,eq_lhs,eq_rhs,1);
% 
% V = V2;
% F = F2;
% VT = VT2;
% end

BE = outline(F);
B = BE(:,1);

n = size(V,1);
nb = numel(B);

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

d01 = mesh.DEC.d01;

G = mesh.GI;

n = mesh.n;
ne = mesh.ne;
f = mesh.f;
ft = ones(1,f);

known = mesh.B;
unknown = mesh.UB;
row_slice_matrix = @(indices) sparse(indices,1:numel(indices),1,n,numel(indices));
R = row_slice_matrix(known);
S = row_slice_matrix(unknown);


Area = doublearea(V,F)/2;



%%


if tp.conj_vmap

    DH = [sparse_diag(zeros(f,1)), -sparse_diag(ones(f,1)); ...
        sparse_diag(ones(f,1)), sparse_diag(zeros(f,1))];

else

    DH = speye(2*f);
    
end


u = V(:,1);
v = V(:,2);
u(known,:) = BC(:,1);
v(known,:) = BC(:,2);
%% setup parameterizations of A

% [s_at2au,s_at2au_mul,s_pdapdt_lmul] = para_fun_old(Area,dim,'diag');

% [s_at2au_v,s_at2au_mul_v,s_pdapdt_lmul_v] = para_fun_old(Area,dim,'invdiag');



para_type = 'diag'; % 'diag-no-mass';

%[s_at2au,s_pdapdt_lmul,s_at2au_v,s_pdapdt_lmul_v] = ...
%    para_fun(Area,dim,para_type);

s_at2au = tp.s_at2au;
s_pdapdt_lmul = tp.s_pdapdt_lmul;

if tp.conj_vmap

    s_at2au_v = tp.s_at2au_v;
    s_pdapdt_lmul_v = tp.s_pdapdt_lmul_v;

else
    s_at2au_v = tp.s_at2au;
    s_pdapdt_lmul_v = tp.s_pdapdt_lmul;
        
end



%%
%ua = da .* [Area;Area];
%va = 1./da .* [Area;Area];

au = s_at2au(reshape(at,[f,numel(at)/f]));

AA = zeros(f,2,2);
IA = zeros(f,2,2);

a11 = au(1:f);
a12 = au(f+1:2*f);
a22 = au(2*f+1:end);

AA(:,1,1) = a11(:);
AA(:,1,2) = a12(:);
AA(:,2,1) = a12(:);
AA(:,2,2) = a22(:);

av = s_at2au_v(reshape(at,[f,numel(at)/f]));

av11 = av(1:f);
av12 = av(f+1:2*f);
av22 = av(2*f+1:end);

IA(:,1,1) = av11(:);
IA(:,1,2) = av12(:);
IA(:,2,1) = av12(:);
IA(:,2,2) = av22(:);

A_u = (G)' * [sparse_diag(AA(:,1,1)), sparse_diag(AA(:,1,2));...
              sparse_diag(AA(:,2,1)), sparse_diag(AA(:,2,2));] * G;
A_v = (DH * G)' * [sparse_diag(IA(:,1,1)), sparse_diag(IA(:,1,2));...
                   sparse_diag(IA(:,2,1)), sparse_diag(IA(:,2,2));] * DH * G;

u(unknown,:) = - A_u(unknown,unknown) \ (A_u(unknown,known) * u(known,:));
v(unknown,:) = - A_v(unknown,unknown) \ (A_v(unknown,known) * v(known,:));

% the span op
d = 2;
span = @(xx) symmetric_tensor_span(xx,d);

g_u = s_pdapdt_lmul( reshape(at,[f,numel(at)/f]), ...
    0.5 * span( G * u)' * G * u -  (span(G*u)' * G * S) * ((S'*A_u*S)\(S'*A_u*u)) ...
    );

%
g_v = s_pdapdt_lmul_v ( reshape(at,[f,numel(at)/f]), ...
        0.5 * span(DH*G * v)' * DH*G * v ...
        -  span(DH*G*v)' *DH * G * S * ((S'*A_v*S)\(S'*A_v*v)) ...
        );

g_bn = 0;
E_bn = 0;
    
if false
assert(tp.conj_vmap==false);


DiagA = [sparse_diag(AA(:,1,1)), sparse_diag(AA(:,1,2));...
              sparse_diag(AA(:,2,1)), sparse_diag(AA(:,2,2));];

RE = G * (A_u * [u,v] - R * BN);
RW = RE - G * S * ((S'*A_u*S)\(S'*G'*(DiagA*RE)));

g_bn = s_pdapdt_lmul( reshape(at,[f,numel(at)/f]), span( G * u)' * RW(:,1)) ... 
     + s_pdapdt_lmul( reshape(at,[f,numel(at)/f]), span( G * v)' * RW(:,2));

E_bn = 0.5 * sum(sum(( A_u * [u,v] - R * BN ).^2)) ; % todo: fill in it
 
PW = (G'*DiagA*G) \ (R*BN); % pseudo inverse, PW is unique up to a constant. 
GPW = G * PW;

g_wbn =  s_pdapdt_lmul( reshape(at,[f,numel(at)/f]), 0.5 * span( G * u)' * G * u -  (span(G*u)' * G * S) * ((S'*A_u*S)\(S'*A_u*u))) ...
       + s_pdapdt_lmul( reshape(at,[f,numel(at)/f]), 0.5 * span( G * v)' * G * v -  (span(G*v)' * G * S) * ((S'*A_u*S)\(S'*A_u*v))) ... % assume A_u==A_v here
      + s_pdapdt_lmul( reshape(at,[f,numel(at)/f]), -0.5 * span( GPW(:,1))' * GPW(:,1)) ... 
      + s_pdapdt_lmul( reshape(at,[f,numel(at)/f]), -0.5 * span( GPW(:,2))' * GPW(:,2));

E_wbn =  0.5 * u' * A_u * u + 0.5 * v' * A_u * v + ...
    0.5 * trace(PW(:,1:2)'*A_u*PW(:,1:2));
end
%
% Hess_v = sparse_diag( (DH * G * v).^2 ./ da.^3 ) ...
%     + sparse_diag(1./da.^2) * sparse_diag(DH * G * v) * ...
%     DH * G * S * inv( S' * (DH * G)' * sparse_diag(1./da) * DH * G * S ) * S' * G' * DH' * ...
%      sparse_diag(DH * G * v) * sparse_diag(1./da.^2);
% 
% Hess_u = - sparse_diag(G * u) * ...
%      G * S * inv( S' * (G)' * sparse_diag(da) * G * S ) * S' * G' * ...
%      sparse_diag( G * u);

E_u = 0.5 * u' * A_u * u;
E_v = 0.5 * v' * A_v * v;


beta = 0;
alpha = 1;

E = alpha * (E_u + E_v) + E_bn * beta;

%
% render_mesh2([u,v],F,'EdgeColor',[0,0,0]);
% axis equal;

% dda = 1e-7 * normrnd(0,1,size(da));
% da = da_old + dda;
%%

% da = da - 0.001 *  (0.001*speye(numel(da)) + Hess_v) \ (g_u + g_v);
% da = da - 0.01 * (g_u + g_v);

%%
g = alpha * (g_u + g_v) + g_bn * beta;

g = g + tp.reg_grad(at);

% Hess = Hess_u + Hess_v;
Hess = [];


out = struct();
out.u = u;
out.v = v;

end
%%

function [s_at2au,s_pdapdt_lmul,s_at2au_v,s_pdapdt_lmul_v] = para_fun(FA,dim,para_type)

    Area = FA;
    f = size(FA,1);
    
    assert(dim==2);

    switch para_type

        case 'diag-no-mass'

            % without mass matrix. 

            s_at2au = @(aatt) [aatt(:,1);zeros(f,1);aatt(:,3)];
            s_pdapdt_lmul = @(aatt,gg) [ones(f,1);zeros(f,1);ones(f,1)].* gg;

            s_at2au_v = @(aatt) [1./aatt(:,1);zeros(f,1);1./aatt(:,3)];
            s_pdapdt_lmul_v = @(aatt,gg)  [-1./aatt(:,1).^2;zeros(f,1);-1./aatt(:,3).^2].* gg;

        case 'diag'
            if false
            %%
            syms aA1(p1,p2,p3);
            syms aA2(p1,p2,p3);
            syms aA3(p1,p2,p3);
                  
            syms bA1(p1,p2,p3);
            syms bA2(p1,p2,p3);
            syms bA3(p1,p2,p3);

            % diag:
            
            aA1(p1,p2,p3) = p1;
            aA2(p1,p2,p3) = 0;
            aA3(p1,p2,p3) = p3;
            
            bA1(p1,p2,p3) = 1/p1;
            bA2(p1,p2,p3) = 0;
            bA3(p1,p2,p3) = 1/p3;

            % diag-sq:
            
            aA1(p1,p2,p3) = p1.^2;
            aA2(p1,p2,p3) = 0;
            aA3(p1,p2,p3) = p3.^2;
            
            bA1(p1,p2,p3) = 1/p1.^2;
            bA2(p1,p2,p3) = 0;
            bA3(p1,p2,p3) = 1/p3.^2;
            
            % complex-sym:

            aA1(p1,p2,p3) = -((p1 - 1)^2 - abs(p3)^2 + p2^2)/(p1^2 - abs(p3 + 1)^2 + p2^2);
            aA2(p1,p2,p3) = (2*p2)/(p1^2 - abs(p3 + 1)^2 + p2^2);
            aA3(p1,p2,p3) = -((p1 + 1)^2 - abs(p3)^2 + p2^2)/(p1^2 - abs(p3 + 1)^2 + p2^2);
 
            bA1(p1,p2,p3) = ((p1^2 - abs(p3 + 1)^2 + p2^2)*(2*p1 - abs(p3)^2 + p1^2 + p2^2 + 1))/(2*abs(p3)^2 - abs(p3)^4 + 2*p1^2*abs(p3)^2 + 2*p2^2*abs(p3)^2 + 2*p1^2 + 2*p2^2 - p1^4 - p2^4 - 2*p1^2*p2^2 - 1);
            bA2(p1,p2,p3) = (2*p2*(p1^2 - abs(p3 + 1)^2 + p2^2))/(2*abs(p3)^2 - abs(p3)^4 + 2*p1^2*abs(p3)^2 + 2*p2^2*abs(p3)^2 + 2*p1^2 + 2*p2^2 - p1^4 - p2^4 - 2*p1^2*p2^2 - 1);
            bA3(p1,p2,p3) = ((p1^2 - abs(p3 + 1)^2 + p2^2)*(p1^2 - abs(p3)^2 - 2*p1 + p2^2 + 1))/(2*abs(p3)^2 - abs(p3)^4 + 2*p1^2*abs(p3)^2 + 2*p2^2*abs(p3)^2 + 2*p1^2 + 2*p2^2 - p1^4 - p2^4 - 2*p1^2*p2^2 - 1);

            % LLT:
            
            aA1(p1,p2,p3) = p1.^2;
            aA2(p1,p2,p3) = p1.*p2;
            aA3(p1,p2,p3) = p2.^2 + p3.^2;
            
            bA1(p1,p2,p3) = (p2^2 + p3^2)/(p1^2*p3^2);
            bA2(p1,p2,p3) = -p2/(p1*p3^2);
            bA3(p1,p2,p3) = 1/p3^2;
            
            mA1 = matlabFunction(aA1);
            mA2 = matlabFunction(aA2);
            mA3 = matlabFunction(aA3);
            
            s_dAdP11 = matlabFunction(diff(aA1,p1));
            s_dAdP12 = matlabFunction(diff(aA1,p2));
            s_dAdP13 = matlabFunction(diff(aA1,p3));

            s_dAdP21 = matlabFunction(diff(aA2,p1));
            s_dAdP22 = matlabFunction(diff(aA2,p2));
            s_dAdP23 = matlabFunction(diff(aA2,p3));

            s_dAdP31 = matlabFunction(diff(aA3,p1));
            s_dAdP32 = matlabFunction(diff(aA3,p2));
            s_dAdP33 = matlabFunction(diff(aA3,p3));
            
            disp('mA1=...'),disp(mA1);
            disp('mA2=...'),disp(mA2);
            disp('mA3=...'),disp(mA3);
            
            disp('s_dAdP11=...'),disp(s_dAdP11);
            disp('s_dAdP12=...'),disp(s_dAdP12);
            disp('s_dAdP13=...'),disp(s_dAdP13);
            
            disp('s_dAdP21=...'),disp(s_dAdP21);
            disp('s_dAdP22=...'),disp(s_dAdP22);
            disp('s_dAdP23=...'),disp(s_dAdP23);
                        
            disp('s_dAdP31=...'),disp(s_dAdP31);
            disp('s_dAdP32=...'),disp(s_dAdP32);
            disp('s_dAdP33=...'),disp(s_dAdP33);
            
            %

            nA1 = matlabFunction(bA1);
            nA2 = matlabFunction(bA2);
            nA3 = matlabFunction(bA3);
            
            t_dAdP11 = matlabFunction(diff(bA1,p1));
            t_dAdP12 = matlabFunction(diff(bA1,p2));
            t_dAdP13 = matlabFunction(diff(bA1,p3));

            t_dAdP21 = matlabFunction(diff(bA2,p1));
            t_dAdP22 = matlabFunction(diff(bA2,p2));
            t_dAdP23 = matlabFunction(diff(bA2,p3));

            t_dAdP31 = matlabFunction(diff(bA3,p1));
            t_dAdP32 = matlabFunction(diff(bA3,p2));
            t_dAdP33 = matlabFunction(diff(bA3,p3));
            
            disp('nA1=...'),disp(nA1);
            disp('nA2=...'),disp(nA2);
            disp('nA3=...'),disp(nA3);
            
            disp('t_dAdP11=...'),disp(t_dAdP11);
            disp('t_dAdP12=...'),disp(t_dAdP12);
            disp('t_dAdP13=...'),disp(t_dAdP13);
            
            disp('t_dAdP21=...'),disp(t_dAdP21);
            disp('t_dAdP22=...'),disp(t_dAdP22);
            disp('t_dAdP23=...'),disp(t_dAdP23);
                        
            disp('t_dAdP31=...'),disp(t_dAdP31);
            disp('t_dAdP32=...'),disp(t_dAdP32);
            disp('t_dAdP33=...'),disp(t_dAdP33);
            else
                
                % res_sym_grad_diag;
                % res_sym_grad_diag_sq;
                % res_sym_grad_complex_sym;
                res_sym_grad_llt;
            
            end

            s_11 = @(aatt) s_dAdP11(aatt(:,1),aatt(:,2),aatt(:,3));
            s_21 = @(aatt) s_dAdP21(aatt(:,1),aatt(:,2),aatt(:,3));
            s_31 = @(aatt) s_dAdP31(aatt(:,1),aatt(:,2),aatt(:,3));

            s_12 = @(aatt) s_dAdP12(aatt(:,1),aatt(:,2),aatt(:,3));
            s_22 = @(aatt) s_dAdP22(aatt(:,1),aatt(:,2),aatt(:,3));
            s_32 = @(aatt) s_dAdP32(aatt(:,1),aatt(:,2),aatt(:,3));

            s_13 = @(aatt) s_dAdP13(aatt(:,1),aatt(:,2),aatt(:,3));
            s_23 = @(aatt) s_dAdP23(aatt(:,1),aatt(:,2),aatt(:,3));
            s_33 = @(aatt) s_dAdP33(aatt(:,1),aatt(:,2),aatt(:,3));
            
            t_11 = @(aatt) t_dAdP11(aatt(:,1),aatt(:,2),aatt(:,3));
            t_21 = @(aatt) t_dAdP21(aatt(:,1),aatt(:,2),aatt(:,3));
            t_31 = @(aatt) t_dAdP31(aatt(:,1),aatt(:,2),aatt(:,3));

            t_12 = @(aatt) t_dAdP12(aatt(:,1),aatt(:,2),aatt(:,3));
            t_22 = @(aatt) t_dAdP22(aatt(:,1),aatt(:,2),aatt(:,3));
            t_32 = @(aatt) t_dAdP32(aatt(:,1),aatt(:,2),aatt(:,3));

            t_13 = @(aatt) t_dAdP13(aatt(:,1),aatt(:,2),aatt(:,3));
            t_23 = @(aatt) t_dAdP23(aatt(:,1),aatt(:,2),aatt(:,3));
            t_33 = @(aatt) t_dAdP33(aatt(:,1),aatt(:,2),aatt(:,3));
            
        
            if false
            idf = 1:f;
            idfau = 1:(3*f);
            s_pdapdt = @(aatt) sparse(...
                        [[idf+0*f;idf+0*f;idf+0*f];[idf+1*f;idf+1*f;idf+1*f];[idf+2*f;idf+2*f;idf+2*f];],...
                        [1:length(idfau),1:length(idfau),1:length(idfau)],...
                        [[s_11(aatt).*Area;s_12(aatt).*Area;s_13(aatt).*Area];...
                         [s_21(aatt).*Area;s_22(aatt).*Area;s_23(aatt).*Area];...
                         [s_31(aatt).*Area;s_32(aatt).*Area;s_33(aatt).*Area];],...
                        3*f,length(idfau));
                    
            s_pdapdt_lmul = @(aatt,gg) s_pdapdt(aatt)' * gg; % not ".*" here, critical to have transpose here. 
            % this somehows gives incorrect result, check it before use!!!
            end
            
            s_at2au = @(aatt) [...
                 Area .* mA1(aatt(:,1),aatt(:,2),aatt(:,3));...
                 Area .* mA2(aatt(:,1),aatt(:,2),aatt(:,3));...
                 Area .* mA3(aatt(:,1),aatt(:,2),aatt(:,3));...
                ];
            
            s_pdapdt_lmul = @(aatt,gg) [...
                Area.*( s_11(aatt).*gg(1:f) + s_21(aatt).*gg(f+1:2*f) + s_31(aatt).*gg(2*f+1:3*f) );...
                Area.*( s_12(aatt).*gg(1:f) + s_22(aatt).*gg(f+1:2*f) + s_32(aatt).*gg(2*f+1:3*f) );...
                Area.*( s_13(aatt).*gg(1:f) + s_23(aatt).*gg(f+1:2*f) + s_33(aatt).*gg(2*f+1:3*f) );...
                ];
            
            % this is for diagonal tensor:
            % s_pdapdt_lmul = @(aatt,gg) [Area;zeros(f,1);Area].* gg;
            
            % s_at2au_v = @(aatt) [Area./aatt(:,1);zeros(f,1);Area./aatt(:,3)];
            % s_pdapdt_lmul_v = @(aatt,gg)  [-Area./aatt(:,1).^2;zeros(f,1);-Area./aatt(:,3).^2].* gg;
            
            s_at2au_v = @(aatt) [...
                 Area .* nA1(aatt(:,1),aatt(:,2),aatt(:,3));...
                 Area .* nA2(aatt(:,1),aatt(:,2),aatt(:,3));...
                 Area .* nA3(aatt(:,1),aatt(:,2),aatt(:,3));...
                ];
            
            s_pdapdt_lmul_v = @(aatt,gg) [...
                Area.*( t_11(aatt).*gg(1:f) + t_21(aatt).*gg(f+1:2*f) + t_31(aatt).*gg(2*f+1:3*f) );...
                Area.*( t_12(aatt).*gg(1:f) + t_22(aatt).*gg(f+1:2*f) + t_32(aatt).*gg(2*f+1:3*f) );...
                Area.*( t_13(aatt).*gg(1:f) + t_23(aatt).*gg(f+1:2*f) + t_33(aatt).*gg(2*f+1:3*f) );...
                ];


        case 'diag2'

            % with mass matrix.
            % this is for diagonal tensor.

            s_at2au = @(aatt) [Area.*aatt(:,1);zeros(f,1);Area.*aatt(:,3)];
            s_pdapdt_lmul = @(aatt,gg) [Area;zeros(f,1);Area].* gg;

            s_at2au_v = @(aatt) [Area./aatt(:,1);zeros(f,1);Area./aatt(:,3)];
            s_pdapdt_lmul_v = @(aatt,gg)  [-Area./aatt(:,1).^2;zeros(f,1);-Area./aatt(:,3).^2].* gg;

        otherwise 
            error('Unsupported para type~\n');

    end

end


function [s_at2au,s_at2au_mul,s_pdapdt_lmul] = para_fun_old(FA,dim,para_type)

f = size(FA,1);
    
base = ones(f,1);

idf_cond = ones(f,1);
idf = find(idf_cond);
idn = find(~idf_cond);
assert(size(idn,1)+size(idf,1)==f);
if dim==2
    idfau = [idf;idf+f;idf+2*f];
    idnau = [idn;idn+f;idn+2*f];
else
    assert(dim==3);
    idfau = [idf;idf+f;idf+2*f;idf+3*f;idf+4*f;idf+5*f;];
    idnau = [idn;idn+f;idn+2*f;idn+3*f;idn+4*f;idn+5*f;];
end


syms aA1(p1,p2,p3);
syms aA2(p1,p2,p3);
syms aA3(p1,p2,p3);



switch para_type
    
    case 'LU'

    delta = 0;
    cond_bound = 0;
    aA1(p1,p2,p3) = p1.^2 + delta + cond_bound * (p2.^2 + p3.^2);
    aA2(p1,p2,p3) = p1.*p2;
    aA3(p1,p2,p3) = p2.^2 + p3.^2 + delta + cond_bound * p1.^2;

    fprintf("LU parameterization of tensor field with delta=%g, cond_bound=%g.\n",delta,cond_bound);

    case 'diag'

    aA1(p1,p2,p3) = p1;
    aA2(p1,p2,p3) = 0;
    aA3(p1,p2,p3) = p3;

    case 'invdiag'

    aA1(p1,p2,p3) = 1/p1;
    aA2(p1,p2,p3) = 0;
    aA3(p1,p2,p3) = 1/p3;
    
    otherwise 
        error('Unsupported para type~\n');
        
end

mA1 = matlabFunction(aA1);
mA2 = matlabFunction(aA2);
mA3 = matlabFunction(aA3);

% sp_mat0 = sparse(idfau,1:length(idfau),ones(length(idfau),1),3*f,length(idfau));
% sp_mat1 = sparse(idnau,1:length(idnau),ones(length(idnau),1),3*f,length(idnau));

s_at2au = @(aatt)  ...
    [(mA1(aatt(:,1),aatt(:,2),aatt(:,3))).*FA(idf).*base(idf);...
     (mA2(aatt(:,1),aatt(:,2),aatt(:,3))).*FA(idf).*base(idf);...
     (mA3(aatt(:,1),aatt(:,2),aatt(:,3))).*FA(idf).*base(idf);]...
    ;

s_at2au_mul_core = @(aatt) ...
    [(mA1(aatt(:,1),aatt(:,2),aatt(:,3))).*FA(idf).*base(idf);...
     (mA2(aatt(:,1),aatt(:,2),aatt(:,3))).*FA(idf).*base(idf);...
     (mA3(aatt(:,1),aatt(:,2),aatt(:,3))).*FA(idf).*base(idf);]...
     ;

s_at2au_mul = @(aatt, g) s_at2au_mul_core(aatt); 
    
%%
    s_dAdP11 = matlabFunction(diff(aA1,p1));
    s_dAdP12 = matlabFunction(diff(aA1,p2));
    s_dAdP13 = matlabFunction(diff(aA1,p3));

    s_dAdP21 = matlabFunction(diff(aA2,p1));
    s_dAdP22 = matlabFunction(diff(aA2,p2));
    s_dAdP23 = matlabFunction(diff(aA2,p3));

    s_dAdP31 = matlabFunction(diff(aA3,p1));
    s_dAdP32 = matlabFunction(diff(aA3,p2));
    s_dAdP33 = matlabFunction(diff(aA3,p3));
    
    % it is critical to use matlabFunction to convert the symbolic function to
    % a function handle. 

    s_11 = @(aatt) (s_dAdP11(aatt(:,1),aatt(:,2),aatt(:,3)));
    s_21 = @(aatt) (s_dAdP21(aatt(:,1),aatt(:,2),aatt(:,3)));
    s_31 = @(aatt) (s_dAdP31(aatt(:,1),aatt(:,2),aatt(:,3)));

    s_12 = @(aatt) (s_dAdP12(aatt(:,1),aatt(:,2),aatt(:,3)));
    s_22 = @(aatt) (s_dAdP22(aatt(:,1),aatt(:,2),aatt(:,3)));
    s_32 = @(aatt) (s_dAdP32(aatt(:,1),aatt(:,2),aatt(:,3)));

    s_13 = @(aatt) (s_dAdP13(aatt(:,1),aatt(:,2),aatt(:,3)));
    s_23 = @(aatt) (s_dAdP23(aatt(:,1),aatt(:,2),aatt(:,3)));
    s_33 = @(aatt) (s_dAdP33(aatt(:,1),aatt(:,2),aatt(:,3)));


%     s_pdapdt_old = @(aatt) sparse(...
%                         [[idf+0*f;idf+1*f;idf+2*f];[idf+1*f;idf+2*f;idf+0*f];[idf+2*f;idf+0*f;idf+1*f];],...
%                         [1:length(idfau),1:length(idfau),1:length(idfau)],... % cannot use ';' must use ',' here!
%                         [[s_11(aatt).*FA(idf);s_22(aatt).*FA(idf);s_33(aatt).*FA(idf)];...
%                          [s_21(aatt).*FA(idf);s_32(aatt).*FA(idf);s_13(aatt).*FA(idf)];...
%                          [s_31(aatt).*FA(idf);s_12(aatt).*FA(idf);s_23(aatt).*FA(idf)];],...
%                         3*f,length(idfau));
    
    % s_pdapdt_old and s_pdapdt are equivalent: s_pdapdt_old(at) -
    % s_pdapdt(at)  yields sparse zeros matrix. 
    
    s_pdapdt = @(aatt) sparse(...
                        [[idf+0*f;idf+0*f;idf+0*f];[idf+1*f;idf+1*f;idf+1*f];[idf+2*f;idf+2*f;idf+2*f];],...
                        [1:length(idfau),1:length(idfau),1:length(idfau)],...
                        [[s_11(aatt).*FA(idf);s_12(aatt).*FA(idf);s_13(aatt).*FA(idf)];...
                         [s_21(aatt).*FA(idf);s_22(aatt).*FA(idf);s_23(aatt).*FA(idf)];...
                         [s_31(aatt).*FA(idf);s_32(aatt).*FA(idf);s_33(aatt).*FA(idf)];],...
                        3*f,length(idfau));
                    
                    
    sp_mat_T1 = sparse(...
                        [idf+0*f;idf+0*f;idf+0*f],...
                        1:length(idfau),...
                        [FA(idf);FA(idf);FA(idf)],...
                        3*f,length(idfau))';
    sp_mat_T2 = sparse(...
                        [idf+1*f;idf+1*f;idf+1*f],...
                        1:length(idfau),...
                        [FA(idf);FA(idf);FA(idf)],...
                        3*f,length(idfau))';                
    sp_mat_T3 = sparse(...
                        [idf+2*f;idf+2*f;idf+2*f],...
                        1:length(idfau),...
                        [FA(idf);FA(idf);FA(idf)],...
                        3*f,length(idfau))';
    
    s_pdapdt_lmul = @(aatt,g) ...
          (sp_mat_T1*g).*[s_11(aatt)+0.0*idf;s_12(aatt)+0.0*idf;s_13(aatt)+0.0*idf] ...
        + (sp_mat_T2*g).*[s_21(aatt)+0.0*idf;s_22(aatt)+0.0*idf;s_23(aatt)+0.0*idf] ...
        + (sp_mat_T3*g).*[s_31(aatt)+0.0*idf;s_32(aatt)+0.0*idf;s_33(aatt)+0.0*idf] ...
        ; % note +0.0*idf is the hack to maintain size when s_ij yields zero. 

end

