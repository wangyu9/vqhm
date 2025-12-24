function [tp] = tensor_para(FA,dim,para_type,varargin)

reg_coeff = 0; % 1;
% this is not yet used anyway. 

nvar = length(varargin);
ii=1;

inv_para_fun = false;

experimental = false; % this option is much slower. 

while(ii<=nvar)
   if(strcmp(varargin{ii},'RegCoeff'))
       reg_coeff = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'InversePara'))
        inv_para_fun = true;
   elseif(strcmp(varargin{ii},'Experimental'))
        experimental = true;
   else
       error('unknown parameter!\n');
   end
   ii = ii + 1;
end

%%
sparse_diag = @(ddd) sparse(1:numel(ddd), 1:numel(ddd), ddd, numel(ddd), numel(ddd));

tp = struct();

tp.para_type = para_type;
tp.test_mark = false;

tp.alpha = 1;
tp.beta = 0;
tp.gamma = 0;

tp.energy_uv = [];
% isempty(tp.energy_uv) is used to determine wheather to add tp.energy_uv

% example: 
% tp.energy_uv = @(uv) deal(...
%     0, ...
%     zeros(numel(uv),1) ); % energy value and grad. 

% regularizer in at
tp.reg_value = @(at) 0;
tp.reg_grad  = @(at) 0;

% regularizer in au
tp.rau_value = @(au) 0;
tp.rau_grad  = @(au) 0;

reg = 0;

tp.pre_cond = [];

tp.append_grad = zeros(0,1);

tp.conj_vmap = true;

    Area = FA;
    f = size(FA,1);
    
    assert(dim==2);

    if false
    
            if false
                % diag without mass matrix. 

                s_at2au = @(aatt) [aatt(:,1);zeros(f,1);aatt(:,3)];
                s_pdapdt_lmul = @(aatt,gg) [ones(f,1);zeros(f,1);ones(f,1)].* gg;

                s_at2au_v = @(aatt) [1./aatt(:,1);zeros(f,1);1./aatt(:,3)];
                s_pdapdt_lmul_v = @(aatt,gg)  [-1./aatt(:,1).^2;zeros(f,1);-1./aatt(:,3).^2].* gg;
            else
                % with mass matrix.
                % this is for diagonal tensor.

                s_at2au = @(aatt) [Area.*aatt(:,1);zeros(f,1);Area.*aatt(:,3)];
                s_pdapdt_lmul = @(aatt,gg) [Area;zeros(f,1);Area].* gg;

                s_at2au_v = @(aatt) [Area./aatt(:,1);zeros(f,1);Area./aatt(:,3)];
                s_pdapdt_lmul_v = @(aatt,gg)  [-Area./aatt(:,1).^2;zeros(f,1);-Area./aatt(:,3).^2].* gg; 
            end
    end



if experimental
        %%
%             syms aA1(p1,p2,p3);
%             syms aA2(p1,p2,p3);
%             syms aA3(p1,p2,p3);
%                   
%             syms bA1(p1,p2,p3);
%             syms bA2(p1,p2,p3);
%             syms bA3(p1,p2,p3);

        np = 3;
        
        switch para_type
        	case 'bbt'
                np = 4;
            case 'bbt-20'
                np = 20;
            otherwise 
                np = 3;
        end

        p = sym('p',[np,1],'real');
        p1 = p(1);
        p2 = p(2);
        p3 = p(3);
        
        p4 = 0;
        if np>=4
           p4 = p(4);
        end
        for ii=5:np
           eval(['p', int2str(ii),' = p(',  int2str(ii),')']);
        end

        if inv_para_fun
            z = sym('z',[3,1],'real');
            z1 = z(1);
            z2 = z(2);
            z3 = z(3);
        end
        
        if true


            switch para_type

                case 'diag-no-mass'
                case 'diag'

                    aA1 = p1;
                    aA2 = 0;
                    aA3 = p3;

                    bA1 = 1/p1;
                    bA2 = 0;
                    bA3 = 1/p3;

                case 'diag-sq'

                    aA1 = p1.^2;
                    aA2 = 0;
                    aA3 = p3.^2;

                    bA1 = 1/p1.^2;
                    bA2 = 0;
                    bA3 = 1/p3.^2;

                case 'complex-sym'

                    aA1 = -((p1 - 1)^2 - abs(p3)^2 + p2^2)/(p1^2 - abs(p3 + 1)^2 + p2^2);
                    aA2 = (2*p2)/(p1^2 - abs(p3 + 1)^2 + p2^2);
                    aA3 = -((p1 + 1)^2 - abs(p3)^2 + p2^2)/(p1^2 - abs(p3 + 1)^2 + p2^2);

                    bA1 = ((p1^2 - abs(p3 + 1)^2 + p2^2)*(2*p1 - abs(p3)^2 + p1^2 + p2^2 + 1))/(2*abs(p3)^2 - abs(p3)^4 + 2*p1^2*abs(p3)^2 + 2*p2^2*abs(p3)^2 + 2*p1^2 + 2*p2^2 - p1^4 - p2^4 - 2*p1^2*p2^2 - 1);
                    bA2 = (2*p2*(p1^2 - abs(p3 + 1)^2 + p2^2))/(2*abs(p3)^2 - abs(p3)^4 + 2*p1^2*abs(p3)^2 + 2*p2^2*abs(p3)^2 + 2*p1^2 + 2*p2^2 - p1^4 - p2^4 - 2*p1^2*p2^2 - 1);
                    bA3 = ((p1^2 - abs(p3 + 1)^2 + p2^2)*(p1^2 - abs(p3)^2 - 2*p1 + p2^2 + 1))/(2*abs(p3)^2 - abs(p3)^4 + 2*p1^2*abs(p3)^2 + 2*p2^2*abs(p3)^2 + 2*p1^2 + 2*p2^2 - p1^4 - p2^4 - 2*p1^2*p2^2 - 1);

                case 'complex-sym-det1'

                    aA1 = -((p1 - 1)^2 + p2^2)/(p1^2 - 1 + p2^2);
                    aA2 = (2*p2)/(p1^2 - 1 + p2^2);
                    aA3 = -((p1 + 1)^2 + p2^2)/(p1^2 - 1 + p2^2);

                    bA1 = aA3;
                    bA2 = -aA2;
                    bA3 = aA1;

                    % dett = [aA1,aA2;aA2,aA3] * [bA1,bA2;bA2,bA3]; 
                    % dett should be eye(2). 
                    % e.g. see if it returns eye(2): eval(dett(0.39,0.43,0.13))

                case 'complex-plane-det1-without-epsilon-cannot-work'
                    
%                     b1 = p1 / sqrt(p1^2+p2^2) * (1 - 1/(1+log(1+sqrt(p1^2+p2^2))));
%                     b2 = p2 / sqrt(p1^2+p2^2) * (1 - 1/(1+log(1+sqrt(p1^2+p2^2))));
%                     b3 = 0;
%                     
%                     -((b1 - 1)^2 + b2^2)/(b1^2 - 1 + b2^2);
%                     (2*b2)/(b1^2 - 1 + b2^2);
%                     -((b1 + 1)^2 + b2^2)/(b1^2 - 1 + b2^2);
                    
                    aA1 = -(((p1*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1))/(p1^2 + p2^2)^(1/2) + 1)^2 + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2))/((p1^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2) + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2) - 1);
                    aA2 = -(2*p2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1))/((p1^2 + p2^2)^(1/2)*((p1^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2) + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2) - 1));
                    aA3 = -(((p1*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1))/(p1^2 + p2^2)^(1/2) - 1)^2 + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2))/((p1^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2) + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2) - 1);
                    
                    bA1 = aA3;
                    bA2 = -aA2;
                    bA3 = aA1;
                    
                case 'complex-plane-det1-tanh'
                    
                    epsilon = 0e-10;
                    
                    b1 = p1 / sqrt(p1^2+p2^2+epsilon) *  tanh(sqrt(p1^2+p2^2)) ;
                    b2 = p2 / sqrt(p1^2+p2^2+epsilon) *  tanh(sqrt(p1^2+p2^2)) ;
                    b3 = 0;
                    
                    %b1 = p1 / sqrt(p1^2+p2^2+epsilon) *  sqrt(p1^2+p2^2) / (1+sqrt(p1^2+p2^2)) ;
                    %b2 = p2 / sqrt(p1^2+p2^2+epsilon) *  sqrt(p1^2+p2^2) / (1+sqrt(p1^2+p2^2)) ;
                    %b3 = 0;
                    
                    -((b1 - 1)^2 + b2^2)/(b1^2 - 1 + b2^2);
                    (2*b2)/(b1^2 - 1 + b2^2);
                    -((b1 + 1)^2 + b2^2)/(b1^2 - 1 + b2^2);

%                     aA1 = -(((p1*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1))/(p1^2 + p2^2 + 1/1000000000000000000)^(1/2) + 1)^2 + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000))/((p1^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000) + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000) - 1);
%                     aA2 = -(2*p2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1))/(((p1^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000) + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000) - 1)*(p1^2 + p2^2 + 1/1000000000000000000)^(1/2));
%                     aA3 = -(((p1*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1))/(p1^2 + p2^2 + 1/1000000000000000000)^(1/2) - 1)^2 + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000))/((p1^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000) + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000) - 1);
%                     
                    aA1 = -((b1 - 1)^2 + b2^2)/(b1^2 - 1 + b2^2);
                    aA2 = (2*b2)/(b1^2 - 1 + b2^2);
                    aA3 = -((b1 + 1)^2 + b2^2)/(b1^2 - 1 + b2^2);
                    
                    bA1 = aA3;
                    bA2 = -aA2;
                    bA3 = aA1;
                    % to double check results using:
                    %  eval(subs(det([aA1,aA2;aA2,aA3]), [p1,p2], [-5,-20] ))
                    
                    reg = reg_coeff * 0.5 *(b1^2 + b2^2);
                  
                    
                case 'complex-plane-det1'
                    
                    epsilon = 0e-10;
                    
                    b1 = p1 / sqrt(p1^2+p2^2+epsilon) * (1 - 1/(1+log(1+sqrt(p1^2+p2^2))));
                    b2 = p2 / sqrt(p1^2+p2^2+epsilon) * (1 - 1/(1+log(1+sqrt(p1^2+p2^2))));
                    b3 = 0;
                    
                    -((b1 - 1)^2 + b2^2)/(b1^2 - 1 + b2^2);
                    (2*b2)/(b1^2 - 1 + b2^2);
                    -((b1 + 1)^2 + b2^2)/(b1^2 - 1 + b2^2);

%                     aA1 = -(((p1*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1))/(p1^2 + p2^2 + 1/1000000000000000000)^(1/2) + 1)^2 + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000))/((p1^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000) + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000) - 1);
%                     aA2 = -(2*p2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1))/(((p1^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000) + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000) - 1)*(p1^2 + p2^2 + 1/1000000000000000000)^(1/2));
%                     aA3 = -(((p1*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1))/(p1^2 + p2^2 + 1/1000000000000000000)^(1/2) - 1)^2 + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000))/((p1^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000) + (p2^2*(1/(log((p1^2 + p2^2)^(1/2) + 1) + 1) - 1)^2)/(p1^2 + p2^2 + 1/1000000000000000000) - 1);
%                     
                    aA1 = -((b1 - 1)^2 + b2^2)/(b1^2 - 1 + b2^2);
                    aA2 = (2*b2)/(b1^2 - 1 + b2^2);
                    aA3 = -((b1 + 1)^2 + b2^2)/(b1^2 - 1 + b2^2);
                    
                    bA1 = aA3;
                    bA2 = -aA2;
                    bA3 = aA1;
                    % to double check results using:
                    %  eval(subs(det([aA1,aA2;aA2,aA3]), [p1,p2], [-5,-20] ))
                    
                    reg = reg_coeff * 0.5 *(b1^2 + b2^2);
                    
                    if inv_para_fun
                        
                    end
                    
                case 'complex-plane-det1-tmp'
                    
                    epsilon = 0;
                    
                    b1 = p1 / sqrt(p1^2+p2^2+epsilon) * (sqrt(p1^2+p2^2)/(1+p1^2+p2^2));
                    b2 = p2 / sqrt(p1^2+p2^2+epsilon) * (sqrt(p1^2+p2^2)/(1+p1^2+p2^2));
                    b3 = 0;
                    
                    aA1 = -((b1 - 1)^2 + b2^2)/(b1^2 - 1 + b2^2);
                    aA2 = (2*b2)/(b1^2 - 1 + b2^2);
                    aA3 = -((b1 + 1)^2 + b2^2)/(b1^2 - 1 + b2^2);
                    
                    bA1 = aA3;
                    bA2 = -aA2;
                    bA3 = aA1;
                    
                case 'both-inverse'

                    aA1 = p3/(- p2^2 + p1*p3);
                    aA2 = -p2/(- p2^2 + p1*p3);
                    aA3 = p1/(- p2^2 + p1*p3);

                    bA1 = p3/(- p2^2 + p1*p3);
                    bA2 = -p2/(- p2^2 + p1*p3);
                    bA3 = p1/(- p2^2 + p1*p3);

                case 'llt'

                    aA1 = p1.^2;
                    aA2 = p1.*p2;
                    aA3 = p2.^2 + p3.^2;

                    bA1 = (p2^2 + p3^2)/(p1^2*p3^2);
                    bA2 = -p2/(p1*p3^2);
                    bA3 = 1/p3^2;
                    
                case 'aa'
                    
                    % S = [p1, p2; p2, p3]
                    % S*S
                    
                    aA1 = p1^2 + p2^2;
                    aA2 = p1*p2 + p2*p3;
                    aA3 = p2^2 + p3^2;

                    bA1 = (p2^2 + p3^2)/(p1^2*p3^2 - 2*p1*p2^2*p3 + p2^4);
                    bA2 = -(p1*p2 + p2*p3)/(p1^2*p3^2 - 2*p1*p2^2*p3 + p2^4);
                    bA3 = (p1^2 + p2^2)/(p1^2*p3^2 - 2*p1*p2^2*p3 + p2^4);
                    
                case 'bbt'
                    
                    % B = [p1, p2; p3, p4]
                    % A = B*B';
                    
                    aA1 = p1^2 + p2^2;
                    aA2 = p1*p3 + p2*p4;
                    aA3 = p3^2 + p4^2;
                    
                    bA1 = (p3^2 + p4^2)/(p1^2*p4^2 - 2*p1*p2*p3*p4 + p2^2*p3^2);
                    bA2 = -(p1*p3 + p2*p4)/(p1^2*p4^2 - 2*p1*p2*p3*p4 + p2^2*p3^2);
                    bA3 = (p1^2 + p2^2)/(p1^2*p4^2 - 2*p1*p2*p3*p4 + p2^2*p3^2);
                    
                case 'bbt-20'
                    
                    % B = reshape(p,[10,2])';
                    % A = B*B';
                    
                    aA1 = p1^2 + p2^2 + p3^2 + p4^2 + p5^2 + p6^2 + p7^2 + p8^2 + p9^2 + p10^2;
                    aA2 = p1*p11 + p2*p12 + p3*p13 + p4*p14 + p5*p15 + p6*p16 + p7*p17 + p8*p18 + p9*p19 + p10*p20;
                    aA3 = p11^2 + p12^2 + p13^2 + p14^2 + p15^2 + p16^2 + p17^2 + p18^2 + p19^2 + p20^2;
                    
                    bA1 = (p11^2 + p12^2 + p13^2 + p14^2 + p15^2 + p16^2 + p17^2 + p18^2 + p19^2 + p20^2)/(p1^2*p12^2 + p1^2*p13^2 + p1^2*p14^2 + p1^2*p15^2 + p1^2*p16^2 + p1^2*p17^2 + p1^2*p18^2 + p1^2*p19^2 + p1^2*p20^2 - 2*p1*p2*p11*p12 - 2*p1*p3*p11*p13 - 2*p1*p4*p11*p14 - 2*p1*p5*p11*p15 - 2*p1*p6*p11*p16 - 2*p1*p7*p11*p17 - 2*p1*p8*p11*p18 - 2*p1*p9*p11*p19 - 2*p1*p10*p11*p20 + p2^2*p11^2 + p2^2*p13^2 + p2^2*p14^2 + p2^2*p15^2 + p2^2*p16^2 + p2^2*p17^2 + p2^2*p18^2 + p2^2*p19^2 + p2^2*p20^2 - 2*p2*p3*p12*p13 - 2*p2*p4*p12*p14 - 2*p2*p5*p12*p15 - 2*p2*p6*p12*p16 - 2*p2*p7*p12*p17 - 2*p2*p8*p12*p18 - 2*p2*p9*p12*p19 - 2*p2*p10*p12*p20 + p3^2*p11^2 + p3^2*p12^2 + p3^2*p14^2 + p3^2*p15^2 + p3^2*p16^2 + p3^2*p17^2 + p3^2*p18^2 + p3^2*p19^2 + p3^2*p20^2 - 2*p3*p4*p13*p14 - 2*p3*p5*p13*p15 - 2*p3*p6*p13*p16 - 2*p3*p7*p13*p17 - 2*p3*p8*p13*p18 - 2*p3*p9*p13*p19 - 2*p3*p10*p13*p20 + p4^2*p11^2 + p4^2*p12^2 + p4^2*p13^2 + p4^2*p15^2 + p4^2*p16^2 + p4^2*p17^2 + p4^2*p18^2 + p4^2*p19^2 + p4^2*p20^2 - 2*p4*p5*p14*p15 - 2*p4*p6*p14*p16 - 2*p4*p7*p14*p17 - 2*p4*p8*p14*p18 - 2*p4*p9*p14*p19 - 2*p4*p10*p14*p20 + p5^2*p11^2 + p5^2*p12^2 + p5^2*p13^2 + p5^2*p14^2 + p5^2*p16^2 + p5^2*p17^2 + p5^2*p18^2 + p5^2*p19^2 + p5^2*p20^2 - 2*p5*p6*p15*p16 - 2*p5*p7*p15*p17 - 2*p5*p8*p15*p18 - 2*p5*p9*p15*p19 - 2*p5*p10*p15*p20 + p6^2*p11^2 + p6^2*p12^2 + p6^2*p13^2 + p6^2*p14^2 + p6^2*p15^2 + p6^2*p17^2 + p6^2*p18^2 + p6^2*p19^2 + p6^2*p20^2 - 2*p6*p7*p16*p17 - 2*p6*p8*p16*p18 - 2*p6*p9*p16*p19 - 2*p6*p10*p16*p20 + p7^2*p11^2 + p7^2*p12^2 + p7^2*p13^2 + p7^2*p14^2 + p7^2*p15^2 + p7^2*p16^2 + p7^2*p18^2 + p7^2*p19^2 + p7^2*p20^2 - 2*p7*p8*p17*p18 - 2*p7*p9*p17*p19 - 2*p7*p10*p17*p20 + p8^2*p11^2 + p8^2*p12^2 + p8^2*p13^2 + p8^2*p14^2 + p8^2*p15^2 + p8^2*p16^2 + p8^2*p17^2 + p8^2*p19^2 + p8^2*p20^2 - 2*p8*p9*p18*p19 - 2*p8*p10*p18*p20 + p9^2*p11^2 + p9^2*p12^2 + p9^2*p13^2 + p9^2*p14^2 + p9^2*p15^2 + p9^2*p16^2 + p9^2*p17^2 + p9^2*p18^2 + p9^2*p20^2 - 2*p9*p10*p19*p20 + p10^2*p11^2 + p10^2*p12^2 + p10^2*p13^2 + p10^2*p14^2 + p10^2*p15^2 + p10^2*p16^2 + p10^2*p17^2 + p10^2*p18^2 + p10^2*p19^2);
                    bA2 = -(p1*p11 + p2*p12 + p3*p13 + p4*p14 + p5*p15 + p6*p16 + p7*p17 + p8*p18 + p9*p19 + p10*p20)/(p1^2*p12^2 + p1^2*p13^2 + p1^2*p14^2 + p1^2*p15^2 + p1^2*p16^2 + p1^2*p17^2 + p1^2*p18^2 + p1^2*p19^2 + p1^2*p20^2 - 2*p1*p2*p11*p12 - 2*p1*p3*p11*p13 - 2*p1*p4*p11*p14 - 2*p1*p5*p11*p15 - 2*p1*p6*p11*p16 - 2*p1*p7*p11*p17 - 2*p1*p8*p11*p18 - 2*p1*p9*p11*p19 - 2*p1*p10*p11*p20 + p2^2*p11^2 + p2^2*p13^2 + p2^2*p14^2 + p2^2*p15^2 + p2^2*p16^2 + p2^2*p17^2 + p2^2*p18^2 + p2^2*p19^2 + p2^2*p20^2 - 2*p2*p3*p12*p13 - 2*p2*p4*p12*p14 - 2*p2*p5*p12*p15 - 2*p2*p6*p12*p16 - 2*p2*p7*p12*p17 - 2*p2*p8*p12*p18 - 2*p2*p9*p12*p19 - 2*p2*p10*p12*p20 + p3^2*p11^2 + p3^2*p12^2 + p3^2*p14^2 + p3^2*p15^2 + p3^2*p16^2 + p3^2*p17^2 + p3^2*p18^2 + p3^2*p19^2 + p3^2*p20^2 - 2*p3*p4*p13*p14 - 2*p3*p5*p13*p15 - 2*p3*p6*p13*p16 - 2*p3*p7*p13*p17 - 2*p3*p8*p13*p18 - 2*p3*p9*p13*p19 - 2*p3*p10*p13*p20 + p4^2*p11^2 + p4^2*p12^2 + p4^2*p13^2 + p4^2*p15^2 + p4^2*p16^2 + p4^2*p17^2 + p4^2*p18^2 + p4^2*p19^2 + p4^2*p20^2 - 2*p4*p5*p14*p15 - 2*p4*p6*p14*p16 - 2*p4*p7*p14*p17 - 2*p4*p8*p14*p18 - 2*p4*p9*p14*p19 - 2*p4*p10*p14*p20 + p5^2*p11^2 + p5^2*p12^2 + p5^2*p13^2 + p5^2*p14^2 + p5^2*p16^2 + p5^2*p17^2 + p5^2*p18^2 + p5^2*p19^2 + p5^2*p20^2 - 2*p5*p6*p15*p16 - 2*p5*p7*p15*p17 - 2*p5*p8*p15*p18 - 2*p5*p9*p15*p19 - 2*p5*p10*p15*p20 + p6^2*p11^2 + p6^2*p12^2 + p6^2*p13^2 + p6^2*p14^2 + p6^2*p15^2 + p6^2*p17^2 + p6^2*p18^2 + p6^2*p19^2 + p6^2*p20^2 - 2*p6*p7*p16*p17 - 2*p6*p8*p16*p18 - 2*p6*p9*p16*p19 - 2*p6*p10*p16*p20 + p7^2*p11^2 + p7^2*p12^2 + p7^2*p13^2 + p7^2*p14^2 + p7^2*p15^2 + p7^2*p16^2 + p7^2*p18^2 + p7^2*p19^2 + p7^2*p20^2 - 2*p7*p8*p17*p18 - 2*p7*p9*p17*p19 - 2*p7*p10*p17*p20 + p8^2*p11^2 + p8^2*p12^2 + p8^2*p13^2 + p8^2*p14^2 + p8^2*p15^2 + p8^2*p16^2 + p8^2*p17^2 + p8^2*p19^2 + p8^2*p20^2 - 2*p8*p9*p18*p19 - 2*p8*p10*p18*p20 + p9^2*p11^2 + p9^2*p12^2 + p9^2*p13^2 + p9^2*p14^2 + p9^2*p15^2 + p9^2*p16^2 + p9^2*p17^2 + p9^2*p18^2 + p9^2*p20^2 - 2*p9*p10*p19*p20 + p10^2*p11^2 + p10^2*p12^2 + p10^2*p13^2 + p10^2*p14^2 + p10^2*p15^2 + p10^2*p16^2 + p10^2*p17^2 + p10^2*p18^2 + p10^2*p19^2);
                    bA3 = (p1^2 + p2^2 + p3^2 + p4^2 + p5^2 + p6^2 + p7^2 + p8^2 + p9^2 + p10^2)/(p1^2*p12^2 + p1^2*p13^2 + p1^2*p14^2 + p1^2*p15^2 + p1^2*p16^2 + p1^2*p17^2 + p1^2*p18^2 + p1^2*p19^2 + p1^2*p20^2 - 2*p1*p2*p11*p12 - 2*p1*p3*p11*p13 - 2*p1*p4*p11*p14 - 2*p1*p5*p11*p15 - 2*p1*p6*p11*p16 - 2*p1*p7*p11*p17 - 2*p1*p8*p11*p18 - 2*p1*p9*p11*p19 - 2*p1*p10*p11*p20 + p2^2*p11^2 + p2^2*p13^2 + p2^2*p14^2 + p2^2*p15^2 + p2^2*p16^2 + p2^2*p17^2 + p2^2*p18^2 + p2^2*p19^2 + p2^2*p20^2 - 2*p2*p3*p12*p13 - 2*p2*p4*p12*p14 - 2*p2*p5*p12*p15 - 2*p2*p6*p12*p16 - 2*p2*p7*p12*p17 - 2*p2*p8*p12*p18 - 2*p2*p9*p12*p19 - 2*p2*p10*p12*p20 + p3^2*p11^2 + p3^2*p12^2 + p3^2*p14^2 + p3^2*p15^2 + p3^2*p16^2 + p3^2*p17^2 + p3^2*p18^2 + p3^2*p19^2 + p3^2*p20^2 - 2*p3*p4*p13*p14 - 2*p3*p5*p13*p15 - 2*p3*p6*p13*p16 - 2*p3*p7*p13*p17 - 2*p3*p8*p13*p18 - 2*p3*p9*p13*p19 - 2*p3*p10*p13*p20 + p4^2*p11^2 + p4^2*p12^2 + p4^2*p13^2 + p4^2*p15^2 + p4^2*p16^2 + p4^2*p17^2 + p4^2*p18^2 + p4^2*p19^2 + p4^2*p20^2 - 2*p4*p5*p14*p15 - 2*p4*p6*p14*p16 - 2*p4*p7*p14*p17 - 2*p4*p8*p14*p18 - 2*p4*p9*p14*p19 - 2*p4*p10*p14*p20 + p5^2*p11^2 + p5^2*p12^2 + p5^2*p13^2 + p5^2*p14^2 + p5^2*p16^2 + p5^2*p17^2 + p5^2*p18^2 + p5^2*p19^2 + p5^2*p20^2 - 2*p5*p6*p15*p16 - 2*p5*p7*p15*p17 - 2*p5*p8*p15*p18 - 2*p5*p9*p15*p19 - 2*p5*p10*p15*p20 + p6^2*p11^2 + p6^2*p12^2 + p6^2*p13^2 + p6^2*p14^2 + p6^2*p15^2 + p6^2*p17^2 + p6^2*p18^2 + p6^2*p19^2 + p6^2*p20^2 - 2*p6*p7*p16*p17 - 2*p6*p8*p16*p18 - 2*p6*p9*p16*p19 - 2*p6*p10*p16*p20 + p7^2*p11^2 + p7^2*p12^2 + p7^2*p13^2 + p7^2*p14^2 + p7^2*p15^2 + p7^2*p16^2 + p7^2*p18^2 + p7^2*p19^2 + p7^2*p20^2 - 2*p7*p8*p17*p18 - 2*p7*p9*p17*p19 - 2*p7*p10*p17*p20 + p8^2*p11^2 + p8^2*p12^2 + p8^2*p13^2 + p8^2*p14^2 + p8^2*p15^2 + p8^2*p16^2 + p8^2*p17^2 + p8^2*p19^2 + p8^2*p20^2 - 2*p8*p9*p18*p19 - 2*p8*p10*p18*p20 + p9^2*p11^2 + p9^2*p12^2 + p9^2*p13^2 + p9^2*p14^2 + p9^2*p15^2 + p9^2*p16^2 + p9^2*p17^2 + p9^2*p18^2 + p9^2*p20^2 - 2*p9*p10*p19*p20 + p10^2*p11^2 + p10^2*p12^2 + p10^2*p13^2 + p10^2*p14^2 + p10^2*p15^2 + p10^2*p16^2 + p10^2*p17^2 + p10^2*p18^2 + p10^2*p19^2);

                                        
                case 'weber'

                    aA1 = sqrt(p1) + 2 * p2;
                    aA2 = 2 * p3;
                    aA3 = sqrt(p1) - 2 * p2;

                    bA1 = (2*p2 - p1^(1/2))/(4*p2^2 + 4*p3^2 - p1);
                    bA2 = (2*p3)/(4*p2^2 + 4*p3^2 - p1);
                    bA3 = -(2*p2 + p1^(1/2))/(4*p2^2 + 4*p3^2 - p1);

                case 'direct'

                    aA1 = p1;
                    aA2 = p2;
                    aA3 = p3;

                    bA1 = p3/(- p2^2 + p1*p3);
                    bA2 = -p2/(- p2^2 + p1*p3);
                    bA3 = p1/(- p2^2 + p1*p3);

                otherwise 
                    error('Unsupported para type~\n');

            end

        else
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

        % complex-sym-det1:

        aA1(p1,p2,p3) = -((p1 - 1)^2 + p2^2)/(p1^2 - 1 + p2^2);
        aA2(p1,p2,p3) = (2*p2)/(p1^2 - 1 + p2^2);
        aA3(p1,p2,p3) = -((p1 + 1)^2 + p2^2)/(p1^2 - 1 + p2^2);

        bA1(p1,p2,p3) = aA3;
        bA2(p1,p2,p3) = -aA2;
        bA3(p1,p2,p3) = aA1;

        % dett = [aA1,aA2;aA2,aA3] * [bA1,bA2;bA2,bA3]; 
        % dett should be eye(2). 
        % e.g. see if it returns eye(2): eval(dett(0.39,0.43,0.13))

        % weber

        aA1(p1,p2,p3) = sqrt(p1) + 2 * p2;
        aA2(p1,p2,p3) = 2 * p3;
        aA3(p1,p2,p3) = sqrt(p1) - 2 * p2;

        bA1(p1,p2,p3) = (2*p2 - p1^(1/2))/(4*p2^2 + 4*p3^2 - p1);
        bA2(p1,p2,p3) = (2*p3)/(4*p2^2 + 4*p3^2 - p1);
        bA3(p1,p2,p3) = -(2*p2 + p1^(1/2))/(4*p2^2 + 4*p3^2 - p1);

        % direct

        aA1(p1,p2,p3) = p1;
        aA2(p1,p2,p3) = p2;
        aA3(p1,p2,p3) = p3;

        bA1(p1,p2,p3) = p3/(- p2^2 + p1*p3);
        bA2(p1,p2,p3) = -p2/(- p2^2 + p1*p3);
        bA3(p1,p2,p3) = p1/(- p2^2 + p1*p3);

        % both-inverse

        aA1(p1,p2,p3) = p3/(- p2^2 + p1*p3);
        aA2(p1,p2,p3) = -p2/(- p2^2 + p1*p3);
        aA3(p1,p2,p3) = p1/(- p2^2 + p1*p3);

        bA1(p1,p2,p3) = p3/(- p2^2 + p1*p3);
        bA2(p1,p2,p3) = -p2/(- p2^2 + p1*p3);
        bA3(p1,p2,p3) = p1/(- p2^2 + p1*p3);

        % LLT:

        aA1(p1,p2,p3) = p1.^2;
        aA2(p1,p2,p3) = p1.*p2;
        aA3(p1,p2,p3) = p2.^2 + p3.^2;

        bA1(p1,p2,p3) = (p2^2 + p3^2)/(p1^2*p3^2);
        bA2(p1,p2,p3) = -p2/(p1*p3^2);
        bA3(p1,p2,p3) = 1/p3^2;

        end
        
        % varlist = {p1,p2,p3};
        % varlist = cell(1,np);
%         varlist = {};
%         for ii=1:np
%            varlist(ii) = p(ii); 
%         end
        
        fun_handle = @(xxx) matlabFunction((xxx)+0*sum(p), 'Vars', p);
        % adding 0*sum(p) as a workaround for 0 for avoid the undefine error 
        % for input arguments of type 'cell'.
        
        mR = fun_handle(reg);

        mA1 = fun_handle(aA1);
        mA2 = fun_handle(aA2);
        mA3 = fun_handle(aA3);
        
        s_dAdP1 = cell(1,np);
        s_dAdP2 = cell(1,np);
        s_dAdP3 = cell(1,np);
    
        for ii=1:np
            s_dAdP1{ii} = fun_handle(diff(aA1,p(ii)));
            s_dAdP2{ii} = fun_handle(diff(aA2,p(ii)));
            s_dAdP3{ii} = fun_handle(diff(aA3,p(ii)));
        end
        
        disp('mA1=...'),disp(mA1);
        disp('mA2=...'),disp(mA2);
        disp('mA3=...'),disp(mA3);
        
        for ii=1:np
            fprintf('s_dAdP1{%d}=\n',ii),disp(s_dAdP1{ii});
        end
        for ii=1:np
            fprintf('s_dAdP2{%d}=\n',ii),disp(s_dAdP2{ii});
        end
        for ii=1:np
            fprintf('s_dAdP3{%d}=\n',ii),disp(s_dAdP3{ii});
        end

        % this hessian only support np==3

        hess_dA1_dP11 = fun_handle(diff(diff(aA1,p1),p1));
        hess_dA1_dP12 = fun_handle(diff(diff(aA1,p1),p2));
        hess_dA1_dP22 = fun_handle(diff(diff(aA1,p2),p2));
        
        hess_dA2_dP11 = fun_handle(diff(diff(aA2,p1),p1));
        hess_dA2_dP12 = fun_handle(diff(diff(aA2,p1),p2));
        hess_dA2_dP22 = fun_handle(diff(diff(aA2,p2),p2));
        
        hess_dA3_dP11 = fun_handle(diff(diff(aA3,p1),p1));
        hess_dA3_dP12 = fun_handle(diff(diff(aA3,p1),p2));
        hess_dA3_dP22 = fun_handle(diff(diff(aA3,p2),p2));

        
        %
        s_dRdP = cell(1,np);
        
        for ii=1:np
            s_dRdP{ii} = fun_handle(diff(reg,p(ii)));
        end
        
        disp('reg=...'),disp(mR);
        
        for ii=1:np
            fprintf('s_dRdP{%d}=\n',ii),disp(s_dRdP{ii});
        end
        %
        

        nA1 = fun_handle(bA1);
        nA2 = fun_handle(bA2);
        nA3 = fun_handle(bA3);
        
        t_dAdP1 = cell(1,3);
        t_dAdP2 = cell(1,3);
        t_dAdP3 = cell(1,3);

        for ii=1:np
            t_dAdP1{ii} = fun_handle(diff(bA1,p(ii)));
            t_dAdP2{ii} = fun_handle(diff(bA2,p(ii)));
            t_dAdP3{ii} = fun_handle(diff(bA3,p(ii)));
        end
        
        disp('nA1=...'),disp(nA1);
        disp('nA2=...'),disp(nA2);
        disp('nA3=...'),disp(nA3);

        for ii=1:np
            fprintf('t_dAdP1{%d}=\n',ii),disp(t_dAdP1{ii});
        end
        for ii=1:np
            fprintf('t_dAdP2{%d}=\n',ii),disp(t_dAdP2{ii});
        end
        for ii=1:np
            fprintf('t_dAdP3{%d}=\n',ii),disp(t_dAdP3{ii});
        end
else
    
    switch para_type

        case 'diag-no-mass'
        case 'diag'
            res_sym_grad_diag;
        case 'diag-sq'
            res_sym_grad_diag_sq;
        case 'complex-sym'
            res_sym_grad_complex_sym;            
        case 'complex-sym-det1'
            res_sym_grad_complex_sym_det1;
        case 'complex-plane-det1'
            ff = @(rr) (1 - 1./(1+log(1+rr)));
            gg = @(rr) 1./(1+log(1+rr)).^2 ./ (1+rr);
            res_sym_grad_complex_plane_det1;
        case 'complex-plane-det1-tanh'
            ff = @(rr) tanh(rr);
            gg = @(rr) 1 - tanh(rr).^2;
            res_sym_grad_complex_plane_det1;
        case 'both-inverse'
            res_sym_grad_both_inverse;
        case 'llt'
            res_sym_grad_llt;
        case 'weber'
            res_sym_grad_weber;
        case 'direct'
            res_sym_grad_direct;
        otherwise 
            error('Unsupported para type~\n');

    end

end

%%
r = cell(1,np);

for ii=1:np
   r{ii} = @(aatt) apf( s_dRdP{ii}, aatt);
end

%%
s_1 = cell(1,np);
s_2 = cell(1,np);
s_3 = cell(1,np);

for ii=1:np
    s_1{ii} = @(aatt) apf( s_dAdP1{ii}, aatt);
    s_2{ii} = @(aatt) apf( s_dAdP2{ii}, aatt);
    s_3{ii} = @(aatt) apf( s_dAdP3{ii}, aatt);
end

t_1 = cell(1,np);
t_2 = cell(1,np);
t_3 = cell(1,np);

for ii=1:np
    t_1{ii} = @(aatt) apf( t_dAdP1{ii}, aatt);
    t_2{ii} = @(aatt) apf( t_dAdP2{ii}, aatt);
    t_3{ii} = @(aatt) apf( t_dAdP3{ii}, aatt);
end

s_at2reg = @(aatt) sum(Area .* apf( mR, aatt));

s_pdrpdt = @(aatt) assemble_pdrpdt(r, aatt, Area);


tp.s_at2reg = s_at2reg;
tp.s_pdrpdt = s_pdrpdt;


s_at2au = @(aatt) [...
     Area .* apf( mA1, aatt);...
     Area .* apf( mA2, aatt);...
     Area .* apf( mA3, aatt);...
    ];

s_pdapdt_lmul = @(aatt,gg) assemble_pdapdt(s_1, s_2, s_3, aatt, gg, Area);


% this one asserts two parameter per triangle. 
tp.s_pdapdt = @(aatt) [...
    sparse_diag(s_1{1}(aatt).*Area) , sparse_diag(s_1{2}(aatt).*Area);...
    sparse_diag(s_2{1}(aatt).*Area) , sparse_diag(s_2{2}(aatt).*Area);...
    sparse_diag(s_3{1}(aatt).*Area) , sparse_diag(s_3{2}(aatt).*Area);...
    ];


% this is for diagonal tensor:
% s_pdapdt_lmul = @(aatt,gg) [Area;zeros(f,1);Area].* gg;

% s_at2au_v = @(aatt) [Area./aatt(:,1);zeros(f,1);Area./aatt(:,3)];
% s_pdapdt_lmul_v = @(aatt,gg)  [-Area./aatt(:,1).^2;zeros(f,1);-Area./aatt(:,3).^2].* gg;

s_at2au_v = @(aatt) [...
     Area .* apf( nA1, aatt);...
     Area .* apf( nA2, aatt);...
     Area .* apf( nA3, aatt);...
    ];

s_pdapdt_lmul_v = @(aatt,gg) assemble_pdapdt(t_1, t_2, t_3, aatt, gg, Area);

% s_pdapdt_lmul_v = @(aatt,gg) [...
%     Area.*( t_1{1}(aatt).*gg(1:f) + t_2{1}(aatt).*gg(f+1:2*f) + t_3{1}(aatt).*gg(2*f+1:3*f) );...
%     Area.*( t_1{2}(aatt).*gg(1:f) + t_2{2}(aatt).*gg(f+1:2*f) + t_3{2}(aatt).*gg(2*f+1:3*f) );...
%     Area.*( t_1{3}(aatt).*gg(1:f) + t_2{3}(aatt).*gg(f+1:2*f) + t_3{3}(aatt).*gg(2*f+1:3*f) );...
%     ];

if experimental

tp.h_dA1_dP11 = @(aatt) Area .* apf( hess_dA1_dP11, aatt); 
tp.h_dA1_dP12 = @(aatt) Area .* apf( hess_dA1_dP12, aatt); 
tp.h_dA1_dP22 = @(aatt) Area .* apf( hess_dA1_dP22, aatt); 

tp.h_dA2_dP11 = @(aatt) Area .* apf( hess_dA2_dP11, aatt); 
tp.h_dA2_dP12 = @(aatt) Area .* apf( hess_dA2_dP12, aatt); 
tp.h_dA2_dP22 = @(aatt) Area .* apf( hess_dA2_dP22, aatt); 

tp.h_dA3_dP11 = @(aatt) Area .* apf( hess_dA3_dP11, aatt); 
tp.h_dA3_dP12 = @(aatt) Area .* apf( hess_dA3_dP12, aatt); 
tp.h_dA3_dP22 = @(aatt) Area .* apf( hess_dA3_dP22, aatt); 

end


tp.s_at2au = s_at2au;
tp.s_pdapdt_lmul = s_pdapdt_lmul;
            
tp.s_at2au_v = s_at2au_v;
tp.s_pdapdt_lmul_v = s_pdapdt_lmul_v;

end 

function [R] = apf(fun, aatt)

    m = size(aatt,2);
    assert(m>=3);

    if(m==3)
    
        % input = {aatt(:,1),aatt(:,2),aatt(:,3)};
        % R = fun(input{:}); % this is super slow. 
        
        R = fun(aatt(:,1),aatt(:,2),aatt(:,3));
        
        
    else
        if(m==4)
        
            input = {aatt(:,1),aatt(:,2),aatt(:,3),aatt(:,4)};
        else
            assert(m==20);
            
            input = { aatt(:,1), aatt(:,2), aatt(:,3), aatt(:,4), aatt(:,5), aatt(:,6), aatt(:,7), aatt(:,8), aatt(:,9),aatt(:,10), ...
                     aatt(:,11),aatt(:,12),aatt(:,13),aatt(:,14),aatt(:,15),aatt(:,16),aatt(:,17),aatt(:,18),aatt(:,19),aatt(:,20)};
            
        end
        R = fun(input{:});
        warning('This is super slow and inefficient. ');
    end
    
    
end

function [rg] = assemble_pdapdt(s_1, s_2, s_3, aatt, gg, Area)

    f = size(Area, 1);
    m = size(aatt, 2);
    
    assert(size(aatt,1)==f);
    
    rg = zeros(f*m, 1);
    
    for ii=1:m
        rg(1+(ii-1)*f:ii*f) = ...
            Area.*( s_1{ii}(aatt).*gg(1:f) + s_2{ii}(aatt).*gg(f+1:2*f) + s_3{ii}(aatt).*gg(2*f+1:3*f) );
    end

end

function [rg] = assemble_pdrpdt(r, aatt, Area)

    f = size(Area, 1);
    m = size(aatt, 2);
    
    assert(size(aatt,1)==f);
    
    rg = zeros(f*m, 1);
    
    for ii=1:m
        rg(1+(ii-1)*f:ii*f) = Area .* r{ii}(aatt);
    end

end