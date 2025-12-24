function [tp] = tensor_para(FA,dim,para_type)

tp = struct();

tp.reg_grad = @(at) 0;
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


    end


tp.s_at2au = s_at2au;
tp.s_pdapdt_lmul = s_pdapdt_lmul;
            
tp.s_at2au_v = s_at2au_v;
tp.s_pdapdt_lmul_v = s_pdapdt_lmul_v;

