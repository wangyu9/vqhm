classdef EnergyClassSymDirichlet < handle
    %Symmetric Dirichlet Energy gradient and quadratic proxies computations
    %class
    properties
        eps=1e-6 % eps for pushing PSD proxy to PD
    end
    properties(SetAccess = private)        
        Ai %elements areas
        AA %diagonal matrix of elements areas repeated 4 times

        D1 %derivative operator in frame direction 1 
        D2 %derivative operator in frame direction 2
        DD %nabla operator 

        x %u,v column stack
        
        fi %energy of terms of SOS objective        
        J
        dS
        ds

        %singular values of J (of the mapping)
        S
        s
        u1
        un
        v1
        vn
        
        %dirichlet proxy (constant w.r.t x)
        Hdirichlet
        
        % wangyu added
        
        V
        F
        
        BFGS_iter = 0
        BFGS_B
        BFGS_Bs
        BFGS_Gs
        BFGS_x_old
        
        H_rest = []
        
        X_his = [];
        F_his = [];
        G_his = [];
        n_his = 10;
    end
    
    methods
       function obj=EnergyClassSymDirichlet()
           
       end
       
       function obj=Init(obj,D1,D2,Ai,V,F)
           %D1,D2 derivative operators in frame directions
           %Ai elements areas
           %vertices and faces: V is Nvx3, F is Nfx3
            obj.D1=D1;
            obj.D2=D2;
            obj.Ai=Ai;
            obj.DD=kron(eye(2),[obj.D1;obj.D2]);
            obj.AA=diag(sparse(repmat(obj.Ai,4,1)));
            obj.Hdirichlet = obj.DD'*obj.AA*obj.DD;           
            NewtonMex('Init',V,F);%used in Energy Class
            
            obj.V = V;
            obj.F = F;
        end
        
        function [f,g] = ComputeEnergy(obj)
            f=obj.Ai'*obj.ComputeEnergyPerElement;
            if  nargout>1
                g = (obj.Ai'*obj.ComputeGradientPerElement)';
            end
        end
        
        %method used for calculated scale for min sym-dirichlet energy
        function [fDirichlet,finDirichlet]=ComputeSymmetricDirichletParts(obj)
            fTerms=obj.ComputeEnergyTerms;
            sqrf = fTerms.^2;
            fDirichlet =    0.5*obj.Ai'*(sqrf(:,1) + sqrf(:,3));
            finDirichlet = 0.5*obj.Ai'*(sqrf(:,2) + sqrf(:,4));
        end      
        
        function UpdatePose(obj,x)
            %x is the current solution in the format [u,v] or [u ; v] (2 columns of column stack of parameterization)
            obj.x=x(:);
            obj.SetSingularValues;
        end
        
        %% Proxies
        function H=ProxySLIM(obj)
            Umat=[obj.u1,obj.un];
            Umat = transposeVec2x2(Umat);
            W=zeros(size(Umat,1),4);
            W(:,4) = sqrt(0.5*(obj.s-obj.s.^-3)./(obj.s-1));
            W(:,1) = sqrt(0.5*(obj.S-obj.S.^-3)./(obj.S-1));
            W(abs(obj.S-1)<1e-8,1)=1;
            W(abs(obj.s-1)<1e-8,4)=1;
            WW=sparse(multiplyVec2x2(multiplyVec2x2(Umat,W),transposeVec2x2(Umat)));
            
            WWWW=[diag([WW(:,1);WW(:,1)]),diag([WW(:,2);WW(:,2)]);diag([WW(:,3);WW(:,3)]),diag([WW(:,4);WW(:,4)])];
            
            WWWDD=WWWW*obj.DD;
            H=WWWDD'*obj.AA*WWWDD;
        end
        
        function H=ProxyCompMajor(obj)
            Js = [obj.dS;obj.ds];
            Hs = [(1+3*obj.S.^-4).*obj.Ai; (1+3*obj.s.^-4).*obj.Ai];
            Hggn = Js'*diag(sparse(Hs))*Js; %generalized gauss newton
            gS = obj.Ai.*(obj.S-obj.S.^-3);
            gs = obj.Ai.*(obj.s-obj.s.^-3);
            walpha = gS+gs;
            wbeta = gS-gs;
            walpha(walpha<0) = 0;

            a1=0.5*[obj.D1, obj.D2]; a2=0.5*[-obj.D2,obj.D1]; %similarity cone coefficients 
            b1=0.5*[obj.D1,-obj.D2]; b2=0.5*[ obj.D2,obj.D1]; %antisimilarity cone coefficients
            ha = ComputeConeHessian(a1,a2,obj.x,walpha);
            hb = ComputeConeHessian(b1,b2,obj.x,wbeta);
            
            H=Hggn+ha+hb;
%             bNoProjection = 1;
%             H2 = NewtonMex('Compute',obj.x,bNoProjection); 
%             norm(H(:)-H2(:));
        end
        
        % This is a reference implementation which is less vectorized but
        % might be easier to read
        function H=ProxyCompMajorLoop(obj)
            
            D1=obj.D1/2;
            D2=obj.D2/2;
            a1=[D1, D2]; a2=[-D2,D1];
            b1=[D1,-D2]; b2=[ D2,D1];
            ha = ComputeConeHessianLoop(a1,a2,obj.x);
            hb = ComputeConeHessianLoop(b1,b2,obj.x);
            
            s=obj.Energy.s; S=obj.Energy.S;
            [dS,ds]=obj.Energy.ComputeSingularDerivatives;
            Js = [dS;ds];
            Hs = [(1+3*S.^-4).*obj.Ai; (1+3*s.^-4).*obj.Ai];
            Hggn = Js'*diag(Hs)*Js;%generalized gauss newton
            gS = obj.Ai.*(S-S.^-3);
            gs = obj.Ai.*(s-s.^-3);

            walpha = gS+gs;
            wbeta  = gS-gs;
            walpha(walpha<0) = 0;

            Hab=zeros(size(Hggn));
            for i=1:length(ha)
                Hab=Hab+walpha(i)*ha{i};
                Hab=Hab+wbeta(i)*hb{i};
            end
            H=Hggn+Hab;
        end
        
        function H=ProxyNewton(obj)
            % bNoProjection - flag mitigating performing per face PSD projection
            %             bNoProjection = 0; % means not simple Newton and perform face projection.
            %             default is 0 (perform face projection)
            bNoProjection = 1;
            H = NewtonMex('Compute',obj.x,bNoProjection);            
        end
        
        function H=ProxyProjectedNewton(obj)
            H = NewtonMex('Compute',obj.x);
            %             For Hessians for each face, do this:
            %             [H,HH] = NewtonMex('Compute',x);
        end        
        
        
        function H=ProxyTemp(obj)
            % H = obj.ProxyIntrinsic;
            H = obj.ProxyQuasi;
        end
        
        function H=ProxyIntrinsic(obj)
            
            H = obj.ComputeEdgeHessians(true);
            
        end
        
        function gI=ComputeIntrinsicGradients(obj)
            n = size(obj.V,1);
            VD = reshape(obj.x,[n,2]);
            [~,g,h] = intrinsic_grad_hessian("symmetric-Dirichlet");
            gI = IntrinsicHessianClass.ComputeIntrinsicGradients(VD, obj.F, obj.V, g);
        end
        
        function HI=ComputeIntrinsicHessians(obj)
            n = size(obj.V,1);
            VD = reshape(obj.x,[n,2]);
            [~,g,h] = intrinsic_grad_hessian("symmetric-Dirichlet");
            HI = IntrinsicHessianClass.ComputeIntrinsicHessians(VD, obj.F, obj.V, h);
        end
        
        function H=HessiansIntrinsic2Extrinsic(obj, HI, gI, assemble, project)
            
            n = size(obj.V,1);
            VD = reshape(obj.x,[n,2]);
            H = IntrinsicHessianClass.HessiansIntrinsic2Extrinsic(VD, obj.F, obj.V, obj.Ai, HI, gI, assemble, project);
        end
        
        function gI=ComputeIntrinsicGradients2(obj)
            % 3 x f
            
            %
            
            E0 = edge_lengths(obj.V,obj.F);
            
            n = size(obj.V,1);
            f = size(obj.F,1);
            
            VD = reshape(obj.x,[n,2]);

            E = edge_lengths(VD,obj.F);
            
            [~,g,h] = intrinsic_grad_hessian("symmetric-Dirichlet");
            
            % [He1,He2,He3] = extrinsic_edge_hessians();
            
%             IDX = [...
%                 obj.F(:,1),n+obj.F(:,1),...
%                 obj.F(:,2),n+obj.F(:,2),...
%                 obj.F(:,3),n+obj.F(:,3),...
%             ];

            %
            
            gI = zeros(3,f);
            
            for i=1:f
            
                a0 = E0(i,1).^2;
                a1 = E(i,1).^2;
                b0 = E0(i,2).^2;
                b1 = E(i,2).^2; 
                c0 = E0(i,3).^2; 
                c1 = E(i,3).^2;

                gi = g(a0,a1,b0,b1,c0,c1);
                % hi = h(a0,a1,b0,b1,c0,c1);

                idx = [...
                    obj.F(i,1),n+obj.F(i,1),...
                    obj.F(i,2),n+obj.F(i,2),...
                    obj.F(i,3),n+obj.F(i,3),...
                    ];

                VTi = reshape(obj.x(idx),[2,3])';
                % KT = grad_edgesq_over_vertex(VTi);
                %KT: deodv, 6 x 3

                gI(:,i) = gi;
                
            end
        end
        
        function HI=ComputeIntrinsicHessians2(obj)
            % 3 x 3 x f
            
            %
            
            E0 = edge_lengths(obj.V,obj.F);
            
            n = size(obj.V,1);
            f = size(obj.F,1);
            
            VD = reshape(obj.x,[n,2]);

            E = edge_lengths(VD,obj.F);
            
            [~,g,h] = intrinsic_grad_hessian("symmetric-Dirichlet");
            
            % [He1,He2,He3] = extrinsic_edge_hessians();
            
%             IDX = [...
%                 obj.F(:,1),n+obj.F(:,1),...
%                 obj.F(:,2),n+obj.F(:,2),...
%                 obj.F(:,3),n+obj.F(:,3),...
%             ];

            %
            
            HI = zeros(3,3,f);
            
            for i=1:f
            
                a0 = E0(i,1).^2;
                a1 = E(i,1).^2;
                b0 = E0(i,2).^2;
                b1 = E(i,2).^2; 
                c0 = E0(i,3).^2; 
                c1 = E(i,3).^2;

                % gi = g(a0,a1,b0,b1,c0,c1);
                hi = h(a0,a1,b0,b1,c0,c1);

                idx = [...
                    obj.F(i,1),n+obj.F(i,1),...
                    obj.F(i,2),n+obj.F(i,2),...
                    obj.F(i,3),n+obj.F(i,3),...
                    ];

                VTi = reshape(obj.x(idx),[2,3])';
                % KT = grad_edgesq_over_vertex(VTi);
                %KT: deodv, 6 x 3

                HI(:,:,i) = hi;
                
            end
        end
        
        function H=HessiansIntrinsic2Extrinsic2(obj, HI, gI, assemble, project)
            % 6 x 6 x f from 3 x 3 x f
                        
            %
            
            E0 = edge_lengths(obj.V,obj.F);
            
            n = size(obj.V,1);
            f = size(obj.F,1);
            
            VD = reshape(obj.x,[n,2]);

            % R = []
            E = edge_lengths(VD,obj.F);
            
            [~,g,h] = intrinsic_grad_hessian("symmetric-Dirichlet");
            
            [He1,He2,He3] = extrinsic_edge_hessians();
            
            IDX = [...
                obj.F(:,1),n+obj.F(:,1),...
                obj.F(:,2),n+obj.F(:,2),...
                obj.F(:,3),n+obj.F(:,3),...
            ];
            %
            
            H_ele = zeros(6,6,f);
            
            for i=1:f
            
%                 a0 = E0(i,1).^2;
%                 a1 = E(i,1).^2;
%                 b0 = E0(i,2).^2;
%                 b1 = E(i,2).^2; 
%                 c0 = E0(i,3).^2; 
%                 c1 = E(i,3).^2;
% 
%                 gi = g(a0,a1,b0,b1,c0,c1);
%                 hi = h(a0,a1,b0,b1,c0,c1);

                hi = HI(:,:,i);
                gi = gI(:,i);

                idx = [...
                    obj.F(i,1),n+obj.F(i,1),...
                    obj.F(i,2),n+obj.F(i,2),...
                    obj.F(i,3),n+obj.F(i,3),...
                    ];

                VTi = reshape(obj.x(idx),[2,3])';
                KT = grad_edgesq_over_vertex(VTi);
                %KT: deodv, 6 x 3

                switch project
                    case 0
                        Hi = ( KT * hi *KT' + He1 * gi(1) + He2 * gi(2) + He3 * gi(3) ) * obj.Ai(i);
                    case 1
                        Hi = ( KT * hi *KT' + He1 * gi(1) + He2 * gi(2) + He3 * gi(3) ) * obj.Ai(i);
                        [EV0,ED0] = eig((Hi+Hi')/2); 
                        Hi = EV0*max(ED0,0)*EV0';
                    case 2
                        Hi = ( KT * hi *KT' + He1 * max(gi(1),0) + He2 * max(gi(2),0) + He3 * max(gi(3),0) ) * obj.Ai(i);
                    case 3
                        EHi = He1 * gi(1) + He2 * gi(2) + He3 * gi(3);
                        [EV,ED] = eig(EHi); 
                        % [EV2,ED2] = eig(hi); 
                        % Hi = ( KT * EV2*max(ED2,0)*EV2' *KT' + EV*max(ED,0)*EV' ) * obj.Ai(i);
                        Hi = ( KT * hi *KT' + EV*max(ED,0)*EV' ) * obj.Ai(i);
                    otherwise
                        error(['unsupport project']);
                end
                                
                H_ele(:,:,i) = Hi;
            
            end
            
            if assemble
                H = sparse([],[],[],n*2,n*2);
                        
                for i=1:6
                   for j=1:6
                      EH = squeeze(H_ele(i,j,:));
                      H = H + sparse(IDX(:,i),IDX(:,j),EH,n*2,n*2);
                   end
                end 
            else
                H = H_ele;
            end
        end
        
        
        function H=ComputeEdgeHessians(obj, assemble)
            % 6 x 6 x f
            
            [~,Hs] = NewtonMex('Compute',obj.x);
            ge = obj.ComputeGradientPerElement;
            
            %
            
            E0 = edge_lengths(obj.V,obj.F);
            
            n = size(obj.V,1);
            f = size(obj.F,1);
            
            VD = reshape(obj.x,[n,2]);

            % R = []
            E = edge_lengths(VD,obj.F);
            
            [~,g,h] = intrinsic_grad_hessian("symmetric-Dirichlet");
            
            [He1,He2,He3] = extrinsic_edge_hessians();
            
            IDX = [...
                obj.F(:,1),n+obj.F(:,1),...
                obj.F(:,2),n+obj.F(:,2),...
                obj.F(:,3),n+obj.F(:,3),...
            ];
            %
            
            H_ele = zeros(6,6,f);
            
            for i=1:f
            
                a0 = E0(i,1).^2;
                a1 = E(i,1).^2;
                b0 = E0(i,2).^2;
                b1 = E(i,2).^2; 
                c0 = E0(i,3).^2; 
                c1 = E(i,3).^2;

                gi = g(a0,a1,b0,b1,c0,c1);
                hi = h(a0,a1,b0,b1,c0,c1);

                idx = [...
                    obj.F(i,1),n+obj.F(i,1),...
                    obj.F(i,2),n+obj.F(i,2),...
                    obj.F(i,3),n+obj.F(i,3),...
                    ];

                VTi = reshape(obj.x(idx),[2,3])';
                KT = grad_edgesq_over_vertex(VTi);
                %KT: deodv, 6 x 3

                Hi = ( KT * hi *KT' + He1 * gi(1) + He2 * gi(2) + He3 * gi(3) ) * obj.Ai(i);
                
%                 % verified it agrees with (Hs(:,:,i))
%                 if(norm((Hs(:,:,i))-Hi)>1e-8*norm(Hs(:,:,i)))
%                    assert(false); 
%                 end
                
                
                Hi = ( KT * hi *KT' + He1 * max(gi(1),0) + He2 * max(gi(2),0) + He3 * max(gi(3),0) ) * obj.Ai(i);

                
                if false
                    EHi = He1 * gi(1) + He2 * gi(2) + He3 * gi(3);
                    [EV,ED] = eig(EHi); 
                    % [EV2,ED2] = eig(hi); 
                    % Hi = ( KT * EV2*max(ED2,0)*EV2' *KT' + EV*max(ED,0)*EV' ) * obj.Ai(i);
                    Hi = ( KT * hi *KT' + EV*max(ED,0)*EV' ) * obj.Ai(i);
                end
                
                if false
                   Hi = ( KT * hi *KT' + He1 * gi(1) + He2 * gi(2) + He3 * gi(3) ) * obj.Ai(i);
                   [EV0,ED0] = eig((Hi+Hi')/2); 
                   Hi = EV0*max(ED0,0)*EV0';
                end
                
                % gei = full(ge(i,idx))'; 
                % verified it agrees with KT * gi
                
                H_ele(:,:,i) = Hi;
            
            end
            
            if assemble
                H = sparse([],[],[],n*2,n*2);
                        
                for i=1:6
                   for j=1:6
                      EH = squeeze(H_ele(i,j,:));
                      H = H + sparse(IDX(:,i),IDX(:,j),EH,n*2,n*2);
                   end
                end 
            else
                H = H_ele;
            end
        end
        
        function H=ProxyQuasi(obj)
            
            %             For Hessians for each face, do this:
            %             [H,HH] = NewtonMex('Compute',x);
            
            %%            
            E0 = edge_lengths(obj.V,obj.F);
            
            n = size(obj.V,1);
            f = size(obj.F,1);
            
            VD = reshape(obj.x,[n,2]);

            % R = []
            E = edge_lengths(VD,obj.F);
            
            [~,g,h] = intrinsic_grad_hessian("symmetric-Dirichlet");
            
            
            [He1,He2,He3] = extrinsic_edge_hessians();

            IDX = [...
                obj.F(:,1),n+obj.F(:,1),...
                obj.F(:,2),n+obj.F(:,2),...
                obj.F(:,3),n+obj.F(:,3),...
            ];
            %%
            
            
            % restart_cond = (mod(obj.BFGS_iter,5)==0);
            % restart_cond = obj.BFGS_iter < 10;
            % restart_cond = obj.BFGS_iter==0;
            % restart_cond = obj.BFGS_iter==0;
            restart_cond = true;
            % restart_cond = obj.BFGS_iter <= 12;
            
            edge_hess = true;
            p = 0;
            if edge_hess
                p = 3;
            else
                p = 6;
            end
            
            if(isempty(obj.F_his))
                obj.X_his = zeros(obj.n_his,p,f);
                obj.F_his = zeros(obj.n_his,f);
                obj.G_his = zeros(obj.n_his,p,f);
            end
            
            Xs = zeros(p,f);
            Fs = zeros(f);
            
            if restart_cond
                
                %  obj.Hdirichlet;
                
                obj.BFGS_Gs = zeros(p,f);

                % note the local hessians are always unprojected.

                % as in the c++ code, the 6 by 6 matrix is indexed as
                % x1,y1,x2,y2,x3,y3, indices 1,2,3 corresponds to F(:,1),F(:,2),F(:,3)
                    
                if true
                    
                    if edge_hess
                        obj.BFGS_Bs = obj.ComputeIntrinsicHessians();
                        % not area weighted
                    else
                        obj.BFGS_B = NewtonMex('Compute',obj.x);
                        [~,obj.BFGS_Bs] = NewtonMex('Compute',obj.x);
                        % is area weighted
                    end
                
                else
                    
                    obj.BFGS_Bs = zeros(p,p,f);
                    
                    for i=1:size(obj.BFGS_Bs,3)
                        obj.BFGS_Bs(:,:,i) = eye(p) * obj.Ai(i);
                    end
                                        
                end
                
            end
            
            
            if true
                
                if edge_hess                    
                    gI = obj.ComputeIntrinsicGradients;
                    EL = edge_length_per_tri(reshape(obj.x,[n,2]), obj.F);
                    
                else
                    ge = obj.ComputeGradientPerElement;
                    % ge has not been weighted with per-element area yet
  
                end
                
                Fs = obj.Ai .* obj.ComputeEnergyPerElement; % the energies are unweighted with area. 
                
                
                for i=1:f
                    
                    idx = [...
                        obj.F(i,1),n+obj.F(i,1),...
                        obj.F(i,2),n+obj.F(i,2),...
                        obj.F(i,3),n+obj.F(i,3),...
                        ];
                    
                    if edge_hess
                        VTi = reshape(obj.x(idx),[2,3])';
                        % JT = grad_vertex_over_edge(VTi);
                        % JT = grad_vertex_over_edgesq(VTi);
                        % gi = JT * gi;
                        gi = gI(:,i);
                        
                        % not area weighted
                    else
                        % gi = full([
%                         ge(i,obj.F(i,1)), ge(i,n+obj.F(i,1)),...
%                         ge(i,obj.F(i,2)), ge(i,n+obj.F(i,2)),...
%                         ge(i,obj.F(i,3)), ge(i,n+obj.F(i,3))]');

                        gi = full(ge(i,idx))';
                        
                        gi = gi * obj.Ai(i);
                        
                        % is area weighted
                    end

                    gi_old = obj.BFGS_Gs(:,i);
                    
                    obj.BFGS_Gs(:,i) = gi;
                    
                    yk = gi - gi_old;

                    if edge_hess
                        Xs(:,i) = EL(i,:);
                    else
                        Xs(:,i) = obj.x(idx);
                    end
                        
                    if false
                    %if ~restart_cond   
                    
                        if edge_hess
                            sk = (EL(i,:) - obj.BFGS_x_old(i,:))';
                        else
                            sk = obj.x(idx) - obj.BFGS_x_old(idx);
                        end
                        
                        Bk = obj.BFGS_Bs(:,:,i);

                        formula = 5;
                        
                        switch formula
                            
                            case 0
                            % BFGS update% note it is not psd for nonconvex
                            % Bkp = Bk + yk*yk'/(yk'*sk) - Bk*sk*(Bk*sk)'/(sk'*Bk*sk);
                            
                                Bkp = Bk;
                                if abs(yk'*sk)>0
                                    Bkp = Bkp + yk*yk'/(yk'*sk);
                                end
                                if abs(sk'*Bk*sk)>0
                                    Bkp = Bkp - Bk*sk*(Bk*sk)'/(sk'*Bk*sk);
                                end
                            
                            case 1
                            % ad-hoc
                            Bkp = Bk + ( yk*yk'/(yk'*sk) - Bk*sk*(Bk*sk)'/(sk'*Bk*sk) ) * ((yk'*sk)>0);

                            case 2
                            % ad-hoc2
                            Bkp = Bk + yk*yk'/abs(yk'*sk) - Bk*sk*(Bk*sk)'/(sk'*Bk*sk);
                            
                            case 5

                            % symmetric-rank-one (SR1) formula
                        
                            Bkp = Bk; 
                            if (((yk-Bk*sk)'*sk)>0)
                                Bkp = Bk + (yk-Bk*sk)*(yk-Bk*sk)' / ((yk-Bk*sk)'*sk);
                            end
                        
                            otherwise
                                error(['undefined formula!']);
                        end
                        
                        if(any(isnan(Bkp(:))))
                            error('Bkp has nan');
                        end

                        obj.BFGS_Bs(:,:,i) = Bkp;
                    else
                        
                        if true
                        [VVV,DDD] = eig((obj.BFGS_Bs(:,:,i)+obj.BFGS_Bs(:,:,i)')/2);
                        e = diag(DDD);
                       	% e(e < 0)= 0;
                        e = abs(e);
                        obj.BFGS_Bs(:,:,i) = VVV*diag(e)*VVV';
                        end
                    end
                    
                end
            end
            
            nh = size(obj.G_his,1);
            obj.X_his(nh+1,:,:) = Xs;
            obj.G_his(nh+1,:,:) = obj.BFGS_Gs;
            obj.F_his(nh+1,:) = Fs;
            
            obj.BFGS_iter = obj.BFGS_iter + 1;
            
            if edge_hess
                obj.BFGS_x_old = EL;
            else
                obj.BFGS_x_old = obj.x;
            end
            
            H = sparse([],[],[],n*2,n*2);
            
            if edge_hess
                
%                 H_ele = zeros(6,6,f);
%                 
%                 for i=1:f
%                     idx = [...
%                         obj.F(i,1),n+obj.F(i,1),...
%                         obj.F(i,2),n+obj.F(i,2),...
%                         obj.F(i,3),n+obj.F(i,3),...
%                         ];
%                     
%                     VX = reshape(obj.x(idx),[2,3])';
%                     KT = grad_edge_over_vertex(VX);
%                     %KT: deodv, 6 x 3
%                     
%                     H_ele(:,:,i) = KT * obj.BFGS_Bs(:,:,i) * KT' + 1e-4 * eye(6);
%                 end
                
                project = 1; % 0: almost no progress.
                H_ele = obj.HessiansIntrinsic2Extrinsic(obj.BFGS_Bs, gI, false, project);

            else
                H_ele = obj.BFGS_Bs;
            end
            
            for i=1:6
               for j=1:6
                  EH = squeeze(H_ele(i,j,:));
                  H = H + sparse(IDX(:,i),IDX(:,j),EH,n*2,n*2);
               end
            end
            
            % checked that at iter 1 
            % sum(sum(abs(H - obj.BFGS_B))) is ~0
            % H = obj.BFGS_B; % should remove this line 
            
            % H = obj.BFGS_B * 0.0 + H * 1;
            
        end       
        
        function H=ProxyProjectedNewtonFull(obj)
           H = obj.ProxyNewton;
           [V,D] = eig(full(0.5*(H+H')));
           e = diag(D);
           e(e < obj.eps)= obj.eps;
           H = V*diag(e)*V';
        end        
    end
    
    methods(Access = private)
        function fe=ComputeEnergyPerElement(obj)
            obj.fi=obj.ComputeEnergyTerms;
            fe=0.5*sum(obj.fi.^2,2);
        end
        
        function f=ComputeEnergyTerms(obj)
            f=[obj.S,1./obj.S,obj.s,1./obj.s];
        end
        
        function ge=ComputeGradientPerElement(obj)
            obj.ComputeSingularDerivatives;
            Nf=numel(obj.S);
            
            dSE=[ones(Nf,1),-1./(obj.S.^2)]; dsE=[ones(Nf,1),-1./(obj.s.^2)];
            
            gS=(dSE(:,1).*obj.fi(:,1)+dSE(:,2).*obj.fi(:,2));
            gs=(dsE(:,1).*obj.fi(:,3)+dsE(:,2).*obj.fi(:,4));
            
            ge = diag(sparse(gS))*obj.dS+...
                diag(sparse(gs))*obj.ds;
        end
 
        function SetSingularValues(obj)
            Nv = length(obj.x)/2;
            U = obj.x(1:Nv);
            V = obj.x(Nv+1:end);
            
            % enetries of Jacobian of all elements (triangles) 
            a=obj.D1*U;
            b=obj.D2*U;
            c=obj.D1*V;
            d=obj.D2*V;
            
            %save for checking flips of elements
            obj.J=a.*d-b.*c;
            
            D1P=obj.D1*[U,V]/2; D2P=obj.D2*[U,V]/2;
            E=D1P(:,1)+D2P(:,2); F=D1P(:,1)-D2P(:,2);
            G=D2P(:,1)+D1P(:,2); H=-D2P(:,1)+D1P(:,2);
            Q=sqrt(E.^2+H.^2); R=sqrt(F.^2+G.^2);
            
            a1=atan2(G,F); a2=atan2(H,E);           
            theta=(a2-a1)/2; phi=(a2+a1)/2;
            
            %save singular values and vectors
            obj.s=Q-R;
            obj.S=Q+R;
            obj.u1=[  cos(phi)  sin(phi)];
            obj.un=[ -sin(phi)  cos(phi)];
            obj.v1=[ cos(theta) -sin(theta)];
            obj.vn=[ sin(theta)  cos(theta)];  
        end
        
        function ComputeSingularDerivatives(obj)
            b=diag(sparse(obj.v1(:,1)))*obj.D1+diag(sparse(obj.v1(:,2)))*obj.D2;
            c=diag(sparse(obj.vn(:,1)))*obj.D1+diag(sparse(obj.vn(:,2)))*obj.D2;
            obj.dS=[diag(sparse(obj.u1(:,1)))*b,diag(sparse(obj.u1(:,2)))*b];
            obj.ds=[diag(sparse(obj.un(:,1)))*c,diag(sparse(obj.un(:,2)))*c];
        end
    end
end

